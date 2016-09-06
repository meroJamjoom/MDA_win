// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: atcheson (Bradley Atcheson)
// Email:   atcheson@cs.ubc.ca
// ==========================================================================

#ifndef BRIDGE_OPENCVBRIDGE_C
#define BRIDGE_OPENCVBRIDGE_C

#include "OpenCVBridge.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;



/** default constructor */
OpenCVReadBridge::OpenCVReadBridge(): 
  width(-1), height(-1), depth(-1), channels(-1)
{
}


bool
OpenCVReadBridge::connect( istream& is )
{
  reader.connect( is );
  if( validateHeader() ) {
    configure();
    return true;
  } else {
    return false;
  }
}


bool
OpenCVReadBridge::connect( const char* filename )
{
  reader.connect( filename );
  if( validateHeader() ) {
    configure();
    return true;
  } else {
    return false;
  }
}


void
OpenCVReadBridge::disconnect()
{
  reader.disconnect();
  width = -1;
  height = -1;
  depth = -1;
  channels = -1;
}


bool
OpenCVReadBridge::validateHeader()
{
  if( !reader.readHeader() )
  {
    cerr << "Cannot read MDA header" << endl;
    return false;
  }
  if( reader.getDim().vec.size() < 2 ) {
    cerr << "Cannot convert images with fewer than 2 dimensions" << endl;
    return false;
  }
  if( (reader.getNumChannels() < 1) || (reader.getNumChannels() > 4) ) {
    cerr << "Can only convert arrays with 1-4 channels" << endl;
    return false;
  }

  return true;
}


unsigned
OpenCVReadBridge::getNumFrames()
{
  // if image is of size n0xn1xn2x...xnN then numFrames = n2*n3*...*nN
  unsigned numFrames = 1;
  for( unsigned i=2; i<reader.getDim().vec.size(); i++ ) {
    numFrames *= reader.getDim().vec[i];
  }
  return numFrames;
}


bool
OpenCVReadBridge::moreFrames()
{
  return reader.getNumScanlinesLeft() > 0;
}


void
OpenCVReadBridge::configure()
{
  // get frame size
  width = reader.getDim().vec[0];
  height = reader.getDim().vec[1];
  channels = reader.getNumChannels();

  // map MDA to OpenCV types
  targetType = reader.getType();
  switch( reader.getType() ) {
    case UByte:
      depth = IPL_DEPTH_8U;
      break;
    case Byte:
      depth = IPL_DEPTH_8S;
      break;
    case Short:
      depth = IPL_DEPTH_16S;
      break;
    case UShort:
      depth = IPL_DEPTH_16U;
      break;
    case Int:
      depth = IPL_DEPTH_32S;
      break;
    case UInt:
      depth = IPL_DEPTH_32S;
      targetType = Int;
      cerr << "Warning: OpenCV does not support 32bit unsigned int. "
           << "Data will be converted to 32bit signed int." << endl;
      break;
    case Float:
      depth = IPL_DEPTH_32F;
      break;
    case Double:
      depth = IPL_DEPTH_64F;
      break;
  }
}


IplImage*
OpenCVReadBridge::createEmptyFrame()
{
  // create the opencv image buffer
  CvSize size = cvSize( width, height );
  IplImage* img = cvCreateImage( size, depth, reader.getNumChannels() );
  if( !img ) {
    cerr << "Could not allocate OpenCV image buffer" << endl;
    return NULL;
  } else {
    return img;
  }
}

IplImage*
OpenCVReadBridge::extractNextFrame( IplImage* img )
{
  // allocate new image if caller didn't already do so
  if( !img ) {
    img = createEmptyFrame();
  } else {
    // check that the caller's img is compatible
    if( (img->width!=width) || (img->height!=height) || (img->depth!=depth) ) {
      cerr << "OpenCV image has incompatible size/type" << endl;
      return NULL;
    }
  } 

  // It would probably be faster to do some kind of direct memory copy (or,
  // ideally to just set the OpenCV buffer's data pointer to the start of
  // the raw MDA data) but the OpenCV internal format is more complex than
  // MDA (scanlines are padded to align nicely in memory, and the channels are
  // stored separately in memory rather then interleaved). So it's only really
  // a very small class of MDA files that could be directly loaded by OpenCV.
  // It's just neater and safer to always run through the typeConvert function
  // although that makes it perhaps unsuitable for real-time stuff.

  // prepare pointers
  unsigned long scanlineSize = reader.getScanlineSize();
  char* scanline = new char[scanlineSize];
  unsigned long numEntries = width * channels;

  // copy and convert each scanline...
  for( unsigned row=0; row<height; row++ ) {
    // copy and convert scanline to temp buffer
    char* scanlineIn = (char*)reader.readScanline();
    typeConvert( scanlineIn,reader.getType(),scanline,targetType,numEntries );
    // copy scanline into output image buffer
    memcpy( (char*)img->imageData+row*img->widthStep, scanline, scanlineSize );
  }

  delete[] scanline;
  return img;
}


void
OpenCVReadBridge::skipNextFrame()
{
  if( reader.getNumScanlinesLeft() > height ) {
    for( unsigned i=0; i<height; i++ ) {
      reader.skipScanline();
    }
  } else {
    cerr << "Cannot skip past end of stream" << endl;
  }
}



OpenCVWriteBridge::OpenCVWriteBridge( IplImage* img, unsigned nFrames )
{
  // set the output MDA size. This will always be a stack of frames, where
  // each frame is the same 2D size as the OpenCV array. The total number 
  // of frames needs to be known ahead of time by the caller.
  dim.vec.push_back( img->width );
  dim.vec.push_back( img->height );
  if( nFrames > 1 ) {
    dim.vec.push_back( nFrames );
  }

  // set the channel depth
  channels = img->nChannels;

  // map the OpenCV types to MDA types
  depth = img->depth;
  switch( img->depth ) {
    case IPL_DEPTH_8U:
      type = UByte;
      break;
    case IPL_DEPTH_8S:
      type = Byte;
      break;
    case IPL_DEPTH_16U:
      type = UShort;
      break;
    case IPL_DEPTH_16S:
      type = Short;
      break;
    case IPL_DEPTH_32S:
      type = Int;
      break;
    case IPL_DEPTH_32F:
      type = Float;
      break;
    case IPL_DEPTH_64F:
      type = Double;
      break;
  }
}


bool
OpenCVWriteBridge::connect( ostream& os )
{
  writer.connect( os );
  if( writer.writeHeader(dim,channels,type) ) {
    return true;
  } else {
    cerr << "Cannot write MDA header" << endl;
    return false;
  }
}


bool
OpenCVWriteBridge::connect( const char* filename )
{
  writer.connect( filename );
  if( writer.writeHeader(dim,channels,type) ) {
    return true;
  } else {
    cerr << "Cannot write MDA header" << endl;
    return false;
  }
}


void
OpenCVWriteBridge::disconnect()
{
  writer.disconnect();
}


bool
OpenCVWriteBridge::writeFrame( const IplImage* img )
{
  // double check that the types are compatible
  if( (img->width!=dim.vec[0]) ||
      (img->height!=dim.vec[1]) || 
      (img->depth!=depth) ) {
    cerr << "OpenCV image has incompatible size/type" << endl;
    return false;
  }

  // copy each row
  CvMat row;
  for( unsigned i=0; i<img->height; i++ ) {
    cvGetRow( img, &row, i );
    writer.writeScanline( row.data.ptr );
  }

  return true;
}


} /* namespace */

#endif /* BRIDGE_OPENCVBRIDGE_C */

