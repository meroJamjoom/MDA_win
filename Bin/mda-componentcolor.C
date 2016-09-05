// ==========================================================================
// $Id: mda-componentcolor.C 260 2008-12-30 09:21:39Z heidrich $
// convert YCbCr to RGB and vice versa
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#include <assert.h>

#include "MDA/Base/Range.hh"
#include "MDA/Base/ChannelList.hh"
#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  unsigned long i, j;
  
  CommandlineParser parser;
  
  // setup options
  DataType	outType= UByte;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  // expose these to commandline options later
  vector<double> params;
  UIntRange lumaRange;
  lumaRange.val.first= 16;
  lumaRange.val.second= 235;
  UIntRange chromaRange;
  chromaRange.val.first= 16;
  chromaRange.val.second= 240;
  
  // forward or inverse conversion
  bool forward= true;
  BoolOption forwardOpt( forward, "\tdirection of the conversion\n",
			 "--YCbCr-to-RGB", NULL, "--RGB-to-YCbCr", NULL );
  parser.registerOption( &forwardOpt );
  
  // video standard
  int standard= 0;
  list<const char *>standards;
  //  standards.push_back( "--custom" );
  standards.push_back( "--bt709" );
  standards.push_back( "--bt601" );
  standards.push_back( "--jpeg" );
  SelectionOption standardOpt( standard,
                               "\tselect parameters according to image/video"
                               " standard:\n"
                               "\tITU-R BT.601 (Digital SD),"
			       "ITU-R BT.709 (HDTV), JPEG\n",
                               standards );
  parser.registerOption( &standardOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], "[options]" );
    exit( 1 );
  }
  
  switch( standard )
  {
  case 0: // BT.709 - High Definition Video (60fps version!)
    params.clear();
    params.push_back( .2126 );
    params.push_back( .7152 );
    params.push_back( .0722 );
    params.push_back( .5389 );
    params.push_back( .6350 );
    lumaRange.val.first= chromaRange.val.first= 16;
    lumaRange.val.second= 235;
    chromaRange.val.second= 240;
    break;
  case 1: // BT.601 - Standard Definition Video (60fps version!)
    params.clear();
    params.push_back( .2990 );
    params.push_back( .5870 );
    params.push_back( .1140 );
    params.push_back( .5643 );
    params.push_back( .7133 );
    lumaRange.val.first= chromaRange.val.first= 16;
    lumaRange.val.second= 235;
    chromaRange.val.second= 240;
    break;
  case 2: // JPEG - same as BT.601 but with full quantization range
    params.clear();
    params.push_back( .2990 );
    params.push_back( .5870 );
    params.push_back( .1140 );
    params.push_back( .5643 );
    params.push_back( .7133 );
    lumaRange.val.first= chromaRange.val.first= 0;
    lumaRange.val.second= chromaRange.val.second= 255;
    break;
  }
  
  if( params.size()!= 5 )
  {
    cerr << argv[0]<< ": Need exactly 5 conversion parameters!\n";
    exit( 1 );
  }

  // setup MDA reader and read input header
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() )
  {
    cerr << "Cannot read file header\n";
    exit( 1 );
  }
  DataType	inType= reader.getType();
  CoordinateVector	dim= reader.getDim();
  unsigned int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long inScanlineSize= reader.getScanlineSize();
  unsigned char *inScanline;
  
  // check the YCbCr data type
  if( (forward && inType!= UByte) || (!forward && outType!= UByte) )
  {
    cerr << "Currently, YCbCr has to be represented as unsigned byte\n";
    exit( 1 );
  }
  
  // check number of channels
  if( numChannels!= 3 )
  {
    cerr << argv[0] << ": Expecting exactly 3 channels!\n";
    exit( 1 );
  }
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannels, outType ) )
  {
    cerr << argv[0] << ": Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  unsigned long outScanlineSize= writer.getScanlineSize();
  char *outScanline= new char[outScanlineSize];
  
  // rgb and yuv triplets for one pixel
  double rgb[3];
  double yuv[3];
  // scale factors for the individual channels
  double lumaScale= 1.0 / (lumaRange.val.second-lumaRange.val.first);
  double chromaScale= 1.0 / (chromaRange.val.second-chromaRange.val.first);
  
  // actually convert the data
  unsigned int bytesPerValue= dataTypeSizes[outType];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    inScanline= (unsigned char *)reader.readScanline();
      
    if( forward )
    {
      // YCbCr to RGB
      for( j= 0 ; j< dim.vec[0] ; j++ )
      {
	// de-quantize YCbCr representation
	yuv[0]= ((double)inScanline[j*3]-lumaRange.val.first) * lumaScale;
	yuv[1]=
	  ((double)inScanline[j*3+1]-chromaRange.val.first) * chromaScale - .5;
	yuv[2]=
	  ((double)inScanline[j*3+2]-chromaRange.val.first) * chromaScale - .5;
	
	// convert to RGB
	rgb[2]= 1.0/params[3] * yuv[1] + yuv[0];
	rgb[0]= 1.0/params[4] * yuv[2] + yuv[0];
	rgb[1]= 1.0/params[1] * (yuv[0] - params[0]*rgb[0] - params[2]*rgb[2]);
	
	// convert triplet from double to whatever is desired
	typeConvert( rgb, Double,
		     outScanline+j*3*dataTypeSizes[outType], outType, 3 );
      }
    }
    else
    {
      // RGB to YCbCr
      for( j= 0 ; j< dim.vec[0] ; j++ )
      {
	// get RGB
	typeConvert( inScanline+j*3*dataTypeSizes[inType], inType,
		     rgb, Double, 3 );
	
	// compute raw YCbCr
	yuv[0]= params[0] * rgb[0] + params[1] * rgb[1] + params[2] * rgb[2];
	yuv[1]= params[3] * (rgb[2] - yuv[0]);
	yuv[2]= params[4] * (rgb[0] - yuv[0]);
	
	// remap range
	outScanline[j*3+0]=
	  (unsigned char)(yuv[0] / lumaScale + lumaRange.val.first);
	outScanline[j*3+1]=
	  (unsigned char)((yuv[1] + 0.5) / chromaScale + chromaRange.val.first);
	outScanline[j*3+2]=
	  (unsigned char)((yuv[2] + 0.5) / chromaScale + chromaRange.val.first);
      }    
    }
    
    // write scanline out
    writer.writeScanline( outScanline );
  }
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << argv[0] << ": Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
