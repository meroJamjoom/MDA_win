// ==========================================================================
// $Id: MDAFileIO.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef MDAFILEIO_C
#define MDAFILEIO_C

#include <string.h>
//#include <stdlib.h>
#if defined (_WIN32) || defined (_WIN64)
#include <io.h>
#include <fcntl.h>
#endif
#include "MDA/Base/Errors.hh"
#include "MDAFileIO.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inlcusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  //
  // MDA_IO methods
  //


/** calculate scanline size and number of scanlines */
void
MDA_IO::setupScanlineInfo()
{
  if( numChannels== 0 || dim.vec.size()== 0 )
  {
    scanlineSize= numScanlinesLeft= 0;
    return;
  }
  
  // the size, in bytes, of one scanline
  scanlineSize= numChannels*dim.vec[0]*dataTypeSizes[dataType];
  
  // number of scanlines is the product of all but the first dimension
  numScanlinesLeft= 1;
  for( int i= 1 ; i< dim.vec.size() ; i++ )
    numScanlinesLeft*= dim.vec[i];
}

  //
  // MDAReader methods
  //

/** default constructor */
MDAReader::MDAReader()
{
#if defined (_WIN32) || defined (_WIN64) 
	_setmode(_fileno(stdin), _O_BINARY);
#endif
  scanline= NULL;
  ownIstream= false;
  inStream= NULL;
  numScanlinesLeft= 0;
}

/** destructor */
MDAReader::~MDAReader()

{
  if( scanline!= NULL )
    delete [] scanline;
  
  disconnect();
} 
    
/** connect to provided istream */
bool
MDAReader::connect( istream &is )
{
  inStream= &is;
  inStream->clear();
  ownIstream= false;
  
  return inStream->good();
}
    
/** connect to open provided file */
bool
MDAReader::connect( const char *fileName )
{
#if defined (_WIN32) || defined (_WIN64)
	myIstream.open(fileName, ifstream::in | ifstream::binary);
#else 
	myIstream.open(fileName, ifstream::in & ifstream::binary);
#endif
 
  ownIstream= myIstream.good();
  if( ownIstream )
    inStream= &myIstream;
  
  return ownIstream;
}


/** disconnect from a stream */
void
MDAReader::disconnect()
{
  numScanlinesLeft= 0;
  if( ownIstream )
  {
    myIstream.close();
    ownIstream= false;
    inStream= NULL;
  }
}

/** read MDA header from the istream we are connected to */
bool
MDAReader::readHeader()
{
  char buf[16];
  
  // skip whitespace and read magic header
  float version;
  *inStream >> ws;
  inStream->width( 16 );
  *inStream >> buf >> version;
  
  if( !warnCond( !inStream->fail() && !strncmp( buf, "MDA", 3 ),
		 "  Not an MDA file!" ) )
    return false;
  warnCond( version== 1.0,
	    "Incompatible MDA file version - I'll do my best\n" );

 
  // read format and byte order
  *inStream >> ws;
  inStream->width( 16 );
  *inStream >> buf;
  *inStream >> dataType;
  *inStream >> byteOrder;
  if( !warnCond( !inStream->fail() && !strncmp( buf, "Format:", 7 ),
		 "  Cannot read data format!" ) )
    return false;
  
  // read dimension
  *inStream >> ws;
  inStream->width( 16 );
  *inStream >> buf;
  *inStream >> dim;
  if( !warnCond( !inStream->fail() && dim.vec.size()> 0 &&
		 !strncmp( buf, "Dimensions:", 11 ),
		 "  Cannot read coordinate vector!" ) )
    return false;
  
  // read number of channels
  numChannels= 0;
  *inStream >> ws;
  inStream->width( 16 );
  *inStream >> buf >> numChannels;
  if( !warnCond( !inStream->fail() && numChannels> 0 &&
		 !strncmp( buf, "Channels:", 9 ),
		 "  Cannot read number of channels!" ) )
    return false;
  
  // read end of header
  *inStream >> ws;
  inStream->read( buf, 4 );
  
  
  if( !warnCond( !inStream->fail() && !strncmp( buf, "###", 3 ),
		 "  Cannot find end of file header!" ) )
    return false;
  
  // compute scanline info based on the header data, and allocate buffer
  setupScanlineInfo();
  scanline= new char[scanlineSize];

  return true;
}
    
/** read one scanline (NULL: error or not ready) */
void *
MDAReader::readScanline()
{
	
	if (numScanlinesLeft == 0) 
    return NULL;
  numScanlinesLeft--;
  
  inStream->read( scanline, scanlineSize );
  
  
  // swap byte order if necessary
  if( byteOrder!= cpuByteOrder() &&
      byteOrder!= UnknownEndian &&
      dataTypeSizes[dataType]> 1 )
    swapArrayByteOrder( scanline, dim.vec[0]*numChannels,
			dataTypeSizes[dataType] );
  
 return inStream->fail() ? NULL : scanline;
 
}

/** skip one scanline */
bool
MDAReader::skipScanline()
{
  if( numScanlinesLeft== 0 )
	  
    return false;
  numScanlinesLeft--;
  
  // ideally, this would use seekg(), but that does not work on pipes!
  //  inStream->seekg( scanlineSize, ios_base::cur );
  //inStream->read( scanline, scanlineSize );
  inStream->ignore( scanlineSize );
  
  return !inStream->fail();
}

/** check if the file has opposite byte order to CPU
 *  This implies a byte swap will need to occur */
bool
MDAReader::endianSwap()
{
  return byteOrder != cpuByteOrder();
}

/** access stream pointer (some bridge code uses this to make conversions
 *  more efficient) */
istream*
MDAReader::getStreamPointer()
{
  return inStream;
}



  //
  // MDAWriter methods
  //

/** default constructor */
MDAWriter::MDAWriter()
{
#if defined (_WIN32) || defined (_WIN64)
	_setmode(_fileno(stdout), _O_BINARY);
#endif

  ownOstream= false;
  outStream= NULL;
  numScanlinesLeft= 0;
}
    
/** denstructor */
MDAWriter::~MDAWriter()
{
  disconnect();
}

/** connect to provided ostream */
bool
MDAWriter::connect( ostream &os )
{
  outStream= &os;
  outStream->clear();
  ownOstream= false;
  
  return outStream->good();
}

/** constructor: open provided file */
bool
MDAWriter::connect( const char *fileName )
{
#if defined (_WIN32) || defined (_WIN64)
  myOstream.open( fileName, ofstream::out | ofstream::binary );
#else
	myOstream.open(fileName, ofstream::out & ofstream::binary);
#endif
  ownOstream= myOstream.good();
  if( ownOstream )
    outStream= &myOstream;
  
  return ownOstream;
}


/** disconnect from a stream (returns false if write incomplete) */
bool
MDAWriter::disconnect()
{
  bool finished= (numScanlinesLeft== 0);

  numScanlinesLeft= 0;
  if( ownOstream )
  {
    myOstream.close();
    ownOstream= false;
    outStream= NULL;
  }
  
  return finished;
}

/** read MDA header from the istream we are connected to */
bool
MDAWriter::writeHeader( const CoordinateVector &dimension,
			unsigned int channels, DataType dType,
			ByteOrder endianness )
{
  // check if we are connected ot an ostream, and ready to go
	if (outStream == NULL || !outStream->good()) {
		cerr << "Output stream is empty." << endl;
		return false;
	}
  
  // copy data to object and setup scanline information
  dim= dimension;
  numChannels= channels;
  dataType= dType;
  if( endianness== UnknownEndian )
    byteOrder= cpuByteOrder();
  else
    byteOrder= endianness;
  
  // write header
  *outStream << "MDA 1.0\n";
  *outStream << "Format:\t\t";
  *outStream << dataType << ' ' << byteOrder << endl;
  *outStream << "Dimensions:\t" << dim << endl;
  *outStream << "Channels:\t" << numChannels << endl
	     << "###\n";
  if( outStream->fail() )
    return false;
  
  // if everything went well so far, we set up the scanline info
  setupScanlineInfo();
  return true;
}
    
/** read one scanline (NULL: not ready) */
bool
MDAWriter::writeScanline( void *scanline )
{
  if( numScanlinesLeft== 0 )
    return false;
  numScanlinesLeft--;
  
  // swap byte order if necessary
  if( byteOrder!= cpuByteOrder() &&
      byteOrder!= UnknownEndian &&
      dataTypeSizes[dataType]> 1 )
    swapArrayByteOrder( scanline, dim.vec[0]*numChannels,
			dataTypeSizes[dataType] );
  // then write the data
  outStream->write( ( char *)scanline, scanlineSize );
  
  return !(outStream->fail());
}


} /* namespace */

#endif /* MDAFILEIO_C */

