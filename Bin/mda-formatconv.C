// ==========================================================================
// $Id: mda-formatconv.C 260 2008-12-30 09:21:39Z heidrich $
// change data type or byte order of an MDA stream
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

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  unsigned long i;
  
  CommandlineParser parser;
  
  // setup options
  DataType	outType= UByte;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  ByteOrder	outOrder= cpuByteOrder();
  ByteOrderOption	orderOption( outOrder );
  parser.registerOption( &orderOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) )
  {
    parser.usage( argv[0] );
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
  char *inScanline;
  
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannels, outType, outOrder ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  unsigned long outScanlineSize= writer.getScanlineSize();
  char *outScanline= new char[outScanlineSize];
  
  
  // actually copy the data
  unsigned int bytesPerValue= dataTypeSizes[outType];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    inScanline= (char *)reader.readScanline();
    
    // do the type conversion
    typeConvert( inScanline, inType, outScanline, outType,
		 dim.vec[0]*numChannels );
    
    // writer takes care of byte order
    writer.writeScanline( outScanline );
  }
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
