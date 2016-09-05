// ==========================================================================
// $Id: mda-fromraw.C 199 2008-04-08 03:58:26Z heidrich $
// extract MDA information from a raw file
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


#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  CommandlineParser parser;
  
  DataType	type= UByte;
  TypeOption	typeOption( type );
  parser.registerOption( &typeOption );
  
  ByteOrder	cpuOrder= cpuByteOrder();
  ByteOrder	order= cpuOrder;
  ByteOrderOption	orderOption( order );
  parser.registerOption( &orderOption );
  
  CoordinateVector	dim;
  CoordinateOption coordOption( dim );
  parser.registerOption( &coordOption );
  
  int	numChannels= 1;
  IntOption	chanOption( numChannels,
			    "\tchannel count (>0!)\n",
			     "--channels", "-c", 1 );
  parser.registerOption( &chanOption );
  
  long offset= 0;
  LongOption offsetOption( offset,
	  "\toffset of raw data in the file (negative: start from back)\n",
			  "--offset", "-off" );
  parser.registerOption( &offsetOption );
  
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "<options> <raw file>" );
    exit( 1 );
  }
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannels, type ) )
  {
    cerr << "Cannot open write file header\n";
    exit( 1 );
  }
  unsigned long numScanlines= writer.getNumScanlinesLeft();
  unsigned long scanlineSize= writer.getScanlineSize();
  char *scanline= new char[scanlineSize];
  cerr << scanlineSize << endl;

  // setup input stream
  ifstream is;
  is.open( argv[argc-1], ifstream::in & ifstream::binary );
  if( is.fail() )
  {
    cerr << "Cannot open " << argv[argc-1] << " for reading\n";
    exit( 1 );
  }
  
  if( offset< 0 )
    is.seekg( -(long)(numScanlines*scanlineSize)+offset+1, ios_base::end );
  else if( offset > 0 )
    is.seekg( offset );
  if( is.fail() )
  {
    cerr << "Cannot seek to offset position.\n";
    exit( 1 );
  }
  
  // actually copy the data
  for( unsigned long i= numScanlines ; i> 0 ; i-- )
  {
    is.read( scanline, scanlineSize );
    if( order!= cpuOrder )
      swapArrayByteOrder( scanline, numChannels*dim.vec[0],
			  dataTypeSizes[type] );
    writer.writeScanline( scanline );
  }
  if( is.fail() || !writer.disconnect() )
  {
    cerr << "Error while copying\n";
    exit( 1 );
  }
  
  return 0;
}
