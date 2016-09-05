// ==========================================================================
// $Id: mda-slice.C 199 2008-04-08 03:58:26Z heidrich $
// split an n-dimensional MDA stream into multiple slices of dimension n-1
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

#include <sstream>
#include <math.h>

#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "<file name template>\n\n\
Cut an n-imensional MDA stream into slices of dimension n-1. The\n\
file names are provided in printf format. I.e. an argument\n\
  test%02d.mda\n\
will result in files test00.mda, test01.mda, etc..\n"

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  char filename[1024];
  
  CommandlineParser parser;
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], USAGE_TEXT );
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
  DataType	type= reader.getType();
  CoordinateVector	dim= reader.getDim();
  unsigned int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  int dimension= dim.vec.size();
  
  // output dimensions...
  CoordinateVector dimOut;
  for( i= 0 ; i< dimension-1 ; i++ )
    dimOut.vec.push_back( dim.vec[i] );
    
  // now go over the individual slices and save them into separate files
  for( i= 0 ; i< dim.vec[dimension-1] ; i++ )
  {
    // create filename
    sprintf( filename, argv[index], i );
    
    // setup output MDA object
    MDAWriter writer;
    writer.connect( filename );
    if( !writer.writeHeader( dimOut, numChannels, type ) )
    {
      cerr << "Cannot write file header\n";
      exit( 1 );
    }
    
    // copy the data for one slice
    while( writer.getNumScanlinesLeft() )
      writer.writeScanline( reader.readScanline() );
  
    if( !writer.disconnect() )
    {
      cerr << "Error while processing data\n";
      exit( 1 );
    }
  }
  
  reader.disconnect();
  
  return 0;
}
