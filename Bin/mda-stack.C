// ==========================================================================
// $Id: mda-stack.C 638 2010-03-08 19:34:24Z heidrich $
// stack up several MDA files of dimension n-1 to one n-dimensional one
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

#include "MDA/Config.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "[<options>] [<slice files>]\n\n\
Stack up several MDA files of dimension n-1 to one n-dimensional one."

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  char filename[BUFFER_SIZE];
  
  CommandlineParser parser;

  // number of slices to expect when reading from a stream
  int numStreamSlices= 1;
  IntOption numStreamSlicesOption( numStreamSlices,
      "\tnumber of slices (if reading from a stream)\n",
      "--num-slices", "-ns", 0 );
  parser.registerOption( &numStreamSlicesOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  bool fromFiles= (index!= argc);
  int numSlices= fromFiles ? argc-index : numStreamSlices;
  
  
  unsigned long numInScanlines, numScanlines, scanlineSize;
  int numChannels;
  DataType type;
  CoordinateVector dimIn;
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );

  // go over each slice and add its contents to the output stream
  i= 0;
  while( (fromFiles) ? index<argc : cin.peek() && cin.good() )
  {
    // setup MDA reader and read input header
    MDAReader reader;
    if( fromFiles ) {
      if( !reader.connect(argv[index]) ) {
        cerr << "Cannot connect to MDA input file" << endl;
        exit( 1 );
      }
    } else {
      if( !reader.connect(cin) ) {
        cerr << "Cannot connect to input stream" << endl;
        exit( 1 );
      }
    }
    if( !reader.readHeader() )
    {
      cerr << "Cannot read file header\n";
      exit( 1 );
    }
    
    if( i== 0 )
    {
      // first slice: extract information on dimension and type, set
      // up output stream
      type= reader.getType();
      dimIn= reader.getDim();
      numChannels= reader.getNumChannels();
      numScanlines= reader.getNumScanlinesLeft();
      scanlineSize= reader.getScanlineSize();
      
      CoordinateVector dimOut= dimIn;
      dimOut.vec.push_back( numSlices );
      if( !writer.writeHeader( dimOut, numChannels, type ) )
      {
        cerr << "Cannot write file header\n";
        exit( 1 );
      }
    }
    else
    {
      // all foolowing slices: make sure they have the same format as
      // the first one
      if( type!= reader.getType() ||
          dimIn.vec!= reader.getDim().vec ||
          numChannels!= reader.getNumChannels() )
      {
        cerr << argv[index] << " has different layout - cannot process\n";
        exit( 1 );
      }
    }
    
    // copy the scanlines for this slice
    for( j= 0 ; j< numScanlines ; j++ )
      writer.writeScanline( reader.readScanline() );
    reader.disconnect();

    i++;
    if( fromFiles ) {
      index++;
    }
  }
  
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
