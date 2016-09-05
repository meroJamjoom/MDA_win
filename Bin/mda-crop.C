// ==========================================================================
// $Id: mda-crop.C 199 2008-04-08 03:58:26Z heidrich $
// crop array in an MDA file
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

#define USAGE_TEXT "<crop range>\n\n\
Crop the MDA stream to the provided range. The length of\n\
the range list must match the dimension of the array.\n"

inline bool
inRange( const CoordinateVector &pos, const ULongRangeList &cropRegion )
{
  for( int i= pos.vec.size()-1 ; i> 0 ; i-- )
    if( pos.vec[i]< cropRegion.vec[i].val.first ||
	pos.vec[i]> cropRegion.vec[i].val.second )
      return false;
  return true;
}

int
main( int argc, char *argv[] )
{
  unsigned long i, j;
  
  CommandlineParser parser;
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // get the crop region as a range list
  ULongRangeList cropRegion;
  istringstream rangeString( argv[index] );
  rangeString >> cropRegion;
  if( rangeString.fail() )
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
  char *scanline;
  if( dim.vec.size()!= cropRegion.vec.size() )
  {
    cerr << "Dimension of crop region " << cropRegion
	 << " does not match dimension of MDA file (" << dim << ")\n";
    exit( 1 );
  }
  
  // now go over the individual ranges and crop them to the dimension
  // of the MDA stream, also compute the new dimensions
  int dimension= dim.vec.size();
  CoordinateVector dimOut;
  for( int i= 0 ; i< dimension ; i++ )
  {
    ULongRange cRange;
    cRange.val.first= 0; cRange.val.second= dim.vec[i]-1;
    cropRegion.vec[i].val.first=
      clampToRange( cropRegion.vec[i].val.first, cRange );
    cropRegion.vec[i].val.second=
      clampToRange( cropRegion.vec[i].val.second, cRange );
    dimOut.vec.push_back( cropRegion.vec[i].val.second-
			  cropRegion.vec[i].val.first+1 );
  }
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dimOut, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  
  
  // actually copy the data
  unsigned long scanlineOffset=
    cropRegion.vec[0].val.first*numChannels*dataTypeSizes[type];
  CoordinateVectorIter pos( dim );
  for( pos.begin() ; !pos.isAtEnd() ; pos.incrComp( 1 ) )
  {
    if( inRange( pos.getPos(), cropRegion ) )
    {
      // read scanline
      scanline= (char *)reader.readScanline();
      
      // and write part of it
      writer.writeScanline( scanline + scanlineOffset );
    }
    else
      // skip scanline
      reader.skipScanline();
  }
  
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
