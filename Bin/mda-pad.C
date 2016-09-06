// ==========================================================================
// $Id: mda-pad.C 319 2009-05-27 21:17:01Z heidrich $
// pad an array with user-provided pixel data
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
#include <string.h>
#include <math.h>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "<options>\n\n\
Pad the MDA array with a (potentially asymmetric) frame\n"

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
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // pixel value
  EXPR::ExpressionParseTree defaultPixel;
  defaultPixel.makeConstant( 0.0 );
  EXPR::ExpressionSequence color;
  color.push_back( &defaultPixel );
  ExpressionOption colorOpt( color,
			     "\tpixel value\n\t(number of expressions must"
			     " match number of channels, or be 1)\n",
                             "--pixel-value", "-p" );
  parser.registerOption( &colorOpt );
  
  // frame thickness
  int frame= 0;
  IntOption frameOption( frame,
			 "\tthickness of a symmetric frame in pixels (>=0)\n",
			 "--frame", "-f", 0 );
  parser.registerOption( &frameOption );
  
  // rows before the current data
  CoordinateVector before;
  CoordinateOption beforeOption( before,
				 "\tnumber of rows to add before the data "
				 "along each dimension\n", "--before", NULL );
  parser.registerOption( &beforeOption );
  
  // rows after the current data
  CoordinateVector after;
  CoordinateOption afterOption( after,
				"\tnumber of rows to add after the data "
				"along each dimension\n", "--after", NULL );
  parser.registerOption( &afterOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
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
  unsigned int dimension= dim.vec.size();
  int numChannels= reader.getNumChannels();  
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  char *scanline;
  
  // process user-provided pixel data such that there is one scalar per channel
  unsigned numColorChannels= color.size();
  for( i= 0 ; i< numColorChannels ; i++ )
    if( color[i]->getMaxVariable( true )>= 0 ||
	color[i]->getMaxVariable( false )>= 0 )
    {
      cerr << "No variable or registers supported\n";
      exit( 1 );
    }
  double *pixel= new double[numChannels];
  if( numColorChannels== 1 )
  {
    // replaicate the one exoression across all channels
    double h= color[0]->eval( NULL, NULL );
    for( i= 0 ; i< numChannels ; i++ )
      pixel[i]= h;
  }
  else if( numColorChannels== numChannels )
    // evaluate a different expression for each channel
    for( i= 0 ; i< numChannels ; i++ )
      pixel[i]= color[i]->eval( NULL, NULL );
  else
  {
    cerr << "The number of channels in the pixel value needs to be either\n"
	 << "one, or identical to the number of channels in the stream\n";
    exit( 1 );
  }
  // convert pixel value to array data type
  char *pixelData= new char[numChannels*dataTypeSizes[type]];
  typeConvert( pixel, Double, pixelData, type, numChannels );
  
  // convert frame option into before/after format
  if( frame> 0 )
  {
    if( before.vec.size()!= 0 || after.vec.size()!= 0 )
    {
      cerr << "Warning: both --frame and --before/after specified:"
	   << "using --frame\n";
      before.vec.clear();
      after.vec.clear();
    }
    for( i= 0 ; i< dimension ; i++ )
    {
      before.vec.push_back( frame );
      after.vec.push_back( frame );
    }
  }
  // fill in before/after values with zero if they have not been provided
  for( i= before.vec.size() ; i< dimension ; i++ )
    before.vec.push_back( 0 );
  for( i= after.vec.size() ; i< dimension ; i++ )
    after.vec.push_back( 0 );
  
  // create pixel range corresponding to original data in the output
  ULongRangeList originalRegion;
  for( i= 0 ; i< dimension ; i++ )
    originalRegion.vec.push_back( ULongRange( before.vec[i],
					      dim.vec[i]+before.vec[i]-1 ) );
  
  // the new dimensions of the padded data
  CoordinateVector dimOut= dim;
  for( i= 0 ; i< dimension ; i++ )
    dimOut.vec[i]+= before.vec[i]+after.vec[i];
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dimOut, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  int outScanlineSize= writer.getScanlineSize();
  
  // create a scanline initialized with the pixel value
  char *emptyScanline= new char[outScanlineSize];
  for( i= 0 ; i< dimOut.vec[0] ; i++ )
    memcpy( emptyScanline+i*numChannels*dataTypeSizes[type], pixelData,
	    numChannels*dataTypeSizes[type] );
  
  // initialize the scanline buffer to also be empty
  char *outScanline= new char[outScanlineSize];
  memcpy( outScanline, emptyScanline, outScanlineSize );
  
  // the actual copy operation
  unsigned long scanlineOffset= before.vec[0]*numChannels*dataTypeSizes[type];  
  CoordinateVectorIter pos( dimOut );
  for( pos.begin() ; !pos.isAtEnd() ; pos.incrComp( 1 ) )
  {
    if( inRange( pos.getPos(), originalRegion ) )
    {
      // read scanline
      scanline= (char *)reader.readScanline();
      
      // and write it together with the padding
      memcpy( outScanline+scanlineOffset, scanline, scanlineSize );
      writer.writeScanline( outScanline );
    }
    else
      writer.writeScanline( emptyScanline );
  }
  
  // close reader, writer
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  // cleanup
  delete [] outScanline;
  delete [] emptyScanline;
  delete [] pixelData;
  delete [] pixel;
  
  return 0;
}
