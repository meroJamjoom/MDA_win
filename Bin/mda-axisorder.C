// ==========================================================================
// $Id: mda-axisorder.C 319 2009-05-27 21:17:01Z heidrich $
// change the order of the axis (fliping and transposing) in an MDA stream
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
#include <string.h>

#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "<options> <axis list>\nRearrange the axis in an MDA stream"


int
main( int argc, char *argv[] )
{
  unsigned long i, j;
  
  CommandlineParser parser;
  
  // a list of exes to reverse
  CoordinateVector revertAxes;
  CoordinateOption revertOpt( revertAxes,
			      "\tList of axes to be reflected\n"
			      "\t(indices with respect to order of"
			      " input axes)\n",
			      "--reflect", "-r" );
  parser.registerOption( &revertOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // parse new order of axes
  CoordinateVector axisOrder;
  istringstream axisString( argv[index] );
  axisString >> axisOrder;
  if( axisString.fail() )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // read input header
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
  unsigned int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  
  // check validity of the axis order the user provided
  if( axisOrder.vec.size()!= dimension )
  {
    cerr << "Number of dimensions does not match dimensions in MDA file\n";
    exit( 1 );
  }
  // check axis list for range and uniqueness
  // this is quadratic, but hopefully the number of dimensions is very small...
  for( i= 0 ; i< dimension ; i++ )
  {
    if( axisOrder.vec[i]< 0 || axisOrder.vec[i]>= dimension )
    {
      cerr << "Axes are labelled 0..." << dimension-1 << endl;
      exit( 1 );
    }
    for( j= 0 ; j< i ; j++ )
      if( axisOrder.vec[i]== axisOrder.vec[j] )
      {
	cerr << "Each axis must be listed exactly once in the axis order!\n";
	exit( 1 );
      }
  }
  // check flip list
  // bit vector for the inversion of direction
  bool *reverseOrder= new bool[dimension];
  for( i= 0 ; i< dimension ; i++ )
    reverseOrder[i]= false;
  for( i= 0 ; i< revertAxes.vec.size() ; i++ )
    if( revertAxes.vec[i]>= dimension )
    {
      cerr << "Axis number out of range: MDA stream has only a dimension of "
	   << dimension << endl;
      exit( 1 );
    }
    else
      reverseOrder[revertAxes.vec[i]]= true;
  
  // if all is fine, read the whole input stream into an array
  char *mdaIn= new char[numScanlines*scanlineSize];
  for( i= 0 ; i< numScanlines ; i++ )
    memcpy( mdaIn+i*scanlineSize, reader.readScanline(), scanlineSize );
  
  // output dimensions
  CoordinateVector dimOut;
  for( i= 0 ; i< dimension ; i++ )
    dimOut.vec.push_back( dim.vec[axisOrder.vec[i]] );
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dimOut, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  unsigned long numOutScanlines= writer.getNumScanlinesLeft();
  unsigned long scanlineSizeOut= writer.getScanlineSize();
  char *scanline= new char[scanlineSizeOut];
  int pixelSize= numChannels*dataTypeSizes[type];

  // the pixel-to-pixel increment along all major axes of the
  // index into the original MDA file
  unsigned long *pixelIncrement= new unsigned long[dimension];
  pixelIncrement[0]= pixelSize;
  for( i= 1 ; i< dimension ; i++ )
    pixelIncrement[i]= pixelIncrement[i-1]*dim.vec[i-1];
#if defined (_WIN32) || defined (_WIN64)
   long incr= reverseOrder[axisOrder.vec[0]] ?

	  -static_cast<long>(pixelIncrement[axisOrder.vec[0]]) : pixelIncrement[axisOrder.vec[0]];
#else
  unsigned long incr= reverseOrder[axisOrder.vec[0]] ?

    -pixelIncrement[axisOrder.vec[0]] : pixelIncrement[axisOrder.vec[0]];
#endif
  // create initial position vector
  CoordinateVector pos;
  for( i= 0 ; i< dimension ; i++ )
    pos.vec.push_back( reverseOrder[axisOrder.vec[i]] ? dimOut.vec[i]-1 : 0 );
  
  // actually copy the data
  for( i= 0 ; i< numOutScanlines ; i++ )
  {
    // for each output scanline
    // calculate base index for first pixel in scanline
    index= 0;
    for( j= 0 ; j< dimension ; j++ )
      index+= pos.vec[j] * pixelIncrement[axisOrder.vec[j]];
    
    // collect all pixels for this output scanline
    for( j= 0 ; j< dimOut.vec[0] ; j++ )
    {
      memcpy( scanline+j*pixelSize, mdaIn+index, pixelSize );
      index+= incr;
    }
    
    // write scanline
    writer.writeScanline( scanline );
    
    // update position to first entry of next scanline
    for( j= 1 ; j< dimension ; j++ )
      if( reverseOrder[axisOrder.vec[j]] )
	if( --pos.vec[j]< 0 )
	  pos.vec[j]= dimOut.vec[j]-1;
	else
	  break;
      else
	if( ++pos.vec[j]>= dimOut.vec[j] )
	  pos.vec[j]= 0;
	else
	  break;
  }
  
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
