// ==========================================================================
// $Id: mda-deinterlace.C 638 2010-03-08 19:34:24Z heidrich $
// read a sequence of interlaced MDAs, and de-interlace them
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

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "MDA/Config.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "[<options>]\n"


bool
loadArray( unsigned long &width, unsigned long &height, int &numChannels,
	   unsigned char *&buffer, bool firstTime= false )
{
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() )
    return false;
  
  // first time setup
  CoordinateVector dim;
  if( firstTime )
  {
    dim= reader.getDim();
    if( dim.vec.size()== 2 )
    {
      width= dim.vec[0];
      height= dim.vec[1];
    }
    numChannels= reader.getNumChannels();
  }
  
  // run a consistency check on the MDA
  dim= reader.getDim();
  if( dim.vec.size()!= 2 )
  {
    cerr << "Expecting 2D array!\n";
    return false;
  }
  if( dim.vec[0]!= width || dim.vec[1]!= height )
  {
    cerr << "Dimensions do not match previous frames\n";
    return false;
  }
  if( reader.getNumChannels()!= numChannels )
  {
    cerr << "Number of channels does not match previous frames\n";
    return false;
  }
  
  // allocate buffer if not yet present
  if( buffer== NULL )
    buffer= new unsigned char[width*height*numChannels];
  
  // read the MDA contents into the buffer
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  for( unsigned long i= 0 ; i< numScanlines ; i++ )
    memcpy( buffer+i*numChannels*width,
	    reader.readScanline(), numChannels*width );
  
  reader.disconnect();
  
  return true;
}


void
interpolateScanline( unsigned char *spatialPrev, unsigned char *spatialNext,
		     unsigned char *temporalPrev, unsigned char *temporalNext,
		     unsigned char *newCurrent, int numChannels,
		     unsigned long width, double sigma )
{
  unsigned long i, j;
  
  for( i= 0 ; i< width ; i++ )
  {
    double diff= 0.0;
    for( j= 0 ; j< numChannels ; j++ )
      diff+=
	(temporalNext[j]-temporalPrev[j]) * (temporalNext[j]-temporalPrev[j]);
    double weight= exp( -diff/(2*sigma*sigma) );
    for( j= 0 ; j< numChannels ; j++ )
      *newCurrent++= (unsigned char)
	(.5 * (weight * (*temporalNext++ + *temporalPrev++) +
	       (1.0-weight) * (*spatialNext++ + *spatialPrev++)));
  }
}

int
main( int argc, char *argv[] )
{
  char	fileName[BUFFER_SIZE];
  char  cmdName[BUFFER_SIZE];
  unsigned long	i;
  
  CommandlineParser parser;
  
  bool even= true;
  BoolOption evenOption( even,
			 "\twhether the first frame is even or odd\n",
			 "--even", NULL, "--odd", NULL );
  parser.registerOption( &evenOption );
  
  float sigma= 1.0;
  FloatOption sigmaOption( sigma,
			   "\tGaussian sigma to use in similarity between "
			   "scanlines\n\t(increase for noisy video streams)\n",
			   "--sigma", NULL );
  parser.registerOption( &sigmaOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  
  unsigned long width;
  unsigned long height;
  int numChannels;
  unsigned char *buffers[3]= {NULL, NULL, NULL};
  if( !loadArray( width, height, numChannels, buffers[0], true ) )
  {
    cerr << argv[0] << ": Cannot load first frame\n";
    exit( 1 );
  }
  if( !loadArray( width, height, numChannels, buffers[1] ) )
  {
    cerr << argv[0] << ": Cannot load second frame\n";
    exit( 1 );
  }
  unsigned long numEntries= width*numChannels;
  
  // output stuff
  unsigned char *scanlineOut= new unsigned char[numEntries];
  CoordinateVector dimOut;
  dimOut.vec.push_back( width );
  dimOut.vec.push_back( height*2 );

  // load as many additional frames as we can, and perform the
  // de-interlacing
  while( loadArray( width, height, numChannels, buffers[2] ) )
  {
    // toggle odd/even
    even= !even;
    
    // set up writer
    MDAWriter writer;
    writer.connect( cout );
    if( !writer.writeHeader( dimOut, numChannels, UByte ) )
    {
      cerr << argv[0] << ": cannot write header\n";
      exit( 1 );
    }
    
    if( even )
      // handle even frames
      for( i= 0 ; i< height ; i++ )
      {
	// write a scanline from the current timestep
	writer.writeScanline( buffers[1]+i*numEntries );
	
	// interpolate the next scanline and write it
	interpolateScanline( buffers[1]+i*numEntries,
			     buffers[1]+(i==height-1 ? i : i+1)*numEntries,
			     buffers[0]+i*numEntries,
			     buffers[2]+i*numEntries,
			     scanlineOut, numChannels, width, sigma );
	writer.writeScanline( scanlineOut );
      }
    else
      // handle odd frames
      for( i= 0 ; i< height ; i++ )
      {
	// interpolate the next scanline and write it
	interpolateScanline( buffers[1]+(i== 0 ? i : i-1)*numEntries,
			     buffers[1]+i*numEntries,
			     buffers[0]+i*numEntries,
			     buffers[2]+i*numEntries,
			     scanlineOut, numChannels, width, sigma );
	writer.writeScanline( scanlineOut );
	
	// write a scanline from the current timestep
	writer.writeScanline( buffers[1]+i*numEntries );
      }
      
    if( !writer.disconnect() )
    {
      cerr << argv[0] << ": error during writing\n";
      exit( 1 );
    }
    
    // swap buffers
    unsigned char *tmp= buffers[0];
    buffers[0]= buffers[1];
    buffers[1]= buffers[2];
    buffers[2]= tmp;
  }
  
  return 0;
}
