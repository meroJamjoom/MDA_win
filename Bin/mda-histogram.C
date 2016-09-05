// ==========================================================================
// $Id: mda-histogram.C 199 2008-04-08 03:58:26Z heidrich $
// compute a histogram of array values
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2008-, UBC
//
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#include <assert.h>
#include <time.h>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // number of bins
  int numBins= 256;
  IntOption binOption( numBins, "\tthe number of bins\n",
		       "--numn-bins", "-n" );
  parser.registerOption( &binOption );
  
  // mapping function
  EXPR::ExpressionSequence map= parse( "%0*#0" );;
  ExpressionOption mapOpt( map,
		   "\ta function mapping channel values to histogram bins\n"
		   "\t\t(%0 is the number of bins, #0 the current channel)\n",
			   "--map" );
  parser.registerOption( &mapOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], "[options]" );
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
  int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long inScanlineSize= reader.getScanlineSize();
  char *inScanline;
  
  // histogram data
  unsigned int *hist= new unsigned int[numBins*numChannels];
  for( i= 0 ; i< numBins*numChannels ; i++ )
    hist[i]= 0;
  
  // make sure the number of mapper expression matches the number of
  // channels, and all register indices are valid
  int numMaps= map.size();
  if( numMaps != 1 && numMaps!= numChannels )
  {
    cerr << "Number of mapping functions needs to match"
	 << " number of channels (or be 1)";
    exit( 1 );
  }
  for( i= 0 ; i< numMaps ; i++ )
    if( map[i]->getMaxVariable( true )> 0 ||
	map[i]->getMaxVariable( false )> 0 )
    {
      cerr << "Only registers #0 (channel value) and %0 (number of bins)"
	   << " are valid\n";
      exit( 1 );
    }
  
  // actually map the data and update the histogram
  unsigned typeSize= dataTypeSizes[inType];
  double inValue;
  double bins= numBins; // for calculations
  double histValue;
  unsigned int bin;
  for( i= numScanlines ; i> 0 ; i-- )
  {
    // read scanline
    inScanline= (char *)reader.readScanline();
    
    // foreach pixel in the scanline
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      // for each channel
      for( k= 0 ; k< numChannels ; k++ )
      {
	// convert pixel value to double
	typeConvert( inScanline+(j*numChannels+k)*typeSize, inType,
		     &inValue, Double, 1 );
	
	// compute histogram bin and update it
	if( numMaps> 1 )
	  histValue= map[k]->eval( &inValue, &bins );
	else
	  histValue= map[0]->eval( &inValue, &bins );
	bin= histValue< 0 ? 0 :
	  (histValue> numBins-1 ? numBins-1 : (int)histValue);
	hist[bin*numChannels+k]++;
      }
    }
  }
  reader.disconnect();

  // setup output MDA object & write histogram
  MDAWriter writer;
  writer.connect( cout );
  CoordinateVector outDim;
  outDim.vec.push_back( numBins );
  if( !writer.writeHeader( outDim, numChannels, UInt ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  writer.writeScanline( hist );
  assert( writer.getNumScanlinesLeft()== 0 );
  unsigned long outScanlineSize= writer.getScanlineSize();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  delete [] hist;
  
  return 0;
}
