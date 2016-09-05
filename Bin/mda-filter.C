// ==========================================================================
// $Id: mda-filter.C 638 2010-03-08 19:34:24Z heidrich $
// apply a filter to an MDA stream
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


// all arrays are kept in the following type
#define TYPE float


#include "MDA/Config.hh"
#include "MDA/Threading/ThreadingOption.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Filters/Filters.hh"

using namespace MDA;
using namespace std;

char usageText[BUFFER_SIZE]=
 "[<options>] <filter> <params...> [[options] <filter> <params...>]...\n\n"
 "Apply one or more filter(s) to the MDA stream. Different filters require\n"
 "different numbers of parameters. Currently supported filters are:\n\n";

int
main( int argc, char *argv[] )
{
  int i;
  Filter<TYPE> *filter;
  FilterParams filterParams;
  FilterFactory<TYPE> factory( &filterParams );
  
  for( i= 0 ; i< (int)UndefinedFiltering ; i++ )
  {
    strncat( usageText, "  ", BUFFER_SIZE );
    strncat( usageText, filterNames[i], BUFFER_SIZE );
    strncat( usageText, "\n", BUFFER_SIZE );
  }
  
  // setup options
  CommandlineParser parser;
  
  // filter parameters first:

  // 1) filter radius
  filterParams.radius= 0;
  IntOption radiusOpt( filterParams.radius,
		       "\tFilter radius"
		       "(if <=0, choose filter-specific default)\n",
		       "--radius", "-r" );
  parser.registerOption( &radiusOpt );
  
  // 2) sigma
  filterParams.sigma= 1.0;
  DoubleOption sigmaOpt( filterParams.sigma,
			 "\tStandard deviation (filters involving Gaussians)\n",
			 "--sigma", "-s" );
  parser.registerOption( &sigmaOpt );
  
  // 2a) edge stop sigma
  filterParams.edgeStopSigma= 0.1;
  DoubleOption edgeOpt( filterParams.edgeStopSigma,
			"\tEdge stop standard deviation (Bilateral etc.)\n",
			"--edge-stop-sigma", "-ess" );
  parser.registerOption( &edgeOpt );
  
  // 3) filter expressions
  ExpressionOption pixelValOpt( filterParams.values,
				"\tFilter values (seq. of expressions)\n",
				"--filter-values", NULL );
  parser.registerOption( &pixelValOpt );
  
  // 4) normalization
  filterParams.normalize= false;
  BoolOption normalizeOpt( filterParams.normalize,
			   "\tWhether or not to normalize certain filters\n",
			   "--normalize", NULL, "--no-normalize", NULL );
  parser.registerOption( &normalizeOpt );
  
  // other options 
  
  // output type
  DataType	outType= UndefinedType;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  // boundary mode
  BoundaryMethod boundary= Clamp;
  BoundaryOption boundaryOpt( boundary );
  parser.registerOption( &boundaryOpt );
  
  // the affected coordinate axes
  AxisList axes;
  ChannelListOption axesOpt( axes,
			     "\tAffected coordinate axes (default: all)\n",
			     "--axislist", "-al" );
  parser.registerOption( &axesOpt );
  
  // the affected channels
  ChannelList channels;
  ChannelListOption channelOpt( channels,
				"\tAffected channels (default: all)\n" );
  parser.registerOption( &channelOpt );

  // how often to repeat each filter
  int repeat= 0;
  IntOption repeatOpt( repeat,
		       "\tHow often to repeat each individual filter\n",
		       "--repeat", NULL );
  parser.registerOption( &repeatOpt );
  
  // whether to use threading or not
  ThreadingOption threading;
  parser.registerOption( &threading );
  
  // read MDA stream
  Array<TYPE> array;
  if( !array.read() )
  {
    cerr << argv[0] << ": Cannot read input MDA stream!\n\n";
    exit( 1 );
  }
  
  // while there are still commandline arguments left,
  // read options and filter names, and apply them in sequence
  int index= 1;
  while( index< argc )
  {
    if( !parser.parse( index, argc, argv ) || index> argc-1 )
    {
      parser.usage( argv[0], usageText );
      exit( 1 );
    }

    // parse the filter
    filter= NULL;
    FilterType filterType= UndefinedFiltering;
    string filterString= argv[index++];
    istringstream filterStream( filterString );
    filterStream >> filterType;
    
    // and create it
    filter= factory.create( filterType );
    
    // fill in channel and axes ranges if not user provided
    if( channels.vec.size()== 0 )
      channels= array.allChannels();
    if( axes.vec.size()== 0 )
      axes= array.allAxes();
    
    // apply the filter, possibly multiple times
    filter->apply( array, boundary, channels, axes );
    for( i= 1 ; i< repeat ; i++ )
      filter->apply( array, boundary, channels, axes );
    
    // clean up
    if( filter!= NULL )
      delete filter;
  }
  
  // write MDA stream
  array.write( cout, outType );
  return 0;
}
