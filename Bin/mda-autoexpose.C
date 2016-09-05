// ==========================================================================
// $Id: mda-resize.C 172 2007-09-09 06:19:09Z heidrich $
// autoexpose a number of channels in an MDA stream
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


#include "MDA/Array/Array.hh"
#include "MDA/Resampling/Resampling.hh"
#include "MDA/GeometricTransform/GeometricTransform.hh"

using namespace MDA;
using namespace std;

#define MAX_TEXTLEN 1024
char usageText[MAX_TEXTLEN]= "[<options>...]";



int
main( int argc, char *argv[] )
{
  double multiplier;
  unsigned long i, j;
  
  // setup options
  CommandlineParser parser;
  
  // setup options
  ChannelList   channels;
  ChannelListOption channelOption( channels );
  parser.registerOption( &channelOption );
  
  // output type
  DataType	outType= UndefinedType;
  TypeOption	typeOption( outType,
			    "\tOutput type (default: same as input)\n" );
  parser.registerOption( &typeOption );
  
  // white balance psoition
  CoordinateVector whiteBalancePos;
  CoordinateOption whiteBalanceOpt( whiteBalancePos,
				    "\timage position used for white"
				    "balancing\n",
				    "--white-balance-pos", "-wb" );
  parser.registerOption( &whiteBalanceOpt );
  
  // exposure correction
  double exposureCorrection= 0.0;
  DoubleOption correctOpt( exposureCorrection,
			   "\texposure correction (in stops)\n",
			   "--correct", NULL );
  parser.registerOption( &correctOpt );
  
  // process commandline options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], usageText );
    exit( 1 );
  }
  
  // read array
  Array<double> array;
  if( !array.read() )
  {
    cerr << argv[0] << ": Cannot read input MDA stream!\n\n";
    exit( 1 );
  }
  CoordinateVector dim= array.getDimension();
  
  // compute number of pixels/channel
  unsigned long numPixels= 1;
  for( i= 0 ; i< dim.vec.size() ; i++ )
    numPixels*= dim.vec[i];
  
  // if no channel list supplied, then we proces all channels
  if( channels.vec.size()== 0 )
    for( i= 0 ; i< array.getNumChannels() ; i++ )
      channels.vec.push_back( i );
  
  // white balance first
  if( whiteBalancePos.vec.size()> 0 )
  {
    if( whiteBalancePos.vec.size()!= dim.vec.size() )
    {
      cerr << "White balancing position does not match array dimensions\n";
      exit( 1 );
    }
    
    // compute array offset for wb position
    unsigned long wbOffset= whiteBalancePos.vec[whiteBalancePos.vec.size()-1];
    for( i= whiteBalancePos.vec.size()-2 ; i>= 0 ; i-- )
      wbOffset= wbOffset*dim.vec[i] + whiteBalancePos.vec[i];
    
    // apply white balance
    for( i= 0 ; i< channels.vec.size() ; i++ )
    {
      multiplier= (*array[i])[wbOffset];
      for( j= 0 ; j< numPixels ; j++ )
	(*array[i])[j]/= multiplier;
    }
  }
  
  // compute mean
  double mean= 0.0;
  for( i= 0 ; i< channels.vec.size() ; i++ )
    for( j= 0 ; j< numPixels ; j++ )
      mean+= (*array[i])[j];
  mean/= numPixels*channels.vec.size();
  
  // multiplier to set mean to 18% + exposure correction
  multiplier= 0.18/mean * pow( 2.0, exposureCorrection );
  
  // apply multiplier
  for( i= 0 ; i< channels.vec.size() ; i++ )
    for( j= 0 ; j< numPixels ; j++ )
      (*array[i])[j]*= multiplier;
  
  // write out array
  array.write( cout, outType );
  
  return 0;
}
