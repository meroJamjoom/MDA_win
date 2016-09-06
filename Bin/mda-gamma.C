// ==========================================================================
// $Id: mda-gamma.C 260 2008-12-30 09:21:39Z heidrich $
// gamma correct some channels in an MDA file
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

#include <assert.h>
#include <math.h>

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


// gamma curve with linear segment for dark tones
inline double
gammaFunc( double value, double gamma, double threshold, double slope,
	   double gain, double bias, bool clamp= true )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( value< 0.0 )
    if( clamp )
      value= 0.0;
    else
    {
      inversion= -1.0;
      value= -value;
    }
  // clamp against 1
  if( value> 1.0 && clamp )
    value= 1.0;
  
  if( value<= threshold )
    return inversion * value * slope;
  else
    return inversion * (gain * pow( value, 1.0/gamma ) + bias);
}


// inverse gamma function with linear segment
inline double
gammaFuncInv( double value, double gamma, double threshold, double slope,
	      double gain, double bias, bool clamp= true )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( value< 0.0 )
    if( clamp )
      value= 0.0;
    else
    {
      inversion= -1.0;
      value= -value;
    }
  // clamp against 1
  if( value> 1.0 && clamp )
    value= 1.0;
  
  if( value<= slope*threshold )
    return inversion * value / slope;
  else
    return inversion * pow( (value - bias) / gain, gamma );
}


int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // setup options
  ChannelList	channels;
  ChannelListOption channelOption( channels );
  parser.registerOption( &channelOption );
  
  // forward or backward?
  bool forward= true;
  BoolOption forwardOption( forward,
          "\tmap linear space to gamma space, or vice versa?\n",
			    "--linear-to-gamma", "-l2g",
			    "--gamma-to-linear", "-g2l" );
  parser.registerOption( &forwardOption );
  
  // gamma value
  float gamma= 2.2;
  FloatOption gammaOption( gamma, "\tgamma value\n", "--gamma", NULL );
  parser.registerOption( &gammaOption );
  
  // threshold
  float threshold= 0.0;
  FloatOption thresholdOption( threshold,
			       "\tthreshold for linear vs. gamma curve\n",
			       "--threshold", NULL );
  parser.registerOption( &thresholdOption );
  
  // slope
  float slope= 1.0;
  FloatOption slopeOption( slope,
			   "\tslope of the linar portion of the curve\n",
			   "--slope", NULL );
  parser.registerOption( &slopeOption );
  
  // gain
  float gain= 1.0;
  FloatOption gainOption( gain, "\tgain for gamma curve\n", "--gain", NULL );
  parser.registerOption( &gainOption );
  
  // bias
  float bias= 0.0;
  FloatOption biasOption( bias, "\tbias for gamma curve\n", "--bias", NULL );
  parser.registerOption( &biasOption );
  
  // forward or backward?
  bool clamp= true;
  BoolOption clampOption( clamp,
			  "\tclamp values to 0...1?\n"
			  "\t(some standards allow "
			  "for values outside this range)\n",
			  "--clamp", NULL,
			  "--no-clamp", NULL );
  parser.registerOption( &clampOption );
  
  int standard= -1;
  list<const char *>standards;
  standards.push_back( "--srgb" );
  standards.push_back( "--bt709" );
  standards.push_back( "--xyYCC" );
  SelectionOption standardOpt( standard,
			       "\tselect parameters according to image/video"
			       " standard:\n"
			       "\tsRGB, ITU-R BT.709 (HDTV), or xyYCC (wide "
			       "gamut HDTV)\n",
			       standards );
  parser.registerOption( &standardOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0],
		  "\n\nApply or un-apply a gamma curve to some channels.\n"
		  "Some standards require a linear segment for dark values.\n"
		  "The overall formula is:\n"
		  "\tL_out = slope * L_in\t\t\t; if 0 <= L_in < threshold\n"
		  "\tL_out =  gain * L_in^(1/gamma) + bias\t; else\n" );
    exit( 1 );
  }
  
  // apply standard curve if specified
  switch( standard )
  {
  case 0: // sRGB
    threshold= 0.00304;
    slope= 12.92;
    bias= -0.055;
    gain= 1.055;
    gamma= 2.4;
    break;
  case 2: // xvYCC: like BT.709 but without clamping
    clamp= false;
  case 1: // BT.709
    threshold= 0.018;
    slope= 4.5;
    bias= -0.099;
    gain= 1.099;
    gamma= 1.0 / 0.45;
    break;
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
  

  // check if all channels in the channel list are vaild
  int numGammaChannels= channels.vec.size();
  if( numGammaChannels== 0 )
    // if no channels are specified, we apply gamma to all of them
    for( ; numGammaChannels< numChannels ; numGammaChannels++ )
      channels.vec.push_back( numGammaChannels );
  else
    for( i= 0 ; i< numGammaChannels ; i++ )
      if( channels.vec[i]>= numChannels )
      {
	cerr << "Input MDA only has " << numChannels << " channels\n";
	exit( 1 );
      }
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  assert( scanlineSize== writer.getScanlineSize() );
  
  
  // actually gamma-correct the data
  unsigned int bytesPerValue= dataTypeSizes[type];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    // read scanline
    scanline= (char *)reader.readScanline();
    
    // for each pixel, gamma correct each selected channel
    for( j= 0 ; j< dim.vec[0] ; j++ )
      for( k= 0 ; k< numGammaChannels ; k++ )
      {
	double value;
	
	// convert current channel value to double
	typeConvert( scanline+ (j*numChannels + channels.vec[k])*bytesPerValue,
		     type, &value, Double, 1 );
	
	// perform gamma correction
	double sign= value< 0.0 ? -1.0 : 1.0;
	if( forward )
	  value= sign *
	    gammaFunc( sign*value, gamma, threshold, slope, gain, bias );
	else
	  value= sign *
	    gammaFuncInv( sign*value, gamma, threshold, slope, gain, bias );
	
	// convert gamma corrected value back to original type
	typeConvert( &value, Double,
		     scanline+ (j*numChannels + channels.vec[k])*bytesPerValue,
		     type, 1 );
      }
    
    writer.writeScanline( scanline );
  }
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
