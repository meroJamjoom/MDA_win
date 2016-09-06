// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef COLOR_WHITEBALANCE_C
#define COLOR_WHITEBALANCE_C

#include "WhiteBalance.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** constructor from named color spaces */
WhiteBalance::WhiteBalance( const Whitepoint &srcWP, const Whitepoint &dstWP,
			    const char *srcSpace, const char *vonKriesSpace )
{
  // create color spaces
  src= new LinearTristimulusSpace( srcSpace );
  vonKries= new LinearTristimulusSpace( vonKriesSpace );
  
  // get XYC coordinates of source and destination whitepoints
  Vector srcXYZ( srcWP.getXYZ() );
  Vector dstXYZ( dstWP.getXYZ() );
  
  // compute multipliers
  mult[0]= dstXYZ[0] / srcXYZ[0];
  mult[1]= dstXYZ[1] / srcXYZ[1];
  mult[2]= dstXYZ[2] / srcXYZ[2];
  cerr << mult[0] << ' '  << mult[1] << ' '  << mult[2] << ' ' << endl;
}

/** constructor from externally created color spaces */
WhiteBalance::WhiteBalance( const Whitepoint &srcWP, const Whitepoint &dstWP,
			    ColorSpace *srcSpace, ColorSpace *vonKriesSpace )
  : src( srcSpace ), vonKries( vonKriesSpace )
{
  // get XYC coordinates of source and destination whitepoints
  Vector srcXYZ( srcWP.getXYZ() );
  Vector dstXYZ( dstWP.getXYZ() );
  
  // compute multipliers
  mult[0]= dstXYZ[0] / srcXYZ[0];
  mult[1]= dstXYZ[1] / srcXYZ[1];
  mult[2]= dstXYZ[2] / srcXYZ[2];
}

/** apply whitebalance (result overwrites input) */
void
WhiteBalance::apply( Vector &color )
{
  Vector XYZ( 3 );
  Vector abc( 3 );
  
  // convert into von Kries space
  src->toXYZ( color, XYZ );
  vonKries->fromXYZ( XYZ, abc );
  // white balance
  abc[0]*= mult[0];
  abc[1]*= mult[1];
  abc[2]*= mult[2];
  // convert back to original space
  vonKries->toXYZ( abc, XYZ );
  src->fromXYZ( XYZ, color );
}


} /* namespace */

#endif /* COLOR_WHITEBALANCE_C */

