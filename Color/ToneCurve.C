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

#ifndef COLOR_TONECURVE_C
#define COLOR_TONECURVE_C

#include "ToneCurve.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** conversion of color vector from linear to non-linear space
    (identical curve for all channels) */
void
ToneCurve::fromLinear( const Vector &linear, Vector &nonlinear )
{
  errorCond( linear.getSize()== nonlinear.getSize(),
	     "  color space dimensions don't match!" );
  
  for( int i= linear.getSize() ; i-->= 0 ; )
    nonlinear[i]= fromLinear( linear[i] );
}
    
/** conversion of color vector from non-linear to linear space
    (identical curve for all channels) */
void
ToneCurve::toLinear( const Vector &nonlinear, Vector &linear )
{
  errorCond( linear.getSize()== nonlinear.getSize(),
	     "  color space dimensions don't match!" );
  
  for( int i= nonlinear.getSize() ; i-->= 0 ; )
    linear[i]= toLinear( nonlinear[i] );
}


} /* namespace */

#endif /* COLOR_TONECURVE_C */

