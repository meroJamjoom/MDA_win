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

#ifndef COLOR_GAMMACORRECTEDSPACE_C
#define COLOR_GAMMACORRECTEDSPACE_C

#include "GammaCorrectedSpace.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** convert from this space to XYZ */
void
GammaCorrectedSpace::toXYZ( const Vector &tristimulus, Vector &XYZ ) const
{
  Vector linearTriStim( 3 );
  gamma->toLinear( tristimulus, linearTriStim );
  space->toXYZ( linearTriStim, XYZ );
}

/** convert from XYZ to this space */
void
GammaCorrectedSpace::fromXYZ( const Vector &XYZ, Vector &tristimulus ) const
{
  Vector linearTriStim( 3 );
  space->fromXYZ( XYZ, linearTriStim );
  gamma->fromLinear( linearTriStim, tristimulus );
}

} /* namespace */

#endif /* COLOR_GAMMACORRECTEDSPACE_C */

