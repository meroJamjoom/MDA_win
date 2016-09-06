// ==========================================================================
// $Id:$
// CIELUV color space (a.k.a. L*u*v*)
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

#ifndef COLOR_CIELUVSPACE_H
#define COLOR_CIELUVSPACE_H

/*! \file  CIELUVSpace.hh
    \brief CIELUV color space (a.k.a. L*u*v*)
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "YuvSpace.hh"
#include "Whitepoint.hh"

namespace MDA {

  /** \class CIELUVSpace CIELUVSpace.hh
   *  CIELUV color space (a.k.a. L*u*v*)
   *
   *  computation of L* uses INTENDED constants rather than the
   *  constants in the CIE stamndard, which would result in a
   *  discontinuity where the two segments meet (see
   *  http://www.brucelindbloom.com/index.html?LContinuity.html)
  */
  
  class CIELUVSpace: public ColorSpace {

  public:

    /** default constructor (uses D65 whitepoint) */
    inline CIELUVSpace()
      : wpYuv( 3 )
    {
      // D65 whitepoint with Y in [0..1]
      Vector XYZ( Whitepoint( "D65", 1.0 ).getXYZ() );
      yuvSpace.fromXYZ( XYZ, wpYuv );
    }
    
    /** constructor with specific whitepoint */
    inline CIELUVSpace( const Whitepoint &wp )
      : wpYuv( 3 )
    {
      Vector XYZ( wp.getXYZ() );
      yuvSpace.fromXYZ( wp.getXYZ(), wpYuv );
    }
    
    /** report dimensionality of space */
    inline virtual unsigned getDimension() const
    {
      return 3;
    }
    
    /** convert from L*u*v* to XYZ */
    inline virtual void toXYZ( const Vector &Luv, Vector &XYZ ) const
    { 
      // Y component
      if( Luv[0]> 216.0/27.0 )
	XYZ[1]= pow( (Luv[0]+16.0) / 116.0, 3.0 ) * wpYuv[0];
      else
	XYZ[1]= Luv[0] * 27.0 / 24389.0 * wpYuv[0];
      
      // coefficients for X, Z
      double a= (52.0*Luv[0] / (Luv[1] + 13.0*Luv[0]*wpYuv[1]) - 1.0) / 3.0;
      double b= -5.0 * XYZ[1];
      double c= -1.0 / 3.0;
      double d= XYZ[1] * (39.0*Luv[0] / (Luv[2] + 13.0*Luv[0]*wpYuv[2]) - 5.0);
      
      XYZ[0]= (d - b) / (a - c);
      XYZ[2]= XYZ[0] * a + b;
    }
    
    /** convert from XYZ to this space */
    inline virtual void fromXYZ( const Vector &XYZ, Vector &Luv ) const
    {
      Vector Yuv( 3 );
      yuvSpace.fromXYZ( XYZ, Yuv ); // Yu'v' representation of the color
      
      double yRatio= Yuv[0] / wpYuv[0];
      if( yRatio> 216.0/24389.0 )
	Luv[0]= 116.0 * pow( yRatio, 1.0/3.0 ) - 16.0;
      else
	Luv[0]= 24389.0 / 27.0 * yRatio;
      
      // u*, v* are straightforward
      Luv[1]= 13.0 * Luv[0] * (Yuv[1] - wpYuv[1]);
      Luv[2]= 13.0 * Luv[0] * (Yuv[2] - wpYuv[2]);
    }    
    
  protected:
  
    /** Yuv color space object */
    YuvSpace yuvSpace;
    
    /** Yu'v' representation of the whitepoint */
    Vector wpYuv;
    
  };


} /* namespace */

#endif /* COLOR_CIELUVSPACE_H */

