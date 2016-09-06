// ==========================================================================
// $Id:$
// CIELAB color space (a.k.a. L*a*b*)
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

#ifndef COLOR_CIELABSPACE_H
#define COLOR_CIELABSPACE_H

/*! \file  CIELABSpace.hh
    \brief CIELAB color space (a.k.a. L*a*b*)
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "YuvSpace.hh"
#include "Whitepoint.hh"

namespace MDA {

  /** \class CIELABSpace CIELABSpace.hh
   *  CIELAB color space (a.k.a. L*a*b*)
   *
   *  computation of L* uses INTENDED constants rather than the
   *  constants in the CIE stamndard, which would result in a
   *  discontinuity where the two segments meet
   *  (see http://www.brucelindbloom.com/index.html?LContinuity.html)
  */
  
  class CIELABSpace: public ColorSpace {

  public:
    
    /** default constructor (uses D65 whitepoint) */
    inline CIELABSpace()
      : wpXYZ( 3 )
    {
      // D65 whitepoint with Y= 1.0
      wpXYZ.copy( Whitepoint( "D65", 1.0 ).getXYZ() );
    }
    
    /** constructor with specific whitepoint */
    inline CIELABSpace( const Whitepoint &wp )
      : wpXYZ( 3 )
    {
      wpXYZ.copy( wp.getXYZ() );
    }
     
    /** report dimensionality of space */
    inline virtual unsigned getDimension() const
    {
      return 3;
    }
    
    /** convert from L*a*b* to XYZ */
    inline virtual void toXYZ( const Vector &Lab, Vector &XYZ ) const
    {
      double fX= Lab[0]/116.0 + Lab[1]/500.0 + 4.0/29.0;
      double fY= Lab[0]/116.0 + 4.0/29.0;
      double fZ= Lab[0]/116.0 - Lab[2]/200.0 + 4.0/29.0;
      
      XYZ[0]= fInv( fX ) * wpXYZ[0];
      XYZ[1]= fInv( fY ) * wpXYZ[1];
      XYZ[2]= fInv( fZ ) * wpXYZ[2];
    }
    
    /** convert from XYZ to L*a*b* */
    inline virtual void fromXYZ( const Vector &XYZ, Vector &Lab ) const
    {
      // non-linear value compression
      double fX= f( XYZ[0] / wpXYZ[0] );
      double fY= f( XYZ[1] / wpXYZ[1] );
      double fZ= f( XYZ[2] / wpXYZ[2] );
      
      Lab[0]= 116.0 * fY - 16.0;
      Lab[1]= 500.0 * (fX - fY);
      Lab[2]= 200.0 * (fY - fZ);
    }
    
    
  protected:

    /** compressive function */
    inline double f( double val ) const
    {
      if( val> 216.0/24389.0 )
	return pow( val, 1.0/3.0 );
      else
	return 24389.0/27.0/116.0 * val + 16.0/116.0;
    }
    
    /** inverse of compressive function */
    inline double fInv( double val ) const
    {
      if( val> 27.0/216.0 )
	return pow( val, 3.0 );
      else
	return (val - 116.0/16.0) * 27.0*116.0/24389.0;
    }
    
    /** Yu'v' representation of the whitepoint */
    Vector wpXYZ;

  };


} /* namespace */

#endif /* COLOR_CIELABSPACE_H */

