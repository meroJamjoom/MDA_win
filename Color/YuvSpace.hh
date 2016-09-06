// ==========================================================================
// $Id:$
// Yu'v' color space
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

#ifndef COLOR_YUVSPACE_H
#define COLOR_YUVSPACE_H

/*! \file  YuvSpace.hh
    \brief Yu'v' color space
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "ColorSpace.hh"

namespace MDA {

  /** \class YuvSpace YuvSpace.hh
      Yu'v' color space */
  
  class YuvSpace: public ColorSpace {

  public:

    /** report dimensionality of space */
    inline virtual unsigned getDimension() const
    {
      return 3;
    }
    
    /** convert from Yu'v' to XYZ */
    inline virtual void toXYZ( const Vector &Yuv, Vector &XYZ ) const
    {
      double mult= 0.25 * Yuv[0] / Yuv[2];
      
      XYZ[1]= Yuv[0]; // Y
      XYZ[0]= mult * 9.0 * Yuv[1]; // X
      XYZ[2]= mult * (12.0 - 3.0*Yuv[1] - 20.0*Yuv[2]); // Z
    }
    
    /** convert from XYZ to this space */
    inline virtual void fromXYZ( const Vector &XYZ, Vector &Yuv ) const
    {
      double denominator= XYZ[0] + 15.0*XYZ[1] + 3.0*XYZ[2];
      
      Yuv[0]= XYZ[1]; // Y
      Yuv[1]= 4.0*XYZ[0] / denominator; // u'
      Yuv[2]= 9.0*XYZ[1] / denominator; // v'
    }    

  };


} /* namespace */

#endif /* COLOR_YUVSPACE_H */

