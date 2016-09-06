// ==========================================================================
// $Id:$
// abstract base class for all color spaces
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

#ifndef COLOR_COLORSPACE_H
#define COLOR_COLORSPACE_H

/*! \file  ColorSpace.hh
    \brief abstract base class for all color spaces
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/LinearAlgebra/Vector.hh"

namespace MDA {

  /** \class ColorSpace ColorSpace.hh
   *   abstract base class for all color spaces
   * 
   *   all color spaces need to be able to convert to/from XYZ  
   */
  
  class ColorSpace {

  public:

    /** destructor */
    virtual ~ColorSpace() {}
    
    /** report dimensionality of space */
    virtual unsigned getDimension() const = 0;
    
    /** convert from this space to XYZ */
    virtual void toXYZ( const Vector &tristimulus, Vector &XYZ ) const= 0;
    
    /** convert from XYZ to this space */
    virtual void fromXYZ( const Vector &XYZ, Vector &tristimulus ) const= 0;
    
  };


} /* namespace */

#endif /* COLOR_COLORSPACE_H */

