// ==========================================================================
// $Id:$
// Baseclass for tone curves
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

#ifndef COLOR_TONECURVE_H
#define COLOR_TONECURVE_H

/*! \file  ToneCurve.hh
    \brief Baseclass for tone curves
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif


#include "MDA/Base/Errors.hh"
#include "MDA/LinearAlgebra/Vector.hh"

namespace MDA {

  /** \class ToneCurve ToneCurve.hh
      Baseclass for tone curves */
  
  class ToneCurve {

  public:
    
    /** destructor */
    virtual ~ToneCurve() {}
    
    /** conversion of scalar from linear to non-linear space */
    virtual double fromLinear( double linearVal )= 0;
    
    /** conversion of scalar from non-linear to linear space */
    virtual double toLinear( double nonlinearVal )= 0;
    
    /** conversion of color vector from linear to non-linear space
	(default implementation uses identical curves for all channels) */
    virtual void fromLinear( const Vector &linear, Vector &nonlinear );
    
    /** conversion of color vector from non-linear to linear space
	(default implementation uses identical curves for all channels) */
    virtual void toLinear( const Vector &nonlinear, Vector &linear );
    
  };


} /* namespace */

#endif /* COLOR_TONECURVE_H */

