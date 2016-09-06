// ==========================================================================
// $Id:$
// Simple, von Kries-style whitebalancing
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

#ifndef COLOR_WHITEBALANCE_H
#define COLOR_WHITEBALANCE_H

/*! \file  WhiteBalance.hh
    \brief Simple, von Kries-style whitebalancing
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "LinearTristimulusSpace.hh"
#include "Whitepoint.hh"

namespace MDA {

  /** \class WhiteBalance WhiteBalance.hh
      Simple, von Kries-style whitebalancing */
  
  class WhiteBalance {

  public:

    /** constructor from named color spaces */
    WhiteBalance( const Whitepoint &srcWP, const Whitepoint &dstWP,
		  const char *srcSpace= "sRGB",
                  const char *vonKriesSpace= "HPE" );
    
    /** constructor from externally created color spaces */
    WhiteBalance( const Whitepoint &srcWP, const Whitepoint &dstWP,
		  ColorSpace *srcSpace, ColorSpace *vonKriesSpace );
    
    /** destructor */
    inline ~WhiteBalance()
    {
      delete src;
      delete vonKries;
    }
    
    /** apply whitebalance (result overwrites input) */
    void apply( Vector &color );
    
  protected:
    
    /** source color space */
    ColorSpace *src;
    
    /** von Kries space */
    ColorSpace *vonKries;
    
    /** multipliers in von Kries space */
    double mult[3];
    
  };


} /* namespace */

#endif /* COLOR_WHITEBALANCE_H */

