// ==========================================================================
// $Id: Rectify.hh $
// image rectification for stereo pairs
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: bradleyd ()
// Email:   bradleyd@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIB_RECTIFY_H
#define CAMERACALIB_RECTIFY_H

/*! \file  Rectify.hh
    \brief image rectification for stereo pairs
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

namespace MDA {
  
  /** compute the 3x3 rectifying transformations for a pair of cameras */
  void getRectifyTransforms( const Camera C1, const Camera C2, Matrix& Rect1, Matrix& Rect2 );
  
  // more functionality to come (for actually performing rectification using the transforms)

} /* namespace */

#endif /* CAMERACALIB_RECTIFY_H */

