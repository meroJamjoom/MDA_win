// ==========================================================================
// $Id: PaperSize.hh 661 2010-03-22 01:19:14Z heidrich $
// A collection of common paper sizes
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIB_PAPERSIZE_H
#define CAMERACALIB_PAPERSIZE_H

/*! \file  PaperSize.hh
    \brief A collection of common paper sizes
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

namespace MDA {
  
  /** \enum  PaperFormat PaperSize.hh
   */
  enum PaperFormat {
    A3,
    A4,
    A5,
    Letter,
    Legal
  };
  
  /** \class PaperSize PaperSize.hh
      simple representation of paper sizes
   */
  class PaperSize {
  public:
    /** constructor for standard paper size */
    inline PaperSize( PaperFormat format= Letter )
      : width( dimensions[format][0] ), height( dimensions[format][1] )
    {}
    /** read-only access to width */
    const double width;
    /** read-only access to height */
    const double height;
  private:
    /** repository of standard paper dimensions */
    static const double dimensions[][2];
  };

} /* namespace */

#endif /* CAMERACALIB_PAPERSIZE_H */

