// ==========================================================================
// $Id: ExtremaDetector.hh 319 2009-05-27 21:17:01Z heidrich $
// Finds local maxima nand minima in the array, and labels them 1 and
// 0, respectively. All other pixels get a value of 0.5.
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

#ifndef FILTERS_EXTREMADETECTOR_H
#define FILTERS_EXTREMADETECTOR_H

/*! \file  ExtremaDetector.hh
    \brief Finds local maxima nand minima in the array, and labels
    them 1 and 0, respectively. All other pixels get a value of 0.5.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"

namespace MDA {

  /** \class ExtremaDetector ExtremaDetector.hh
      Finds local maxima nand minima in the array, and labels them 1
      and 0, respectively. All other pixels get a value of 0.5. */
  template<class T>
  class ExtremaDetector: public Filter<T> {

  public:

    /** default constructor */
    ExtremaDetector( unsigned _radius= 1 )
      : radius( _radius )
    {}
    
    /** apply filter */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
    
  protected:
    
    /** filter radius */
    unsigned radius;
    
  };


} /* namespace */

#endif /* FILTERS_EXTREMADETECTOR_H */

