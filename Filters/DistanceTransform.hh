// ==========================================================================
// $Id: DistanceTransform.hh 382 2009-09-25 05:16:27Z heidrich $
// distance transform for binary images
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

#ifndef FILTERS_DISTANCETRANSFORM_H
#define FILTERS_DISTANCETRANSFORM_H

/*! \file  DistanceTransform.hh
    \brief distance transform for binary images
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"

namespace MDA {

  /** \class DistanceTransform DistanceTransform.hh
      distance transform for binary images */
  template<class T>
  class DistanceTransform: public Filter<T> {

  public:
    
    /** default constructor */
    inline DistanceTransform( double p= 2.0 )
      : metric( p )
    {}
    
    /** apply flood fill to a number of dimensions and channels
        \bug the function does support operations on a subset of axes,
        but right now the axis 0 always needs to be part of the axis
        list */
#if defined (_WIN32) || defined (_WIN64)
	bool apply( Array<T> &a, BoundaryMethod boundary,
		ChannelList &channels, AxisList &axes );
#else
     virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes );
#endif
    
   protected:
    
    /** p, specifying an LP distance metric */
    double metric;
  };


} /* namespace */

#endif /* FILTERS_DISTANCETRANSFORM_H */

