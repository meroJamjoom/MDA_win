// ==========================================================================
// $Id: Filter.hh 549 2010-01-08 08:42:25Z felixhei $
// An n-dimensional filter for any number ofchannels
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef FILTERS_FILTER_H
#define FILTERS_FILTER_H

/*! \file  Filter.hh
    \brief An n-dimensional filter for any number of channels
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/Array.hh"
#include "MDA/Array/Boundary.hh"


namespace MDA {
  
  /** \class Filter Filter.hh
      An n-dimensional filter for any number of channels */
  template<class T>
  class Filter {
    
  public:

    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )= 0;
    
    /** destructor */
    virtual ~Filter() {}
  };


} /* namespace */

#endif /* FILTERS_FILTER_H */

