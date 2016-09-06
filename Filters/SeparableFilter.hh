// ==========================================================================
// $Id: SeparableFilter.hh 248 2008-11-25 05:45:35Z heidrich $
// A separable filter composed of 1D filters
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

#ifndef FILTERS_SEPARABLEFILTER_H
#define FILTERS_SEPARABLEFILTER_H

/*! \file  SeparableFilter.hh
    \brief A separable filter composed of 1D filters
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter1D.hh"

namespace MDA {

  /** \class SeparableFilter SeparableFilter.hh
      A separable filter composed of 1D filters */
  template<class T>
  class SeparableFilter: public Filter<T> {

  public:
    
    /** default constructor */
    SeparableFilter( Filter1D<T> *_filter )
      : filter( _filter )
    {}
    
    /** destructur (cleans up the 1D filter as well!) */
    ~SeparableFilter()
    {
      delete filter;
    }

    /** apply the filter to a number of dimensions and channels */
    bool apply( Array<T> &a, BoundaryMethod boundary,
                ChannelList &channels, AxisList &axes );
    
  protected:

    /** the 1D filter */
    Filter1D<T> *filter;
    
  };


} /* namespace */

//
// this class is using explicit instantiation, so the folowing does not apply
//
// include templates in header file (GNU compilers)
//#if defined(__GNUG__) && !defined(FILTERS_SEPARABLEFILTER_C)
//#define FILTERS_SEPARABLEFILTER_TEMPLATES
//#include "SeparableFilter.C"
//#endif


#endif /* FILTERS_SEPARABLEFILTER_H */

