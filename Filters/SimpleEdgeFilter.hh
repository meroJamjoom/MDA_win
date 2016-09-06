// ==========================================================================
// $Id: SimpleEdgeFilter.hh 248 2008-11-25 05:45:35Z heidrich $
// A simple (divided difference) gradient based edge filter
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

#ifndef FILTERS_SIMPLEEDGEFILTER_H
#define FILTERS_SIMPLEEDGEFILTER_H

/*! \file  SimpleEdgeFilter.hh
    \brief A simple (divided difference) gradient based edge filter
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"
#include "Linear1DFilter.hh"

namespace MDA {

  /** \class SimpleEdgeFilter SimpleEdgeFilter.hh
      A simple (divided difference) gradient based edge filter */
  template<class T>
  class SimpleEdgeFilter: public Filter<T> {

  public:

    /** default constructor */
    SimpleEdgeFilter( int baseline= 1 )
      : derivative( FirstDerivative1D<T>( baseline ) )
    {}

    /** apply the filter to a number of dimensions and channels */
    bool apply( Array<T> &a, BoundaryMethod boundary,
                ChannelList &channels, AxisList &axes );
    
  protected:
    
    /** derivative kernel */
    FirstDerivative1D<T> derivative;
    
  };


} /* namespace */

//
// this class is using explicit instantiation, so the folowing does not apply
//
// include templates in header file (GNU compilers)
//#if defined(__GNUG__) && !defined(FILTERS_SIMPLEEDGEFILTER_C)
//#define FILTERS_SIMPLEEDGEFILTER_TEMPLATES
//#include "SimpleEdgeFilter.C"
//#endif


#endif /* FILTERS_SIMPLEEDGEFILTER_H */

