// ==========================================================================
// $Id: SobelFilter.hh 248 2008-11-25 05:45:35Z heidrich $
// Sobel edge detection filter
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

#ifndef FILTERS_SOBELFILTER_H
#define FILTERS_SOBELFILTER_H

/*! \file  SobelFilter.hh
    \brief Sobel edge detection filter
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"
#include "Linear1DFilter.hh"

namespace MDA {

  /** \class SobelFilter SobelFilter.hh
      Sobel edge detection filter */
  template<class T>
  class SobelFilter: public Filter<T> {

  public:

    /** apply the filter to a number of dimensions and channels */
    bool apply( Array<T> &a, BoundaryMethod boundary,
                ChannelList &channels, AxisList &axes );
    
  protected:
    
    /** a smoothing kernel */
    static HatFilter1D<T> smoothen;
    
    /** derivative kernel */
    static FirstDerivative1D<T> derivative;
  };


} /* namespace */


#endif /* FILTERS_SOBELFILTER_H */

