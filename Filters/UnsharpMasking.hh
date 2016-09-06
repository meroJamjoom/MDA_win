// ==========================================================================
// $Id: UnsharpMasking.hh 249 2008-11-25 06:07:08Z heidrich $
// an unsharp-masking filter
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

#ifndef FILTERS_UNSHARPENMASK_H
#define FILTERS_UNSHARPENMASK_H

/*! \file  UnsharpMasking.hh
    \brief an unsharpen-mask filter
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"


namespace MDA {

  /** \class UnsharpMasking UnsharpMasking.hh
      an unsharpen-mask filter */
  template<class T>
  class UnsharpMasking: public Filter<T> {

  public:
    
    /** constructor */
    UnsharpMasking( double _sigma= 1 )
      : sigma( _sigma )
    {}
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
    
  protected:

    /** sigma of the Gaussian to be used */
    double sigma;
  };


} /* namespace */


#endif /* FILTERS_UNSHARPENMASK_H */

