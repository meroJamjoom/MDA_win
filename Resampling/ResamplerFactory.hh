// ==========================================================================
// $Id: ResamplerFactory.hh 319 2009-05-27 21:17:01Z heidrich $
// factory for resamplers
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

#ifndef RESAMPLING_RESAMPLERFACTORY_H
#define RESAMPLING_RESAMPLERFACTORY_H

/*! \file  ResamplerFactory.hh
    \brief factory for resamplers
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "ResamplerOption.hh"

namespace MDA {

  /** \class ResamplerFactory ResamplerFactory.hh
      factory for resamplers */
  
  template<class T>
  class ResamplerFactory {

  public:
    
    /** default constructor */
    ResamplerFactory()
    {}
    
    /** create a new resampler */
    Resampler<T> *create( ResampleMode mode );

  };


} /* namespace */



#endif /* RESAMPLING_RESAMPLERFACTORY_H */

