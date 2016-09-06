// ==========================================================================
// $Id: ResamplerFactory.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef RESAMPLING_RESAMPLERFACTORY_C
#define RESAMPLING_RESAMPLERFACTORY_C

#include "MDA/Base/Errors.hh"
#include "ResamplerFactory.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** create a new resampler */
template<class T>
Resampler<T> *
ResamplerFactory<T>::create( ResampleMode mode )
{
  switch( mode )
  {
  case NearestNeighborResampling:
    return new NearestNeighborResampler<T>;
  case BoxResampling:
    return new BoxResampler<T>;
  case LinearResampling:
    return new LinearResampler<T>;
  case CubicResampling:
    return new CubicResampler<T>;
  case GaussResampling:
    return new GaussResampler<T>;
  case MinResampling:
    return new MinResampler<T>;
  case MaxResampling:
    return new MaxResampler<T>;
  default:
    warning( "  unsupported resampling mode!" );
  };

  // if we get here, things have gone wrong...
  return NULL;
}


//
// template instantiation code
//
template class ResamplerFactory<float>;
template class ResamplerFactory<double>;


} /* namespace */

#endif /* RESAMPLING_RESAMPLERFACTORY_C */

