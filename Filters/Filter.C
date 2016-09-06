// ==========================================================================
// $Id: Filter.C 87 2007-05-22 12:56:39Z heidrich $
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

#ifndef FILTERS_FILTER_C
#define FILTERS_FILTER_C

#include "Filter.hh"

namespace MDA {


  /* we use explicit instantiation for this class */
  template class Filter<float>;
  template class Filter<double>;


} /* namespace */

#endif /* FILTERS_FILTER_C */

