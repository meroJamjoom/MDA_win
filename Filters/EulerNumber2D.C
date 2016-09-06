// ==========================================================================
// $Id: EulerNumber2D.C 372 2009-09-20 18:57:38Z heidrich $
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

#ifndef FILTERS_EULERNUMBER2D_C
#define FILTERS_EULERNUMBER2D_C

#include "EulerNumber2D.hh"

namespace MDA {

/** table for 4-neighborhood */
template <class T>
const double EulerNumber2D<T>::fourNeighborTable[16]= {
  .0, .25, .25, .0, .25, .0, .5, -.25, .25, .5, .0, -.25, .0, -.25, -.25, .0
};


/** table for 8-neighborhood */
template <class T>
const double EulerNumber2D<T>::eightNeighborTable[16]= {
  .0, .25, .25, .0, .25, .0, -.5, -.25, .25, -.5, .0, -.25, .0, -.25, -.25, .0
};


// template instantiation
template class EulerNumber2D<float>;
template class EulerNumber2D<double>;



} /* namespace */

#endif /* FILTERS_EULERNUMBER2D_C */

