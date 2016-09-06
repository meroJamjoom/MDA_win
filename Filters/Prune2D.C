// ==========================================================================
// $Id: Prune2D.C 374 2009-09-21 21:13:13Z heidrich $
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

#ifndef FILTERS_PRUNE2D_C
#define FILTERS_PRUNE2D_C

#include "Prune2D.hh"

namespace MDA {

/** structuring element 1 */
template<class T>
bool Prune2D<T>::se1[9]={false,true ,false,false,true ,false,false,false,false};
  
/** mask 1 */
template<class T>
bool Prune2D<T>::ma1[9]={false,true ,false,true ,true ,true ,true ,true ,true };
  
/** structuring element 2 */
template<class T>
bool Prune2D<T>::se2[9]={true ,false,false,false,true ,false,false,false,false};
  
/** mask 2 */
template<class T>
bool Prune2D<T>::ma2[9]={true ,true ,true ,true ,true ,true ,true ,true ,true };



// explicit template instation code
template class Prune2D<float>;
template class Prune2D<double>;
template class Endpoints2D<float>;
template class Endpoints2D<double>;

} /* namespace */

#endif /* FILTERS_PRUNE2D_C */

