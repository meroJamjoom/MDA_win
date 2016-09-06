// ==========================================================================
// $Id: SeparableFilter.C 248 2008-11-25 05:45:35Z heidrich $
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

#ifndef FILTERS_SEPARABLEFILTER_C
#define FILTERS_SEPARABLEFILTER_C

#include "SeparableFilter.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions


/** apply the filter to a number of dimensions and channels */
template<class T>
bool
SeparableFilter<T>::apply( Array<T> &a, BoundaryMethod boundary,
			   ChannelList &channels, AxisList &axes )
{
  bool status= true;
  for( unsigned i= 0 ; i< channels.vec.size() ; i++ )
    for( unsigned j= 0 ; j< axes.vec.size() ; j++ )
      status&= filter->apply( a, boundary, axes.vec[j],
			      channels.vec[i], channels.vec[i] );
  return status;
}


/* we use explicit instantiation for this class */
template class SeparableFilter<float>;
template class SeparableFilter<double>;


} /* namespace */

#endif /* FILTERS_SEPARABLEFILTER_C */

