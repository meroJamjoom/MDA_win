// ==========================================================================
// $Id: SimpleEdgeFilter.C 248 2008-11-25 05:45:35Z heidrich $
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

#ifndef FILTERS_SIMPLEEDGEFILTER_C
#define FILTERS_SIMPLEEDGEFILTER_C

#include "SimpleEdgeFilter.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions



/** apply the filter to a number of dimensions and channels */
template<class T>
bool
SimpleEdgeFilter<T>::apply( Array<T> &a, BoundaryMethod boundary,
			    ChannelList &channels, AxisList &axes )
{
  unsigned long i, j, k;
  
  // need two temporary channels:
  // one for accumulation
  unsigned accumulator= a.addChannel();
  // and one for a directional derivative
  unsigned partialDeriv= a.addChannel();
  
  bool status= true;
  CoordinateVector dim= a.getDimension();
  unsigned long numEntries= 1;
  for( int l= dim.vec.size()-1 ; l>= 0 ; l-- )
    numEntries*= dim.vec[l];
  
  for( i= 0 ; i< channels.vec.size() ; i++ )
  {
    // initialize the accumulator for this channel
    for( j= 0 ; j< numEntries ; j++ )
      (*a[accumulator])[j]= 0.0;
    
    for( j= 0 ; j< axes.vec.size() ; j++ )
    {
      // first, compute derivative along requested axis (parallel)
      derivative.apply( a, boundary, axes.vec[j],
			channels.vec[i], partialDeriv );
      
      // then add square of result to accumulator (sequential for now)
      for( k= 0 ; k< numEntries ; k++ )
	(*a[accumulator])[k]+= (*a[partialDeriv])[k] * (*a[partialDeriv])[k];
    }
    
    // simply swap the accumulator and the source channel
    a.swapChannels( i, accumulator );
  }    
  
  // clean up temp channels
  a.deleteChannel( partialDeriv );
  a.deleteChannel( accumulator );  
  return status;
}
    



/* we use explicit instantiation for this class */
template class SimpleEdgeFilter<float>;
template class SimpleEdgeFilter<double>;

} /* namespace */

#endif /* FILTERS_SIMPLEEDGEFILTER_C */

