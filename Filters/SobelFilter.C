// ==========================================================================
// $Id: SobelFilter.C 248 2008-11-25 05:45:35Z heidrich $
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

#ifndef FILTERS_SOBELFILTER_C
#define FILTERS_SOBELFILTER_C

#include "SobelFilter.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions

// a smoothing kernel
template<class T>
HatFilter1D<T> SobelFilter<T>::smoothen= HatFilter1D<T>( 1, 2.0 );
    
// derivative kernel
template<class T>
FirstDerivative1D<T> SobelFilter<T>::derivative;


/** apply the filter to a number of dimensions and channels */
template<class T>
bool
SobelFilter<T>::apply( Array<T> &a, BoundaryMethod boundary,
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
  
  // this method performs redundant smoothings
  // for high dimensional datasets
  for( i= 0 ; i< channels.vec.size() ; i++ )
  {
    // initialize the accumuulator for this channel
    for( j= 0 ; j< numEntries ; j++ )
      (*a[accumulator])[j]= 0.0;
    
    for( j= 0 ; j< axes.vec.size() ; j++ )
    {
      // first, compute derivative along requested axis (parallel)
      derivative.apply( a, boundary, axes.vec[j],
			channels.vec[i], partialDeriv );
      
      // then, smoothen the result along all dimensions BUT the
      // requested one (parallel)
      for( k= 0 ; k< dim.vec.size() ; k++ )
	if( k!= axes.vec[j] )
	  smoothen.apply( a, boundary, k, partialDeriv, partialDeriv );
      
      // add square of result to accumulator (sequential)
      for( k= 0 ; k< numEntries ; k++ )
	(*a[accumulator])[k]+= (*a[partialDeriv])[k] * (*a[partialDeriv])[k];
    }
    
    // swap original channel with new data
    a.swapChannels( i, accumulator );
  }    
  
  // clean up temp channel
  a.deleteChannel( partialDeriv );
  a.deleteChannel( accumulator );  
  return status;
}
    

/* we use explicit instantiation for this class */
template class SobelFilter<float>;
template class SobelFilter<double>;

} /* namespace */

#endif /* FILTERS_SOBELFILTER_C */

