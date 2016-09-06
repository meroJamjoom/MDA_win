// ==========================================================================
// $Id: ImprovedHarrisCorner.C 259 2008-12-30 09:15:14Z heidrich $
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

#ifndef FILTERS_IMPROVEDHARRISCORNER_C
#define FILTERS_IMPROVEDHARRISCORNER_C

#include <math.h>

#include "MorphologicalOps.hh"
#include "Linear1DFilter.hh"
#include "ImprovedHarrisCorner.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** apply the filter to a number of dimensions and channels */
template <class T>
bool
ImprovedHarrisCorner<T>::apply( Array<T> &a, BoundaryMethod boundary,
				ChannelList &channels, AxisList &axes )
{
  unsigned long i, k;
  
  // for now, we can only handle 2D corners
  errorCond( axes.vec.size()== 2, "  only 2D corners supported so far" );
  
  // basic dimension stats
  CoordinateVector dim= a.getDimension();
  unsigned dimension= dim.vec.size();
  unsigned long numElements= 1;
  for( i= 0 ; i< dimension ; i++ )
    numElements*= dim.vec[i];
  
  // Gauss smoothing of afected channels (in ALL dimensions)
  Gaussian1D<T> gauss( sigma );
  for( k= 0 ; k< channels.vec.size() ; k++ )
    for( i= 0 ; i< dimension ; i++ )
      gauss.apply( a, boundary, i, channels.vec[k], channels.vec[k] );
  
  // we need two new channels for intermediate results
  unsigned tmpCh1= a.addChannel();
  unsigned tmpCh2= a.addChannel();
  
  // compute second derivatives and Harris formula for each channel
  FirstDerivative1D<T> firstDeriv;
  SecondDerivative1D<T> secondDeriv( 2 );
  Dilate1D<T> dilate;
  for( k= 0 ; k< channels.vec.size() ; k++ )
  {
    unsigned ch= channels.vec[k];
    // temp channels 1 and 2 are Ixx, Iyy (second deriv in x, y)
    secondDeriv.apply( a, boundary, axes.vec[0], ch, tmpCh1 );
    secondDeriv.apply( a, boundary, axes.vec[1], ch, tmpCh2 );
    // original channel is mixed derivative Ixy
    firstDeriv.apply( a, boundary, axes.vec[0], ch, ch );
    firstDeriv.apply( a, boundary, axes.vec[1], ch, ch );
    
    // set original channel to det(M) - traceWeight *trace(M)^2
    // where M is Harris' matrix:
    // (Ixx Ixy)
    // (Ixy Iyy)
    T* tmp1= &((*a[tmpCh1])[0]);
    T* tmp2= &((*a[tmpCh2])[0]);
    T* tmp3= &((*a[ch])[0]);
    for( i= 0 ; i< numElements ; i++ )
    {
      double a= tmp1[i];
      double b= tmp2[i];
      double c= tmp3[i];
      tmp3[i]= a*b-c*c - traceWeight*(a+b)*(a+b);
    }
    
    // compute the max of the Harris value over a radius of width 3*sigma
    dilate.apply( a, boundary, axes.vec[0], ch, tmpCh1 );
    dilate.apply( a, boundary, axes.vec[1], tmpCh1, tmpCh1 );
    
    for( i= 0 ; i< numElements ; i++ )
      if( tmp1[i]< cornerThreshold )
	tmp3[i]= 0.0;
      else
	tmp3[i]= (tmp1[i] - tmp3[i]) / tmp1[i];
  } 
  
  // remove temp channels
  a.deleteChannel( tmpCh2 );
  a.deleteChannel( tmpCh1 );
  
  return true;
}


// explicit template instation code

template class ImprovedHarrisCorner<float>;
template class ImprovedHarrisCorner<double>;



} /* namespace */

#endif /* FILTERS_IMPROVEDHARRISCORNER_C */

