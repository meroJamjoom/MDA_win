// ==========================================================================
// $Id: UnsharpMasking.C 249 2008-11-25 06:07:08Z heidrich $
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

#ifndef FILTERS_UNSHARPENMASK_C
#define FILTERS_UNSHARPENMASK_C

#include "Linear1DFilter.hh"
#include "UnsharpMasking.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions


/** apply the filter to a number of dimensions and channels */
template<class T>
bool
UnsharpMasking<T>::apply( Array<T> &a, BoundaryMethod boundary,
			 ChannelList &channels, AxisList &axes )
{
  unsigned long i, j;
  bool status= true;
  
  // compute number of pixels in a channel
  unsigned long numPixels= 1;
  for( i= 0 ; i< a.getDimension().vec.size() ; i++ )
    numPixels*= a.getDimension().vec[i];
  
  // define the blur filter (use fast Gaussian for large radii, normal
  // Gaussian for small ones)
  Filter1D<T> *filter;
  if( sigma> 3.0 && boundary!= Renormalize )
    filter= new FastGaussian1D<T>( sigma );
  else
    filter= new Gaussian1D<T>( sigma );
  
  // a temporary channel for a blurred intermediate result
  int tmp= a.addChannel();
  
  // go over each channel, blur it, and then mask it out of the original
  for( i= 0 ; i< channels.vec.size() ; i++ )
  {
    // blur along the relevant axes
    // - blur along first axis also copies to tmp array
    // - this is done in parallel, since the blur filters are parallel
    status&= filter->apply( a, boundary, axes.vec[0], channels.vec[i], tmp );
    for( j= 1 ; j< axes.vec.size() ; j++ )
      status&= filter->apply( a, boundary, axes.vec[j], tmp, tmp );
    
    // now do the masking, per pixel
    // - this is not currently parallel; it is probably not worth the effort
    T* orig= &((*a[channels.vec[i]])[0]);
    T* blur= &((*a[tmp])[0]);
    for( j= 0 ; j< numPixels ; j++, orig++, blur++ )
      *orig+= *orig - *blur;
  }
  
  // clean up
  a.deleteChannel( tmp );
  delete filter;
  
  return status;
}
  
/* we use explicit instantiation for this class */
template class UnsharpMasking<float>;
template class UnsharpMasking<double>;

} /* namespace */

#endif /* FILTERS_UNSHARPENMASK_C */

