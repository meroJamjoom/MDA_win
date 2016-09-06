// ==========================================================================
// $Id: ExtremaDetector.C 319 2009-05-27 21:17:01Z heidrich $
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

#ifndef FILTERS_EXTREMADETECTOR_C
#define FILTERS_EXTREMADETECTOR_C

#include "MorphologicalOps.hh"
#include "ExtremaDetector.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** apply the filter to a number of dimensions and channels */
template <class T>
bool
ExtremaDetector<T>::apply( Array<T> &a, BoundaryMethod boundary,
			   ChannelList &channels, AxisList &axes )
{
  unsigned long i, k;
  
  // basic dimension stats
  CoordinateVector dim= a.getDimension();
  unsigned dimension= dim.vec.size();
  unsigned long numElements= 1;
  for( i= 0 ; i< dimension ; i++ )
    numElements*= dim.vec[i];
  
  // two temp channels for min, max
  unsigned minCh= a.addChannel();
  unsigned maxCh= a.addChannel();
  T* minData= &((*a[minCh])[0]);
  T* maxData= &((*a[maxCh])[0]);
  
  // apply filter to each channel
  Erode1D<T> erode( radius );
  Dilate1D<T> dilate( radius );
  for( k= 0 ; k< channels.vec.size() ; k++ )
  {
    // first compute min/max through erosion/dilation along requested axes
    erode.apply( a, boundary, axes.vec[0], channels.vec[k], minCh );
    dilate.apply( a, boundary, axes.vec[0], channels.vec[k], maxCh );
    for( i= 1 ; i< axes.vec.size() ; i++ )
    {
      erode.apply( a, boundary, axes.vec[i], minCh, minCh );
      dilate.apply( a, boundary, axes.vec[i], maxCh, maxCh );
    }
    
    // now find and label extrema
    T* dstData= &((*a[channels.vec[k]])[0]);
    for( i= 0 ; i< numElements ; i++ )
      if( dstData[i]== minData[i] )
	if( dstData[i]== maxData[i] )
	  dstData[i]= 0.5;
	else
	  dstData[i]= 0.0;
      else
	if( dstData[i]== maxData[i] )
	  dstData[i]= 1.0;
	else
	  dstData[i]= 0.5;
  }
  
  // clean up
  a.deleteChannel( maxCh );
  a.deleteChannel( minCh );
  
  return true;
}


// explicit template instation code

template class ExtremaDetector<float>;
template class ExtremaDetector<double>;



} /* namespace */

#endif /* FILTERS_EXTREMADETECTOR_C */

