// ==========================================================================
// $Id: DistanceTransform.C 383 2009-09-25 05:55:26Z heidrich $
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

#ifndef FILTERS_DISTANCETRANSFORM_C
#define FILTERS_DISTANCETRANSFORM_C

#include <MDA/Array/Neighborhood.hh>
#include <MDA/Array/PaddedChannelTraverser.hh>
#include <MDA/LinearAlgebra/Vector.hh>
#include "DistanceTransform.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** apply filter to selected channels */
template<class T>
bool
DistanceTransform<T>::apply(Array<T> &a, BoundaryMethod boundary,
ChannelList &channels, AxisList &axes)
{
	unsigned long i, j;

	// determine basic dimensions
	CoordinateVector dim = a.getDimension();
	unsigned dimension = dim.vec.size();
	unsigned long numPixels = size(dim);
	CoordinateVector paddedDim = pad(dim, 1);
	unsigned long numPaddedPixels = size(paddedDim);

	// define a 1-neighborhood according to the metric, and a traverser
	// for the padded and unpadded grid
	Neighborhood n(dim, 1, axes, 1, metric);
	PaddedChannelTraverser traverser(dim, 1);

	// data pointers for padded and unpadded channels
	T *uData;

	// allocate Vector array of closest points and set to point far away
#if defined (_WIN32) || defined (_WIN64)
	Vector far_v( dimension, 1e10 );
#else
  Vector far( dimension, 1e10 );
#endif
  Vector *closestPoints= new Vector[numPaddedPixels];
  for( i= 0 ; i< numPaddedPixels ; i++ )
#if defined (_WIN32) || defined (_WIN64)
	  closestPoints[i]= far_v;
#else
    closestPoints[i]= far;
#endif
  // process all channels independently
  for( j= 0 ; j< channels.vec.size() ; j++ )
  {
    // get pointer to raw data (unpadded)
    uData= &((*a[channels.vec[j]])[0]);
    
    //
    // 1st pass: causal neighborhood
    // (top-left to bottom-right using first half of pixels in neighborhood)
    //
    Vector current( dimension, 0.0 );
    for( traverser.begin() ; !traverser.isAtEnd() ; ++traverser )
    {
      unsigned long uOff= traverser.uOffset();
      unsigned long pOff= traverser.pOffset();
      double currDist, d;
      
      if( uData[uOff]> 0.0 )
      {
	// case 1: if current position has an object pixel, the closest
	// point is the pixel itself
	closestPoints[pOff].copy( current );
	uData[uOff]= 0.0;
      }
      else
      {
	// case 2: closest point is also the closest point for one of
	// the previously visited pixels (i.e. one of the pixels in
	// the first half of the 1 neighborhood)
	currDist= dist( closestPoints[pOff], current, metric );
	for( i= n.getNumPixels()/2, n.begin() ; i> 0 ; i--, ++n )
	{
	  d= dist( closestPoints[pOff+(*n).pOff], current, metric );
	  if( d< currDist )
	  {
	    closestPoints[pOff]= closestPoints[pOff+(*n).pOff];
	    currDist= d;
	  }
	}
	uData[uOff]= currDist;
      }
      
      // advance current position
      for( i= 0 ; i< dimension ; i++ )
	if( (current[i]+= 1.0)>= dim.vec[i] )
	  current[i]= 0.0;
	else
	  break;
    }
    
    //
    // 2nd pass: anti-causal neighborhood
    // (bottom-right to top-left using second half of pixels in neighborhood)
    //
    for( i= 0 ; i< dimension ; i++ )
      current[i]= dim.vec[i]-1;
    for( traverser.begin() ; !traverser.isAtEnd() ; ++traverser )
    {
      unsigned long uOff= numPixels-1 - traverser.uOffset();
      unsigned long pOff= numPaddedPixels-1 - traverser.pOffset();
      double currDist, d;
      
      // only one case now - object pixels are already initialized
      currDist= uData[uOff];
      for( n.begin( n.getNumPixels()/2+1 ) ; !n.isAtEnd() ; ++n )
      {
	d= dist( closestPoints[pOff+(*n).pOff], current, metric );
	if( d< currDist )
	{
	  closestPoints[pOff]= closestPoints[pOff+(*n).pOff];
	  currDist= d;
	}
      }
      uData[uOff]= currDist;
      
      // advance current position
      for( i= 0 ; i< dimension ; i++ )
	if( (current[i]-= 1.0)< 0.0 )
	  current[i]= dim.vec[i]-1;
	else
	  break;
    }
  }
  
  //
  // clean up
  //
  delete [] closestPoints;
  
  return true;
}


// explicit template instation code

template class DistanceTransform<float>;
template class DistanceTransform<double>;



} /* namespace */

#endif /* FILTERS_DISTANCETRANSFORM_C */

