// ==========================================================================
// $Id: MorphologicalOps.C 247 2008-11-25 05:23:41Z heidrich $
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

#ifndef FILTERS_MORPHOLOGICALOPS_C
#define FILTERS_MORPHOLOGICALOPS_C

#include "MDA/Array/Boundary.hh"
#include "MorphologicalOps.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions


/** apply erosion filter to one channel of an array
 **
 ** This code is mostly copied from Linear1DFilter<T>::apply, but with
 ** a min operator in the center loop.
 */
template<class T>
bool
Erode1D<T>::apply( Array<T> &array,
		   BoundaryMethod boundary,
		   unsigned int axis,
		   unsigned int inChannel,
		   unsigned int outChannel )
{
  // renormalize is identical to clamp for this filter
  if( boundary== Renormalize )
    boundary= Clamp;
  
  // call superclass method for actually applying the filter
  return Filter1D<T>::apply( array, boundary, axis, inChannel, outChannel );
}

/** apply filter to a single line in the array */
template <class T>
void
Erode1D<T>::apply( T *startPos, unsigned long incr,
		   unsigned long numElements, BoundaryMethod boundary,
		   T background, T *startPosOut )
{
  unsigned long j, k, l;
  
  // allocate temporary memory for one line
  T* lineBuf= new T[numElements+2*radius];
  
  // fetch line into temp buffer
  fetchLine<T>( lineBuf, startPos, incr, numElements, radius,
		boundary, background );
  
  // apply filter:
  // find minimum over first 2*radius+1 elements
  T minVal= lineBuf[0];
  for( j= 1 ; j< 2*radius+1 ; j++ )
    minVal= minVal<= lineBuf[j] ? minVal : lineBuf[j];
  
  // application to the full line
  for( j= 0, k= 2*radius+1 ; j< numElements ; j++, k++ )
  {
    startPosOut[j*incr]= minVal;
    // update minval
    if( j< numElements-1 )
      if( lineBuf[k]< minVal )
	// if the new value added to the sliding window is less
	  // than current min, then that is the new min
	minVal= lineBuf[k];
      else
	if( minVal== lineBuf[j] )
	{
	  // if the minimum fell out of the sliding window, we
	  // have to find it from scratch
	  minVal= lineBuf[j+1];
	  for( l= j+2 ; l<= k ; l++ )
	    minVal= minVal<= lineBuf[l] ? minVal : lineBuf[l];
	}
        // else, the new min is the same as the old one
  }
  
  // clean up
  delete [] lineBuf;
}



/** apply dilation filter to one channel of an array
 **
 ** This code is mostly copied from Linear1DFilter<T>::apply, but with
 ** a min operator in the center loop.
 */
template<class T>
bool
Dilate1D<T>::apply( Array<T> &array,
		    BoundaryMethod boundary,
		    unsigned int axis,
		    unsigned int inChannel,
		    unsigned int outChannel )
{
  // renormalize is identical to clamp for this filter
  if( boundary== Renormalize )
    boundary= Clamp;
  
  // call superclass method for actually applying the filter
  return Filter1D<T>::apply( array, boundary, axis, inChannel, outChannel );
}

/** apply filter to a single line in the array */
template <class T>
void
Dilate1D<T>::apply( T *startPos, unsigned long incr,
		    unsigned long numElements, BoundaryMethod boundary,
		    T background, T *startPosOut )
{
  unsigned long j, k, l;
  
  // allocate temporary memory for one line
  T* lineBuf= new T[numElements+2*radius];
  
  // fetch line into temp buffer
  fetchLine<T>( lineBuf, startPos, incr, numElements, radius,
		boundary, background );
  
  // apply filter:
  // find maximum over first 2*radius+1 elements
  T maxVal= lineBuf[0];
  for( j= 1 ; j< 2*radius+1 ; j++ )
    maxVal= maxVal>= lineBuf[j] ? maxVal : lineBuf[j];
  
  // simple application to the full line
  for( j= 0, k= 2*radius+1 ; j< numElements ; j++, k++ )
  {
    startPosOut[j*incr]= maxVal;
    // update maxVal
    if( j< numElements-1 )
      if( lineBuf[k]> maxVal )
	// if the new value added to the sliding window is greater
	// than current max, then that is the new max
	maxVal= lineBuf[k];
      else
	if( maxVal== lineBuf[j] )
	{
	  // if the maximum fell out of the sliding window, we
	  // have to find it from scratch
	  maxVal= lineBuf[j+1];
	  for( l= j+2 ; l<= k ; l++ )
	    maxVal= maxVal>= lineBuf[l] ? maxVal : lineBuf[l];
	}
        // else, the new max is the same as the old one
  }
  
  // clean up
  delete [] lineBuf;
}


/* we use explicit instantiation for classes in this file */
template class Erode1D<float>;
template class Erode1D<double>;
template class Dilate1D<float>;
template class Dilate1D<double>;


} /* namespace */

#endif /* FILTERS_MORPHOLOGICALOPS_C */

