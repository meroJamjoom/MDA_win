// ==========================================================================
// $Id: Scaling.C 712 2010-04-23 09:16:34Z martinle $
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

#ifndef GEOMETRICTRANSFORM_SCALING_C
#define GEOMETRICTRANSFORM_SCALING_C

#include <string.h>

#include "MDA/Base/Errors.hh"
#include "MDA/Threading/SMPJobManager.hh"
#include "MDA/Resampling/NormalizedDeviceCoordinates.hh"

#include "Scaling.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** apply the transform to a single channel in one array */
template<class T>
bool
Scaling<T>::apply( Array<T> &srcArray, int srcChannel,
		   Array<T> &dstArray, int dstChannel )
{
  unsigned long i, j;
  
  // consistency checks
  CoordinateVector srcDim= srcArray.getDimension();
  CoordinateVector dstDim= dstArray.getDimension();
  int dimension= srcDim.vec.size();
  if( !warnCond( dimension== dstDim.vec.size(),
		 "  input and output dimensions do not match" ) )
    return false;
  if( !warnCond( dimension== scaling.getNumRows()-1,
		 "  array dimensions don't match scaling coefficients" ) )
    return false;
  
  // mapping from destination pixel space to source pixel space
  // (first map dst->NDC, then scale, then map NDC->dst)
  NormalizedDeviceCoordinates srcNDC( srcDim );
  NormalizedDeviceCoordinates dstNDC( dstDim );
  Matrix tmpMat( dimension+1, dimension+1 );
  Matrix pixelXform( dimension+1, dimension+1 );
  multMatrixMatrix( scaling, dstNDC.getPixelToNDCMatrix(), tmpMat );
  multMatrixMatrix( srcNDC.getNDCToPixelMatrix(), tmpMat, pixelXform );
  
  // set up buffer pointers
  CoordinateVector currDim, nextDim;
  T *srcBuf= &((*srcArray[srcChannel])[0]);
  T *dstBuf= &((*dstArray[dstChannel])[0]);
  T *currBuf, *nextBuf;
  unsigned long nextSize;
  unsigned long stride= 1;
  
  // scale along each axis
  for( i= 0, currDim= srcDim, currBuf= srcBuf ;
       i< dimension ;
       i++, currDim= nextDim, currBuf= nextBuf )
  {
    nextDim= currDim;
    nextDim.vec[i]= dstDim.vec[i];
    
    // determine destination buffer
    if( i== dimension-1 )
      // last axis: write to destination array
      nextBuf= dstBuf;
    else
    {      
      // compute size of next array, and allocate buffer
      for( j= dimension, nextSize= 1 ; j> 0 ; )
	nextSize*= nextDim.vec[--j];
      nextBuf= new T[nextSize];
    }
      
    // number of scanlines for this axis
    unsigned long numScanlines= 1;
    for( j= 0 ; j< dimension ; j++ )
      if( j!= i )
	numScanlines*= currDim.vec[j];
    
    // scale a and offset b for the resamplers along axis i
    double a= pixelXform[i][i];
    double b= pixelXform[i][dimension];
    
    // scale all scanlines: put all scanlines into a job list, and
    // batch process it
    SMPJobList jobs;
    Resampler<T> *resampler= GeometricTransformation<T>::resampler;
    unsigned long currOff, nextOff;
    for( j= currOff= nextOff= 0 ; j< numScanlines ; j++ )
    {
      jobs.push_back( new LinearResamplerJob<T>( resampler,
						 currBuf+currOff,
						 currDim.vec[i],
						 stride,
						 nextBuf+nextOff,
						 nextDim.vec[i],
						 stride,
						 a, b ) ); 
      
      // update offsets
      currOff++;
      nextOff++;
      if( currOff % stride== 0 )
      {
	currOff+= stride * (currDim.vec[i]-1);
	nextOff+= stride * (nextDim.vec[i]-1);
      }
    }
    SMPJobManager::getJobManager()->batch( jobs );
    
    // update stride
    stride*= nextDim.vec[i];
    
    // clean up temp buffers that are no longer required
    if( currBuf!= srcBuf )
      delete [] currBuf;
  }
  
  return true;
}



//
// template instantiations
//

template class Scaling<float>;
template class Scaling<double>;
template class Scaling<unsigned char>;

} /* namespace */

#endif /* GEOMETRICTRANSFORM_SCALING_C */

