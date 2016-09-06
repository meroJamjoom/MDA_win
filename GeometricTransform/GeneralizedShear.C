// ==========================================================================
// $Id:$
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

#ifndef GEOMETRICTRANSFORM_GENERALIZEDSHEAR_C
#define GEOMETRICTRANSFORM_GENERALIZEDSHEAR_C

#include "MDA/Base/Errors.hh"
#include "MDA/Threading/SMPJobManager.hh"
#include "MDA/Resampling/Resampling.hh"

#include "GeneralizedShear.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** apply the transform to a single channel in one array */
template<class T>
bool
GeneralizedShear<T>::apply( Array<T> &srcArray, int srcChannel,
			    Array<T> &dstArray, int dstChannel )
{
  unsigned long i, j;
  
  //
  // consistency checks
  //
  CoordinateVector srcDim= srcArray.getDimension();
  CoordinateVector dstDim= dstArray.getDimension();
  int dimension= srcDim.vec.size();
  if( !warnCond( dimension== dstDim.vec.size(),
                "  input/output dimensions do not match" ) )
    return false;
  
  if( !warnCond( axis>= dimension || coefficients.size() == (dimension+1),
                "  shear parameters do not match array" ) )
    return false;
  
  //
  // setup
  //
  unsigned long stride= 1;
  for( i= 0 ; i< axis ; i++ )
    stride*= srcDim.vec[i];
  unsigned long srcOffset= 0;
  unsigned long dstOffset= 0;
  T *srcPtr= &(*srcArray[srcChannel])[0];
  T *dstPtr= &(*dstArray[dstChannel])[0];
  
  CoordinateVector pos;
  for( i= 0 ; i< dimension ; i++ )
    pos.vec.push_back( 0 );
  
  // actual transformation:
  // we stuff all the scanlines into a job list, and batch-process it
  // at the end
  SMPJobList jobs;
  Resampler<T> *resampler= GeometricTransformation<T>::resampler;
  while( 1 )
  {
    double a= coefficients[axis];
    double b= coefficients[dimension];
    // compute shear equations
    for( i= 0 ; i< dimension ; i++ )
      b+= coefficients[i] * pos.vec[i];
    
    jobs.push_back( new LinearResamplerJob<T>( resampler,
					       srcPtr+srcOffset,
					       srcDim.vec[axis],
					       stride,
					       dstPtr+dstOffset,
					       dstDim.vec[axis],
					       stride,
					       a, b ) );
    
    // update position & offsets
    for( i= 0 ; i< dimension ; i++ )
      if( i!= axis )
	if( ++pos.vec[i]>= srcDim.vec[i] )
	  pos.vec[i]= 0;
	else
	  break;
    // terminate loop if we are through the whole array
    if( i== dimension )
      break;
    srcOffset++;
    dstOffset++;
    if( srcOffset % stride == 0 )
    {
      srcOffset+= stride * (srcDim.vec[axis]-1);
      dstOffset+= stride * (dstDim.vec[axis]-1);
    }
  }
  SMPJobManager::getJobManager()->batch( jobs );
  return true;
}


//
// template instantiations
//

template class GeneralizedShear<float>;
template class GeneralizedShear<double>;



} /* namespace */

#endif /* GEOMETRICTRANSFORM_GENERALIZEDSHEAR_C */

