// ==========================================================================
// $Id: Filter1D.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef FILTERS_FILTER1D_C
#define FILTERS_FILTER1D_C

#include "MDA/Base/Errors.hh"
#include "MDA/Threading/SMPJobManager.hh"

#include "Filter1D.hh"

namespace MDA {

/** apply filter to one channel of an array (pure virtual) */
template <class T>
bool
Filter1D<T>::apply( Array<T> &array, BoundaryMethod boundary,
		    unsigned int axis, unsigned int inChannel,
		    unsigned int outChannel )
{
  unsigned long i;
  unsigned long startPos;
  
  //
  // checks
  //
  
  // get & check channel information
  typename Array<T>::Channel *in=  array[inChannel];
  typename Array<T>::Channel *out= array[outChannel];
  if( !warnCond(  in!= NULL, "  input channel out of range" ) ||
      !warnCond( out!= NULL, "  output channel out of range" ) )
    return false;
  
  // get & check dimension information
  CoordinateVector dim=  array.getDimension();
  if( !warnCond( axis< dim.vec.size(), "  axis out of range" ) )
    return false;
  unsigned long numElements= dim.vec[axis];
  
  // intra-channel array increment from one entry to the next along the axis 
  unsigned long incr= 1;
  for( i= 0 ; i< axis ; i++ )
    incr*= dim.vec[i];
  // total number of lines
  unsigned long numLines= 1;
  for( i= 0 ; i< dim.vec.size() ; i++ )
    if( i!= axis )
      numLines*= dim.vec[i];
  
  //
  // apply the filter to all lines
  //

  // first put all jobs in a job queue
  SMPJobList jobs;
  for( i= startPos= 0 ; i< numLines ; i++ )
  {
    // push the line job into the job list
    SMPJob *job= new Filter1DLineJob<T>( (Filter1D<T> *)this,
					 &(*in)[startPos], incr,
					 numElements, boundary,
					 in->getBackground(),
					 &(*out)[startPos] );
    job->timeEstimate= getLineCost( numElements );
    jobs.push_back( job );
  
    // update starting position
    // the increment is 1 (i.e next element along axis 0, then axis 1 etc
    // EXCEPT if the increment is along "axis"
    if( ++startPos % incr == 0 )
      startPos+= (numElements-1)*incr;
  }
  
  // then execute the jobs in parallel
  SMPJobManager::getJobManager()->batch( jobs );
  
  return true;
}
  

/** provides an estimate of the computational affort involved in
    processing a certain number of elements (used to optimize
    threading) */
template <class T>
double
Filter1D<T>::getLineCost( unsigned long numElements )
{
  return numElements * 10;
}

  
/** execute job (pure virtual)
 * \param threadID is an int that identifies individual threads
 * primarily for debugging
 */
template <class T>
void
Filter1DLineJob<T>::execute( int threadID )
{
  filter->apply( startPos, incr, numElements, boundary,
		 background, startPosOut );
}

/* we use explicit instantiation for this class */
template class Filter1D<float>;
template class Filter1D<double>;
template class Filter1DLineJob<float>;
template class Filter1DLineJob<double>;

} /* namespace */

#endif /* FILTERS_FILTER1D_C */

