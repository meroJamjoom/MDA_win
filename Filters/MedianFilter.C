// ==========================================================================
// $Id: MedianFilter.C 259 2008-12-30 09:15:14Z heidrich $
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

#ifndef FILTERS_MEDIANFILTER_C
#define FILTERS_MEDIANFILTER_C


#include "MDA/Threading/SMPJobManager.hh"

#include "MedianFilter.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;



//
// MedianFilter members
//


/** apply the filter to a number of dimensions and channels */
template<class T>
bool
MedianFilter<T>::apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )
{
  long i, j, k, l;

  // dimensions, and dimensional pointer increments
  CoordinateVector dim=  a.getDimension();
  const int dimension= dim.vec.size();
  const int numAxes= axes.vec.size();
  CoordinateVector paddedDim;
  for( i= 0, j= 1 ; i< dimension ; i++ )
  {
    paddedDim.vec.push_back( dim.vec[i] + 2*radius );
    j*= paddedDim.vec[i];
  }
  const unsigned long paddedSize= j;
  T *paddedData= new T[paddedSize]; // temp buffer for boundary-padded channel
  short *hashes= new short[paddedSize];
  CoordinateVector incr;
  CoordinateVector paddedIncr;
  incr.vec.push_back( 1 );
  paddedIncr.vec.push_back( 1 );
  for( i= 1 ; i< dimension ; i++ )
  {
    incr.vec.push_back( incr.vec[i-1] * dim.vec[i-1] );
    paddedIncr.vec.push_back( paddedIncr.vec[i-1] * paddedDim.vec[i-1] );
  }

  // compute maximum number of iterators required, and allocate them
  unsigned long maxNumIter= 1;
  for( i= 1 ; i< numAxes ; i++ )
    maxNumIter*= 2*radius+1;
  
  // number of scanlines to process
  int numScanlines= 1;
  for( i= 0 ; i< dimension ; i++ )
    if( i!= axes.vec[0] )
      numScanlines*= dim.vec[i];
  
  // process all channels
  for( k= 0 ; k< channels.vec.size() ; k++ )
  {
    // get channel data, then copy and pad it according to boundary
    // mode. that allows us to overwrite the source channel with the
    // result directly
    T *channelData= &((*a[channels.vec[k]])[0]);
    T background= a[channels.vec[k]]->getBackground();
    fetchChannel( paddedData, channelData, dim, radius, boundary, background );
    
    // precompute hash values
    MedianTable<T> *hash;
    if( quantized )
      hash= new QuantizedMedianTable<T>( 0.0, 1.0, numBins );
    else
      hash= new ContinuousMedianTable<T>( 0.0, 1.0, numBins );
    for( i= 0 ; i< paddedSize ; i++ )
      hashes[i]= hash->getHash( paddedData[i] );
    delete hash;
    
    // initial position relative to original resolution
    CoordinateVector pos;
    for( i= 0 ; i< dimension ; i++ )
      pos.vec.push_back( 0 );
    
    // create a job list all scanlines in the channel
    SMPJobList jobs;
    for( l= 0 ; l< numScanlines ; l++ )
    {
      // set index to first pixel of destination line
      unsigned long index= 0;
      for( i= 0 ; i< dimension ; i++ )
	index+= incr.vec[i] * pos.vec[i];

      // process one scanline
      jobs.push_back( new LineMedianJob<T>( (MedianFilter<T> *)this,
					    paddedData, hashes, pos, axes,
					    paddedIncr, maxNumIter,
					    paddedIncr.vec[axes.vec[0]],
					    dim.vec[axes.vec[0]],
					    boundary, background, 
					    channelData+index,
					    incr.vec[axes.vec[0]] ) );
      
      // update position to next scanline
      for( i= 0 ; i< dimension ; i++ )
	if( i!= axes.vec[0] )
	  if( ++pos.vec[i]>= dim.vec[i] )
	    pos.vec[i]= 0;
	  else
	    break;
    }
    // batch-process the job list
    SMPJobManager::getJobManager()->batch( jobs );
  }
  
  delete [] hashes;
  delete [] paddedData;
  return true;
}

/** apply median filter to a single line of the array */
template <class T>
void
MedianFilter<T>::apply( T *data, short *hashes, CoordinateVector &pos,
			AxisList &axes, CoordinateVector &paddedIncr,
			unsigned long maxNumIter, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut, unsigned long incrOut )
{
  unsigned long i, j;
  
  // set addOff and delOff arrays as all pixels in neighborhood of
  // start pixel, orthogonal to axes[0]
  unsigned dimension= pos.vec.size();
  unsigned long *addOff= new unsigned long[maxNumIter];
  unsigned long *delOff= new unsigned long[maxNumIter];
  unsigned long numIter;
  CoordinateVector tmpPos= pos;
  for( i= 0 ; i< dimension ; i++ )
    tmpPos.vec[i]+= radius;
  unsigned numAxes= axes.vec.size();
  for( i= 0 ; i< numAxes ; i++ )
    tmpPos.vec[axes.vec[i]]-= radius;
  for( i= numIter= 0 ; i< maxNumIter ; i++ )
  {
    // compute index into (padded) source array
    unsigned long srcIndex= 0;
    for( j= 0 ; j< dimension ; j++ )
      srcIndex+= paddedIncr.vec[j] * tmpPos.vec[j];
    // set pointers
    delOff[numIter]= addOff[numIter]= srcIndex;
    numIter++;
    // update tmpPos
    for( j= 1 ; j< numAxes ; j++ )
    {
      if( ++tmpPos.vec[axes.vec[j]]> pos.vec[axes.vec[j]]+2*radius )
	tmpPos.vec[axes.vec[j]]= pos.vec[axes.vec[j]];
      else
	break;
    }
  }
  
  // create hash table & prime it a with the front 2*radius+1 elements
  MedianTable<T> *hash;
  if( quantized )
    hash= new QuantizedMedianTable<T>( 0, 1, numBins );
  else
    hash= new ContinuousMedianTable<T>( 0, 1, numBins );
  hash->clear();
  for( j= 2*radius+1 ; j> 0 ; j-- )
    for( i= 0 ; i< numIter ; i++ )
    {
      hash->add( data[addOff[i]], hashes[addOff[i]] );
      addOff[i]+= incr;
    }
  
  // foreach pixel in the scanline/column, get the median, then
  // update the median table by advancing the pointer array
  for( j= numElements ; j> 0 ; j--, startPosOut+= incrOut )
  {
    *startPosOut= hash->getMedian();
    for( i= 0 ; i< maxNumIter ; i++ )
    {
      hash->replace( data[delOff[i]], hashes[delOff[i]],
		     data[addOff[i]], hashes[addOff[i]] );
      addOff[i]+= incr;
      delOff[i]+= incr;
    }
  }
  
  // clean up
  delete hash;
  delete [] delOff;
  delete [] addOff;
}


/******************************************************************************/

/** apply the filter to a number of dimensions and channels */
template<class T>
bool
MedianFilterMasked<T>::apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )
{
  long i, j, k, l;

  // get mask data
  CoordinateVector dim=  a.getDimension();
  const int dimension= dim.vec.size();
  CoordinateVector paddedDim;
  for( i= 0, j= 1 ; i< dimension ; i++ )
  {
    paddedDim.vec.push_back( dim.vec[i] + 2*MedianFilter<T>::radius );
    j*= paddedDim.vec[i];
  }
  const unsigned long paddedSize= j;
  paddedMaskData= new T[paddedSize]; 
  int maskChannel=channels.vec.size()-1;
  T* maskChannelData= &((*a[channels.vec[maskChannel]])[0]);
  T background= a[channels.vec[maskChannel]]->getBackground();
  fetchChannel( paddedMaskData, maskChannelData, dim, MedianFilter<T>::radius, boundary, background );

  // strip off the last channel (the mask)
  channels.vec.pop_back();

  return MedianFilter<T>::apply(a, boundary, channels, axes);

}

/** apply median filter to a single line of the array */
template <class T>
void
MedianFilterMasked<T>::apply( T *data, short *hashes, CoordinateVector &pos,
			AxisList &axes, CoordinateVector &paddedIncr,
			unsigned long maxNumIter, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut, unsigned long incrOut )
{
  unsigned long i, j;
  
  // set addOff and delOff arrays as all pixels in neighborhood of
  // start pixel, orthogonal to axes[0]
  unsigned dimension= pos.vec.size();
  unsigned long *addOff= new unsigned long[maxNumIter];
  unsigned long *delOff= new unsigned long[maxNumIter];
  unsigned long numIter;
  CoordinateVector tmpPos= pos;
  for( i= 0 ; i< dimension ; i++ )
    tmpPos.vec[i]+= MedianFilter<T>::radius;
  unsigned numAxes= axes.vec.size();
  for( i= 0 ; i< numAxes ; i++ )
    tmpPos.vec[axes.vec[i]]-= MedianFilter<T>::radius;
  

  for( i= numIter= 0 ; i< maxNumIter ; i++ )
  {
    // compute index into (padded) source array
    unsigned long srcIndex= 0;
    for( j= 0 ; j< dimension ; j++ )
      srcIndex+= paddedIncr.vec[j] * tmpPos.vec[j];

    // set pointers
    delOff[numIter]= addOff[numIter]= srcIndex;
    numIter++;
    // update tmpPos
    for( j= 1 ; j< numAxes ; j++ )
    {
      if( ++tmpPos.vec[axes.vec[j]]> pos.vec[axes.vec[j]]+2*MedianFilter<T>::radius )
	tmpPos.vec[axes.vec[j]]= 0;
      else
	break;
    }
  }
  

  // create hash table & prime it a with the front 2*radius+1 elements
  MedianTable<T> *hash;
  if( MedianFilter<T>::quantized )
    hash= new QuantizedMedianTable<T>( 0, 1, MedianFilter<T>::numBins );
  else
    hash= new ContinuousMedianTable<T>( 0, 1, MedianFilter<T>::numBins );
  hash->clear();
  for( j= 2*MedianFilter<T>::radius+1 ; j> 0 ; j-- )
    {
      for( i= 0 ; i< numIter ; i++ )
	{
	  if (paddedMaskData[addOff[i]]>0.0)
	    hash->add( data[addOff[i]], hashes[addOff[i]] );
	  addOff[i]+= incr;
	}
    }

  unsigned long maskIndex = 0;
  for( j= 0 ; j< dimension ; j++ )
    maskIndex+= paddedIncr.vec[j] * (pos.vec[j] + MedianFilter<T>::radius);

  // foreach pixel in the scanline/column, get the median, then
  // update the median table by advancing the pointer array
  for( j= numElements ; j> 0 ; j--, startPosOut+= incrOut )
  {
    if (paddedMaskData[maskIndex]>0.0)
      *startPosOut= hash->getMedian();
    maskIndex += paddedIncr.vec[0];

    for( i= 0 ; i< maxNumIter ; i++ )
    {
      if (paddedMaskData[delOff[i]]>0.0)
	hash->remove(data[delOff[i]], hashes[delOff[i]]);
      if (paddedMaskData[addOff[i]]>0.0)
	hash->add(data[addOff[i]], hashes[addOff[i]]);

      addOff[i]+= incr;
      delOff[i]+= incr;
    }
  }
  
  // clean up
  delete hash;
  delete [] delOff;
  delete [] addOff;
}


/** execute job (pure virtual)
 * \param threadID is an int that identifies individual threads
 * primarily for debugging
 */
template <class T>
void
LineMedianJob<T>::execute( int thredID )
{
  filter->apply( data, hashes, pos, axes, paddedIncr, maxNumIter, incr,
		 numElements, boundary, background, startPosOut, incrOut );
}


// template instantiation code
template class MedianFilter<float>;
template class MedianFilter<double>;
template class MedianFilterMasked<float>;
template class MedianFilterMasked<double>;
template class LineMedianJob<float>;
template class LineMedianJob<double>;

} /* namespace */

#endif /* FILTERS_MEDIANFILTER_C */

