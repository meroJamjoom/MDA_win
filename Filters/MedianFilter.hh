// ==========================================================================
// $Id: MedianFilter.hh 256 2008-12-09 02:41:53Z bradleyd $
// a median filter
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

#ifndef FILTERS_MEDIANFILTER_H
#define FILTERS_MEDIANFILTER_H

/*! \file  MedianFilter.hh
    \brief a median filter
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>

#include "MDA/Threading/SMPJob.hh"
#include "Filter.hh"
#include "MedianTable.hh"

namespace MDA {

  using namespace std;
  
  // forward declaration
  template <class T> class LineMedianJob;

  /** \class MedianFilter MedianFilter.hh
      a median filter */
  template<class T>
  class MedianFilter: public Filter<T> {

  public:

    /** default constructor */
    MedianFilter( unsigned rad= 1, unsigned bins= 256, bool quant= true )
      : radius( rad ), numBins( bins ), quantized( quant )
    {}
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
    
  protected:
    
    friend class LineMedianJob<T>;
    
    /** apply median filter to a single line of the array */
    virtual void apply(  T *data, short *hashes, CoordinateVector &pos,
		 AxisList &axes,
		 CoordinateVector &paddedIncr, unsigned long maxNumIter, 
		 unsigned long incr, unsigned long numElements,
		 BoundaryMethod boundary, T background,
		 T *startPosOut, unsigned long incrOut );

    /** whether we compute a quantized or a continuous mean */
    bool quantized;
    
    /** number of bins */
    unsigned numBins;
    
    /** radius of the median */
    unsigned radius;
  };




  /** \class MedianFilterMasked MedianFilter.hh
      a masked median filter */
  template<class T>
  class MedianFilterMasked: public MedianFilter<T> {

  public:
    
    /** default constructor */
    MedianFilterMasked( unsigned rad= 1, unsigned bins= 256, bool quant= true )
      : MedianFilter<T>(rad,bins,quant), paddedMaskData(NULL)
    {}
    
    ~MedianFilterMasked()
    {
      if (paddedMaskData)
	delete [] paddedMaskData;
      paddedMaskData = NULL;
    }

    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
  protected:
    
    friend class LineMedianJob<T>;
    
    /** apply median filter to a single line of the array */
    virtual void apply(  T *data, short *hashes, CoordinateVector &pos,
		 AxisList &axes,
		 CoordinateVector &paddedIncr, unsigned long maxNumIter, 
		 unsigned long incr, unsigned long numElements,
		 BoundaryMethod boundary, T background,
		 T *startPosOut, unsigned long incrOut );

  private:
    T* paddedMaskData;
  };


  
  /** \class LineMedianJob MedianFilter.hh
      Application of a MedianFilter to a single line in the array  */
  template<class T>
  class LineMedianJob: public SMPJob {
    
  public:
    
    /** constructor */
    LineMedianJob( MedianFilter<T> *f,
		   T* d, short *h, CoordinateVector &p, AxisList &a,
		   CoordinateVector &pI, unsigned long nI,
		   unsigned long inc, unsigned long nE, BoundaryMethod b,
		   T ba, T* so, unsigned long io )
      : filter( f ), data( d ), hashes( h ), pos( p ), axes( a ),
	paddedIncr( pI ), maxNumIter( nI ), incr( inc ), numElements( nE ),
	boundary( b ), background( ba ), startPosOut( so ), incrOut( io )
    {}
    
    /** execute job (pure virtual)
     * \param threadID is an int that identifies individual threads
     * primarily for debugging
     */
    virtual void execute( int threadID );

  protected:
    
    /** pointer to Filter1D with all the details */
    MedianFilter<T> *filter;
    
    /** pointer to padded channel data */
    T *data;
    
    /** an array of precomputed hashes for each element */
    short *hashes;
    
    /** current position (beginning of line) */
    CoordinateVector pos;
    
    /** list of axes along which to filter */
    AxisList axes;
    
    /** array increments along all major coordinate axes */
    CoordinateVector paddedIncr;
    
    /** the number of iterators required to describe the neighborhood
	of a pixel */
    unsigned long maxNumIter;
    
    /** increment along axis[0] in padded array */
    unsigned long incr;
    
    /** number of elements in line */
    unsigned long numElements;
    
    /** boundary mode */
    BoundaryMethod boundary;
    
    /** background value */
    T background;
    
    /** pointer to first element of line (output) */
    T* startPosOut;
    
    /** element increment in output array */
    unsigned long incrOut;
  };



} /* namespace */



#endif /* FILTERS_MEDIANFILTER_H */

