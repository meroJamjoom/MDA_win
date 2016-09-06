// ==========================================================================
// $Id: HitOrMiss2D.hh 372 2009-09-20 18:57:38Z heidrich $
// hit-or-miss transform
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

#ifndef FILTERS_HITORMISS_H
#define FILTERS_HITORMISS_H

/*! \file  HitOrMiss2D.hh
    \brief hit-or-miss transform
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <string.h>

#include "MDA/Threading/SMPJob.hh"
#include "Filter.hh"

namespace MDA {

  // forward declaration
  template <class T> class HitOrMissLineJob;


  /** \class HitOrMiss2D HitOrMiss2D.hh
      A 2D hit-and-miss transform (with radius 1)
      
      The class can be used as the usual 3x3 HOM transform, or as a
      general 3x3 pattern matching, with an arbitrary double-valued
      result depending on the neighborhood configuration.
  */
  template<class T>
  class HitOrMiss2D: public Filter<T> {

  public:

#if defined(_WIN32) || defined(_WIN64)
#define bzero(p, l) memset(p, 0, l)
#endif


    /** constructor from provided case table */
    HitOrMiss2D( double *cases= NULL, bool _preserveMisses= false,
		 T _hitValue= 1.0, T _missValue= 0.0 )
      : preserveMisses( _preserveMisses ),
	hitValue( _hitValue ), missValue( _missValue )
    {
      caseTable= new double[512];
      if( cases!= NULL )
	memcpy( caseTable, cases, 512*sizeof( double ) );
      else
	bzero( caseTable, 512*sizeof( double ) );
    }
    
    /** constructor from a traditional structural element (&mask)
     *  this constructor can automatically generate symmetry cases if
     *  desired. */
    HitOrMiss2D( bool structureElem[9], bool mask[9],
		 bool rotate= false, bool reflect= false, 
		 bool preserveMisses= false,
		 T hitValue= 1.0, T missValue= 0.0 );
    
    /** destructor */
    ~HitOrMiss2D()
    {
      delete [] caseTable;
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
    
    /** combine with another hit&miss transform by or-ing the case table */
    inline HitOrMiss2D<T> &operator|=( const HitOrMiss2D<T> &other )
    {
      for( unsigned i= 0 ; i< 512 ; i++ )
      {
	caseTable[i]= (caseTable[i]> 0.0 || other.caseTable[i]> 0.0) ? 1.0:0.0;
      }
      return *this;
    }
    
    /** combine with another hit&miss transform by and-ing the case table */
    inline HitOrMiss2D<T> &operator&=( const HitOrMiss2D<T> &other )
    {
      for( unsigned i= 0 ; i< 512 ; i++ )
	caseTable[i]= (caseTable[i]> 0.0 && other.caseTable[i]> 0.0) ? 1.0:0.0;
      return *this;
    }
    
    /** negate the case table */
    inline HitOrMiss2D<T> &operator!()
    {
      for( unsigned i= 0 ; i< 512 ; i++ )
	caseTable[i]= (caseTable[i] > 0.0) ? 0.0 : 1.0;
      return *this;
    }
    
  protected:
    
    /** convert a neighborhood bit vector to a table index */
    unsigned neighborhoodToIndex( const bool neighborhood[9] );
    
    /** convert a table index to a neighborhood bit vector */
    void indexToNeighborhood( unsigned index, bool neighborhood[9] );
    
    /** rotate a table index by 90 degrees ccw */
    unsigned rotate( unsigned ind );
    
    /** reflect a table index horizontally */
    unsigned reflect( unsigned ind );
    
    /** whether to preserve the original array content for a miss, or
        repace it with the miss value */
    bool preserveMisses;
    
    /** new pixel value for hits */
    T hitValue;
    
    /** new pixel value for misses (unless preserveMisses is true) */
    T missValue;
    
    /** case table */
    double *caseTable;
    
    /** apply to a single scanline */
    void apply( bool* currLine, T* dst,	unsigned long numElements );

    /** HitOrMissLineJob can execute apply method */
    friend class HitOrMissLineJob<T>;
  };



  /** \class HitOrMissLineJob HitOrMiss2D.hh
      A single line in an HMT */
  template<class T>
  class HitOrMissLineJob: public SMPJob {

  public:
    
    /** constructor */
    inline HitOrMissLineJob( HitOrMiss2D<T> *_filter, bool* _currLine,
			     T* _dst, unsigned long _numElements )
      : SMPJob( _numElements*3 ), filter( _filter ), currLine( _currLine ),
	dst( _dst ), numElements( _numElements )
    {}
    
    /** execute line job */
    virtual void execute( int jobID );

  protected:
    
    /** the filter */
    HitOrMiss2D<T> *filter;
    
    /** pointer to a buffer containign the current line as bools */
    bool *currLine;
    
    /** result pointer */
    T* dst;
    
    unsigned long numElements;
  };

} /* namespace */

#endif /* FILTERS_HITORMISS_H */

