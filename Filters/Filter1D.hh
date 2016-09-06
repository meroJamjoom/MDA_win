// ==========================================================================
// $Id: Filter1D.hh 247 2008-11-25 05:23:41Z heidrich $
// A 1D filter abstract baseclass
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

#ifndef FILTERS_FILTER1D_H
#define FILTERS_FILTER1D_H

/*! \file  Filter1D.hh
    \brief A 1D filter abstract baseclass
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Threading/SMPJob.hh"
#include "MDA/Array/Array.hh"
#include "Filter.hh"


namespace MDA {
  
  // forward declarartion
  template <class T> class Filter1DLineJob;
  
  /** \class Filter1D Filter1D.hh
      A 1D filter abstract baseclass */
  template<class T>
  class Filter1D {

  public:

    /** apply filter to one channel of an array (pure virtual) */
    virtual bool apply( Array<T> &array, BoundaryMethod boundary,
			unsigned int axis, unsigned int inChannel,
			unsigned int outChannel );
    
    /** destructor */
    virtual ~Filter1D() {};
    
  protected:
    
    friend class Filter1DLineJob<T>;
    
    /** provides an estimate of the computational affort involved in
	processing a certain number of elements (override to optimize
	threading) */
    virtual double getLineCost( unsigned long numElements );

    /** apply filter to a single line in the array */
    virtual void apply( T *startPos, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut )= 0;
    
  };
  
  
  /** \class Filter1DLineJob Filter1D.hh
      Application of a Filter1D to a single line in the array  */
  template<class T>
  class Filter1DLineJob: public SMPJob {
    
  public:
    
    /** constructor */
    Filter1DLineJob( Filter1D<T> *f, T* s, unsigned long i,
		     unsigned long n, BoundaryMethod b, T ba, T* so )
      : SMPJob( n*10 ), filter( f ), startPos( s ), incr( i ), numElements( n ), boundary( b ),
	background( ba ), startPosOut( so )
    {}
    
    /** execute job (pure virtual)
     * \param threadID is an int that identifies individual threads
     * primarily for debugging
     */
    virtual void execute( int threadID );

  protected:
    
    /** pointer to Filter1D with all the details */
    Filter1D<T> *filter;
    
    /** pointer to first element of line (input) */
    T* startPos;
    
    /** pointer increment between each element along the line */
    unsigned long incr;
    
    /** number of elements along the line */
    unsigned long numElements;
    
    /** which method to use for boundary padding */
    BoundaryMethod boundary;
    
    /** pointer to first element of line (output) */
    T* startPosOut;
    
    /** background value */
    T background;
    
  };
  
} /* namespace */


#endif /* FILTERS_FILTER1D_H */

