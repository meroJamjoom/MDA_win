// ==========================================================================
// $Id: Linear1DFilter.hh 319 2009-05-27 21:17:01Z heidrich $
// A regularly sampled linear 1D filter class
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

#ifndef FILTERS_LINEAR1DFILTER_H
#define FILTERS_LINEAR1DFILTER_H

/*! \file  Linear1DFilter.hh
    \brief A regularly sampled linear 1D filter class
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <string.h>
#include <vector>

#include "MDA/Array/Array.hh"
#include "Filter1D.hh"


namespace MDA {

  /** \class Linear1DFilter Linear1DFilter.hh
      A regularly sampled linear 1D filter class */
  
  template<class T>
  class Linear1DFilter: public Filter1D<T> {
    
  public:
    
    /** constructor */
    Linear1DFilter( unsigned int radius )
    {
      Linear1DFilter<T>::radius= radius;
      filter= new double[2*radius+1];
    }
    
    /** destructor */
    virtual ~Linear1DFilter()
    {
      delete [] filter;
    }
    
    /** apply filter to one channel of an array */
    virtual bool apply( Array<T> &array, BoundaryMethod boundary,
			unsigned int axis, unsigned int inChannel,
			unsigned int outChannel );
    
    /** set filter contents. */
    inline void setFilter( unsigned radius, double *filter )
    {
      Linear1DFilter<T>::radius= radius;
      memcpy( Linear1DFilter<T>::filter, filter, (2*radius+1)*sizeof(double) );
    }
    
    /** get filter contents */
    inline const double *getFilter() const
    {
      return filter;
    }
    
  protected:
    
    /** provides an estimate of the computational affort involved in
	processing a certain number of elements (i.e #elem*(radius*2+1)) */
    virtual double getLineCost( unsigned long numElements );

    /** apply filter to a single line in the array */
    virtual void apply( T *startPos, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut );
    
    /** filter radius */
    unsigned int radius;
    
    /** the filter array */
    double *filter;
  };

  
  /** \class BoxFilter1D Linear1DFilter.hh
      A 1D box filter */
  template<class T>
  class BoxFilter1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    BoxFilter1D( unsigned radius= 1, double multiplier= 1.0 );
  };
  
  
  /** \class HatFilter1D Linear1DFilter.hh
      A 1D piecewise linear (i.e. hat) filter */
  template<class T>
  class HatFilter1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    HatFilter1D( unsigned radius= 1, double multiplier= 1.0 );
  };
  
  
  /** \class Gaussian1D Linear1DFilter.hh
      A 1D Gaussian filter */
  template<class T>
  class Gaussian1D: public Linear1DFilter<T> {
  public:
    /** constructor
	\param sigma: standard deviation
	\param _radius: radius of filter (if -1: use 2*sigma)
	\param multiplier: integral over the Gaussian */
    Gaussian1D( double sigma, int radius= -1,
		double multiplier= 1.0 );
  };

  
  /** \class FastGaussian1D Linear1DFilter.hh
      Algorithm from Young/van Viet,
      "Recursive Implementation of the Gaussian Filter",
      Signal Processing 44 (1995), 139-151 */
  template<class T>
  class FastGaussian1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    FastGaussian1D( double sigma, double multiplier= 1.0 );
 
    /** apply filter to array channel */
    virtual bool apply( Array<T> &array, BoundaryMethod boundary,
			unsigned int axis, unsigned int inChannel,
			unsigned int outChannel );
    
  protected:

    /** apply filter to a single line in the array */
    virtual void apply( T *startPos, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut );
    
    /** standard deviation */
    double sigma;
  };
  
  
  /** \class LaplacianOfGaussian1D Linear1DFilter.hh
      A 1D Laplacian-of-Gaussian filter */
  template<class T>
  class LaplacianOfGaussian1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    LaplacianOfGaussian1D( double sigma, int radius= -1,
			   double multiplier= 1.0 );
  };
  
  
  /** \class FirstDerivative1D Linear1DFilter.hh
      A 1D (divided difference) derivative filter */
  template<class T>
  class FirstDerivative1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    FirstDerivative1D( unsigned baseline= 1 );
  };
  
  
  /** \class SecondDerivative1D Linear1DFilter.hh
      A 1D (divided difference) filter of the second derivative */
  template<class T>
  class SecondDerivative1D: public Linear1DFilter<T> {
  public:
    /** constructor */
    SecondDerivative1D( unsigned baseline= 1 );
  };
  
  
} /* namespace */

#endif /* ARRAY_LINEAR1DFILTER_H */

