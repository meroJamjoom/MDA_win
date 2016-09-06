// ==========================================================================
// $Id: PointSampler.hh 386 2009-10-03 04:07:57Z heidrich $
// interface for sampling an array at individual points
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

#ifndef RESAMPLING_POINTSAMPLER_H
#define RESAMPLING_POINTSAMPLER_H

/*! \file  PointSampler.hh
    \brief interface for sampling an array at individual points
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CoordinateVector.hh"
#include <MDA/Array/Boundary.hh>
#include <MDA/LinearAlgebra/Vector.hh>

namespace MDA {

  /** \class PointSampler PointSampler.hh
      interface for sampling an array at individual points */
  template<class T>
  class PointSampler {

  public:

    /** constructor */
    PointSampler( CoordinateVector &arrayDimension,
		  unsigned long supportWidth );
    
    /** detructor */
    ~PointSampler();
    
    /** set the location of the point sampling process
     *  (if radius is 0, the standard magnification-case radius of the
     *  individual filter is used)
     */
    virtual void setSampleLocation( Vector pos,
				    BoundaryMethod boundary );
    
    /** return sample value as extracted from specific array data */
    inline T getSample( T *data ) const
    {
      T result= 0.0;
      for( unsigned i= 0 ; i< supSize ; i++ )
	result+= supWeights[i] * data[supOffsets[i]];
      return result;
    }
    
    /** splat a value into a specific array at the current sample position
     *  (the splat can either be additive, or replace all pixel values in
     *  the support)
     */
    inline void splatSample( T value, T *data, bool add= true ) const
    {
      for( unsigned i= 0 ; i< supSize ; i++ )
      {
	long off= supOffsets[i];
	if( off>= 0 )
	  data[off]= add*data[off] + supWeights[i] * value;
      }
    }
    
    /** pure virtual methods for doing sequences of 1D interpolations
     *	for determining the filter coefficients.
     *  \param numSamples: number of current values in pixel array (should be
     *  power of support radius)
     *  \param x: fractional part of coordinate along the current axis
     *  \param radius: filter radius (0 means use default)
     *  \returns new number of values in pixel array
     */
    virtual unsigned interp( unsigned numSamples, double x )=0;
    
  protected:
    
    /** calculate the array offset for a given integer position
     *  (including boundary effects)
     *  \param pos: position vector
     *  \param boundary: boundary method
     */
    unsigned long calculateOffset( CoordinateVector &pos,
				   BoundaryMethod boundary );
    
    /** dimension of array */
    CoordinateVector dim;
    
    /** scalar dimension */
    unsigned dimension;
    
    /** support width along one dimension */
    unsigned supWidth;
    
    /** number of pixels in a support (= supWidth ^ dimension) */
    unsigned supSize;
    
    /** array indices for each pixel in support (signed, because they
	can be -1 to represent background color) */
    long *supOffsets;
    
    /** weights for each pixel in support */
    T *supWeights;

    /** transient position variable (avoids repeated memory reallocation) */
    CoordinateVector currentPos;
  };


  
  /** \class NearestNeighborPointSampler PointSampler.hh
      point sampling with nearest neighbor "interpolation" */
  template<class T>
  class NearestNeighborPointSampler: public PointSampler<T> {

  public:

    /** constructor */
    inline NearestNeighborPointSampler( CoordinateVector &arrayDimension )
      : PointSampler<T>( arrayDimension, 1 )
    {
      PointSampler<T>::supWeights[0]= 1.0;
    }
    
    /** set the location of the point sampling process - special
	implementation for nearest neighbor */
    virtual void setSampleLocation( Vector pos,
				    BoundaryMethod boundary );
    
    /** this method does not get called for Nearest neighbor... */
    virtual unsigned interp( unsigned numSamples, double x )
    {
      return numSamples;
    }
    
  };
  
  
  
  /** \class LinearPointSampler PointSampler.hh
      point sampling with linear interpolation */
  template<class T>
  class LinearPointSampler: public PointSampler<T> {

  public:

    /** constructor */
    inline LinearPointSampler( CoordinateVector &arrayDimension )
      : PointSampler<T>( arrayDimension, 2 )
    {}
    
    /** create interpolation weights by interpolating along one axis... */
    virtual unsigned interp( unsigned numSamples, double x );
    
  };

  

  /** \class CubicPointSampler PointSampler.hh
      point sampling with linear interpolation */
  template<class T>
  class CubicPointSampler: public PointSampler<T> {

  public:

    /** constructor */
    inline CubicPointSampler( CoordinateVector &arrayDimension )
      : PointSampler<T>( arrayDimension, 4 )
    {}
    
    /** create interpolation weights by interpolating along one axis... */
    virtual unsigned interp( unsigned numSamples, double x );
    
  };

  
  
  /** \class GaussPointSampler PointSampler.hh
      point sampling with linear interpolation */
  template<class T>
  class GaussPointSampler: public PointSampler<T> {
    
  public:
    
    /** constructor */
    inline GaussPointSampler( CoordinateVector &arrayDimension )
      : PointSampler<T>( arrayDimension, 4 )
    {}
    
    /** create inteprolation weights by interpolating along one axis... */
    virtual unsigned interp( unsigned numSamples, double x );
    
  private:
    
    /** the exponential scale factor for a sigma that best
	approximates a sinc with a Gaussian */
    static const double scaleFactor;
    
  };

  

} /* namespace */

#endif /* RESAMPLING_POINTSAMPLER_H */

