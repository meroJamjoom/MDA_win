// ==========================================================================
// $Id: Resampler.hh 999 2014-05-28 15:07:31Z heidrich $
// Virtual base class for resampling filters
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

#ifndef RESAMPLING_RESAMPLER_H
#define RESAMPLING_RESAMPLER_H

/*! \file  Resampler.hh
 \brief Virtual base class for resampling filters
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#if defined (_WIN32) || defined (_WIN64)
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <math.h>
#endif


#include "MDA/Array/Boundary.hh"

namespace MDA {
  
  /** \class Resampler Resampler.hh
   Virtual base class for resampling filters */
  template<class T>
  class Resampler {
    
  public:
    
    /** default constructor */
    Resampler( BoundaryMethod boundaryMethod= Background, T backgroundValue= 0)
    {
      // enforce vitual table lookup
      this->setBoundaryMethod( boundaryMethod, backgroundValue );
    }
    
    /** destructor */
    virtual ~Resampler() {}
    
    /** set the background color */
    inline void setBackground( T backgroundValue )
    {
      background= backgroundValue;
    }
    
    /** set boundary handling */
    virtual void setBoundaryMethod( BoundaryMethod boundaryMethod,
				    T backgroundValue= 0 )
    {
      boundary= boundaryMethod;
      background= backgroundValue;
    }
    
    /** return boundary method */
    inline BoundaryMethod getBoundaryMethod()
    {
      return boundary;
    }
    
    /** get a pixel from the boundary */
    T getBoundaryPixel( T *inLine, long inCount,
			long inStride, long pos );
    
    /** resample a scanline/column with linear function oldX= a*newX + b
     
     (the default implementation is to use resampleIrregular, but
     subclasses can override with more efficient specialized
     implementations)
     */
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** resample a scanline/column with rational function
     oldX = (a*newX + b) / (c*newX + d)
     
     (the default implementation is to use resampleIrregular, but
     subclasses can override with more efficient specialized
     implementations)
     */
    virtual void resampleRational( T *inLine, unsigned long inCount,
				   unsigned long inStride, T *outLine,
				   unsigned long outCount,
				   unsigned long outStride,
				   double a, double b, double c, double d );
    
    /** resample a scanline/column at a specified set of sample points */
    virtual void resampleIrregular( T* inLine, unsigned long inCount,
				    unsigned long inStride, T* outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul )= 0;
    
  protected:
    
    /** mode for dealing with boundaries */
    BoundaryMethod boundary;
    
    /** value to be used as background in Background boundary mode */
    T background;
    
  };
  
  
  
  /** \class NearestNeighborResampler Resampler.hh
   class for resampling with nearest neighbor filter */
  template<class T>
  class NearestNeighborResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    NearestNeighborResampler( BoundaryMethod boundary= Background,
                             T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with NN reconstruction */
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** resample a scanline/column at a specified set of sample points */
    virtual void resampleIrregular( T* inLine, unsigned long inCount,
				    unsigned long inStride, T* outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
  };
  
  
  
  /** \class BoxResampler Resampler.hh
   class for resampling with box filter
   (similar to NearestNeighbor for magnification, but with better
   minification)*/
  template<class T>
  class BoxResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    BoxResampler( BoundaryMethod boundary= Background, T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with box reconstruction */
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** resample a scanline/column at a specified set of sample points */
    virtual void resampleIrregular( T* inLine, unsigned long inCount,
				    unsigned long inStride, T* outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
  };
  
  
  
  /** \class LinearResampler Resampler.hh
   class for resampling with linear (hat) filter */
  template<class T>
  class LinearResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    LinearResampler<T>( BoundaryMethod boundary= Background, T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with linear (hat) reconstruction*/
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** rational resampling (oldX = (a*newX + b) / (c*newX + d)),
     linear (hat) reconstruction*/
    virtual void resampleIrregular( T *inLine, unsigned long inCount,
				    unsigned long inStride, T *outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
    
  };
  
  
  
  /** \class CubicResampler Resampler.hh
   class for resampling with linear (hat) filter */
  template<class T>
  class CubicResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    CubicResampler<T>( BoundaryMethod boundary= Background, T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with linear (hat) reconstruction*/
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** rational resampling (oldX = (a*newX + b) / (c*newX + d)),
     linear (hat) reconstruction*/
    virtual void resampleIrregular( T *inLine, unsigned long inCount,
				    unsigned long inStride, T *outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
    
  };
  
  
  
  /** \class GaussResampler Resampler.hh
   class for resampling with Gauss filter */
  template<class T>
  class GaussResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    GaussResampler( double sigma= -1.0, BoundaryMethod boundary= Background,
                   T background= 0 )
    : Resampler<T>( boundary, background )
    {
      setSigma( sigma );
    }
    
    /** initialize the standard deviation of the Gaussian */
    inline void setSigma( double _sigma= -1.0 )
    {
      if( _sigma<= 0 )
        // least squares fit of a sinc
        sigma= sqrt( 3.0 * log( 2.0 ) ) / M_PI;
      else
        sigma= _sigma;
      
      multiplier= -1.0 / (2.0 * sigma * sigma);
    }
    
    /** linear resampling (oldX= a*newX + b) with Gaussian reconstruction*/
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** rational resampling (oldX = (a*newX + b) / (c*newX + d)),
     Gaussian reconstruction*/
    virtual void resampleIrregular( T *inLine, unsigned long inCount,
				    unsigned long inStride, T *outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
    
  protected:
    
    /** standard deviation */
    double sigma;
    
    /** multiplier used in the Gaussian exponent (i.e. -1/2/sigma^2) */
    double multiplier;
    
  };
  
  
  /** \class MinResampler Resampler.hh
   class for resampling using minimum of all affected samples */
  template<class T>
  class MinResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    MinResampler( BoundaryMethod boundary= Background,
                 T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with Gaussian reconstruction*/
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** rational resampling (oldX = (a*newX + b) / (c*newX + d)),
     Gaussian reconstruction*/
    virtual void resampleIrregular( T *inLine, unsigned long inCount,
				    unsigned long inStride, T *outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
    
  };
  
  
  
  /** \class MaxResampler Resampler.hh
   class for resampling using maximum of all affected samples */
  template<class T>
  class MaxResampler: public Resampler<T> {
    
  public:
    
    /** constructor */
    MaxResampler( BoundaryMethod boundary= Background,
                 T background= 0 )
    : Resampler<T>( boundary, background )
    {}
    
    /** linear resampling (oldX= a*newX + b) with Gaussian reconstruction*/
    virtual void resampleLinear( T *inLine, unsigned long inCount,
				 unsigned long inStride, T *outLine,
				 unsigned long outCount,
				 unsigned long outStride,
				 double a, double b );
    
    /** rational resampling (oldX = (a*newX + b) / (c*newX + d)),
     Gaussian reconstruction*/
    virtual void resampleIrregular( T *inLine, unsigned long inCount,
				    unsigned long inStride, T *outLine,
				    unsigned long outCount,
				    unsigned long outStride,
				    double *outSamples,
				    unsigned long sampleStride= 1ul );
    
  };
  
  
  
} /* namespace */




#endif /* RESAMPLING_RESAMPLER_H */
