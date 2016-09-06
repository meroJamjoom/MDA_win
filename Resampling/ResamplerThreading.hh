// ==========================================================================
// $Id: ResamplerThreading.hh 750 2010-09-04 02:14:32Z heidrich $
// threading support for resamplers
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

#ifndef RESAMPLING_RESAMPLERTHREADING_H
#define RESAMPLING_RESAMPLERTHREADING_H

/*! \file  ResamplerThreading.hh
    \brief threading support for resamplers
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Threading/SMPJob.hh"
#include "Resampler.hh"

namespace MDA {

  /** \class LinearResamplerJob Resampler.hh
      multithreading job for linear resampling of a data line */
  template <class T>
  class LinearResamplerJob: public SMPJob {
    
  public:
    
    /** constructor */
    LinearResamplerJob( Resampler<T> *_resampler, T* _inLine,
			unsigned long _inCount, unsigned long _inStride,
			T* _outLine, unsigned long _outCount,
			unsigned long _outStride, double _a, double _b )
      : resampler( _resampler ), inLine( _inLine ), inCount( _inCount ),
	inStride( _inStride ), outLine( _outLine ), outCount( _outCount ),
	outStride( _outStride ), a( _a ), b( _b )
    {}
    
    /** resample the line */
    virtual void execute( int threadID );
    
  public:
    
    /** the resampler */
    Resampler<T> *resampler;
    
    /** input line */
    T* inLine;

    /** input element count */
    unsigned long inCount;
    
    /** input stride */
    unsigned long inStride;
    
    /** output line */
    T *outLine;
    
    /** output element count */
    unsigned long outCount;
    
    /** output stride */
    unsigned long outStride;
    
    /** linear scaling parameter */
    double a;
    
    /** linear offset parameter */
    double b;
  };
  
  
  /** \class RationalResamplerJob Resampler.hh
   *  multithreading job for rational resampling of a data line
   *  (for convenience, this is just subclassed from the linear case)
   */
  template <class T>
  class RationalResamplerJob: public LinearResamplerJob<T> {
    
  public:
    
    /** constructor */
    RationalResamplerJob( Resampler<T> *_resampler, T* _inLine,
			  unsigned long _inCount, unsigned long _inStride,
			  T* _outLine, unsigned long _outCount,
			  unsigned long _outStride,
			  double _a, double _b, double _c, double _d )
      : LinearResamplerJob<T>( _resampler, _inLine, _inCount, _inStride,
			       _outLine, _outCount, _outStride, _a, _b ),
	c( _c ), d( _d )
    {}
    
    /** resample the line */
    virtual void execute( int threadID );
    
  public:
    
    /** rational scaling parameter */
    double c;
    
    /** rational offset parameter */
    double d;
  };
  
  
  /** \class IrregularResamplerJob Resampler.hh
   *  multithreading job for rational resampling of a data line
   *  (for convenience, this is just subclassed from the linear case)
   */
  template <class T>
  class IrregularResamplerJob: public LinearResamplerJob<T> {
    
  public:
    
    /** constructor */
    IrregularResamplerJob( Resampler<T> *_resampler, T* _inLine,
			   unsigned long _inCount, unsigned long _inStride,
			   T* _outLine, unsigned long _outCount,
			   unsigned long _outStride, double *_samples,
			   unsigned long _sampleStride= 1ul )
      : LinearResamplerJob<T>( _resampler, _inLine, _inCount, _inStride,
			       _outLine, _outCount, _outStride, 0.0, 0.0 ),
	samples( _samples ), sampleStride( _sampleStride )
    {}
    
    /** resample the line */
    virtual void execute( int threadID );
    
  public:
    
    /** rational scaling parameter */
    double* samples;
    
    /** rational offset parameter */
    unsigned long sampleStride;
  };
  
  
  
} /* namespace */

#endif /* RESAMPLING_RESAMPLERTHREADING_H */

