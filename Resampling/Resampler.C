// ==========================================================================
// $Id: Resampler.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef RESAMPLING_RESAMPLER_C
#define RESAMPLING_RESAMPLER_C

#include <iostream>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "MDA/Config.hh"
#include "MDA/Base/Errors.hh"
#include "Resampler.hh"

namespace MDA {
  
  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;
  
  
  //
  // Resampler methods
  //
  
  
  /** get a pixel from the boundary */
  template<class T>
  T
  Resampler<T>::getBoundaryPixel( T *inLine, long inCount,
                                 long inStride, long pos )
  {
    // the following test may seem redundant since it deals with
    // interior pixels, but supporting pixels inside the inLine makes it
    // easier to code resampling filters with wide support.
    if( pos< 0 || pos>= inCount )
    {
      // special case: samples very far outside the center range, just
      // return the background to avoid singularities
      if( boundary== Background || pos< -10000*inCount || pos> 10000*inCount )
        return Resampler<T>::background;
      
      // most cases should end up here: handling of pixels outside the
      // boundaries
      switch( boundary )
      {
      case Clamp:
      case Renormalize:
        if( pos< 0 )
          return *inLine;
        else
          return inLine[(inCount-1) * inStride];
        break;
      case Cyclic:
        if( pos< 0 )
          return inLine[(pos%inCount+inCount) * inStride];
        else
          return inLine[(pos%inCount) * inStride];
        break;
      case Mirror:
        if( pos< 0 )
          pos= pos%(inCount*2) + inCount*2;
        pos= pos%(inCount*2);
        if( pos>= inCount )
          return inLine[(2*inCount-1 - pos) * inStride];
        else
          return inLine[pos * inStride];
      case Background:
	// case covered above; should never get here
       return Resampler<T>::background;
      }
    }
    
    // all thats left over are the interior pixels
    return inLine[pos * inStride];
  }
  
  
  /** resample a scanline/column with linear function oldX= a*newX + b
   
   (the default implementation is to use resampleIrregular, but
   subclasses can override with more efficient specialized
   implementations)
   */
  template <class T>
  void
  Resampler<T>::resampleLinear( T *inLine, unsigned long inCount,
				unsigned long inStride, T *outLine,
				unsigned long outCount, unsigned long outStride,
				double a, double b )
  {
    double *samples= new double[outCount];
    for( unsigned long i= 0ul ; i< outCount ; i++ )
      samples[i]= a*i + b;
    this->resampleIrregular( inLine, inCount, inStride,
                            outLine, outCount, outStride,
                            samples );
    delete [] samples;
  }
  
  
  /** resample a scanline/column with rational function
   oldX = (a*newX + b) / (c*newX + d)
   
   (the default implementation is to use resampleIrregular, but
   subclasses can override with more efficient specialized
   implementations)
   */
  template <class T>
  void
  Resampler<T>::resampleRational( T *inLine, unsigned long inCount,
				  unsigned long inStride, T *outLine,
				  unsigned long outCount,
				  unsigned long outStride,
				  double a, double b, double c, double d )
  {
    double *samples= new double[outCount];
    double denom;
    
    for( unsigned long i= 0ul ; i< outCount ; i++ )
    {
      denom= (c*i + d);
      if( fabs( denom )< NUM_ZERO_THRESHOLD )
	samples[i]= DBL_MAX;
      else
	samples[i]= (a*i + b) / denom;
    }
    this->resampleIrregular( inLine, inCount, inStride,
			     outLine, outCount, outStride,
			     samples );
    delete [] samples;
  }
  
  
  //
  // NearestNeighborResampler methods
  //
  
  /** linear resampling with NN filter */
  template<class T>
  void
  NearestNeighborResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                              unsigned long inStride,
                                              T *outLine, unsigned long outCount,
                                              unsigned long outStride,
                                              double a, double b )
  {
    double srcXContinuous;
    long srcX, dstX;
    
    for( dstX= 0, srcXContinuous= b+.5 ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride  )
    {
      srcX= (long)srcXContinuous;
      if( srcX< 0 || srcX>= inCount )
        *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                 inStride, srcX );
      else
        *outLine= inLine[srcX*inStride];
    }
  }
  
  
  /** irregular resampling with NN filter */
  template<class T>
  void
  NearestNeighborResampler<T>::resampleIrregular( T *inLine,
						  unsigned long inCount,
						  unsigned long inStride,
						  T *outLine,
						  unsigned long outCount,
						  unsigned long outStride,
						  double *samples,
						  unsigned long sampleStride )
  {
    long srcX, dstX;
    
    for( dstX= 0 ; dstX< outCount ;
	 dstX++, outLine+= outStride, samples+=sampleStride )
    {
      srcX= (long)(*samples + .5);
      if( srcX< 0 || srcX>= inCount || (*samples) == DBL_MAX)
        *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
						  inStride, srcX );
      else
        *outLine= inLine[srcX*inStride];
    }
  }
  
  
  
  //
  // BoxResampler methods
  //
  
  /** linear resampling with box filter */
  template<class T>
  void
  BoxResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                  unsigned long inStride, T *outLine,
                                  unsigned long outCount,
                                  unsigned long outStride,
                                  double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    unsigned long numSamples;
    double suppRadius= fabs( a )<= 1.001 ? 0.5005 : fabs( a ) / 2.0;
    
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      srcX= (long)ceil( srcXContinuous-suppRadius );
      // average over source pixel within support
      for( value= 0.0, numSamples= 0 ;
          srcX<= srcXContinuous + suppRadius ;
          srcX++, numSamples++ )
        if( srcX< 0 || srcX>= inCount )
          value+=Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                inStride, srcX );
        else
          value+= inLine[srcX*inStride];
      
      // if the math is right, this should always hold...
      errorCond( numSamples> 0,
                "  Internal error: number of samples per pixel not positive!" );
      
      *outLine= value / numSamples;
    }
  }
  
  
  /** rational resampling with box filter */
  template<class T>
  void
  BoxResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
				      unsigned long inStride, T *outLine,
				      unsigned long outCount,
				      unsigned long outStride,
				      double *samples,
				      unsigned long sampleStride )
  {
    double suppRadius, srcXContinuous, srcXContinuousOld;
    T value;
    long srcX, dstX;
    unsigned long numSamples;
    bool haveAverage= false;
    T average= 0.0;
    
    srcXContinuousOld= samples[1]; // at start use a right-sided difference
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      suppRadius= fabs( srcXContinuous - srcXContinuousOld ) / 2.0;
      // adjust support radius for magnification regions
      if( suppRadius< 0.5005 )
        suppRadius= 0.5005;
      
      // beyond a certain minification level, we just give up and use
      // the background color or average, depending on boundary mode...
      if( suppRadius > inCount/2 || 
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX)
      {
        if( Resampler<T>::boundary== Background )
          *outLine= Resampler<T>::background;
        else
        {
          // compute average over whole line if not already done
          if( !haveAverage )
          {
            for( unsigned long i= 0 ; i< inCount ; i++ )
              average+= inLine[i*inStride];
            average/= inStride;
            haveAverage= true;
          }
          *outLine= average;
        }
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        // average over pixels in support
        for( value= 0.0, numSamples= 0 ;
            srcX<= srcXContinuous + suppRadius ;
            srcX++, numSamples++ )
          if( srcX< 0 || srcX>= inCount )
            value+= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                   inStride, srcX );
          else
            value+= inLine[srcX*inStride];
        
        *outLine= value / numSamples;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  
  //
  // LinearResampler methods
  //
  
  
  /** linear resampling with linear (hat) filter */
  template<class T>
  void
  LinearResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                     unsigned long inStride, T *outLine,
                                     unsigned long outCount,
                                     unsigned long outStride,
                                     double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    double weight, weightSum;
    double suppRadius= fabs( a )< 1.0 ? 1.0 : fabs( a );
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      srcX= (long)ceil( srcXContinuous - suppRadius );
      // weighted average over support
      for( value= weightSum= 0.0 ; srcX<= srcXContinuous + suppRadius ; srcX++ )
      {
        weight= 1.0 - fabs( srcXContinuous-srcX ) / suppRadius;
        weightSum+= weight;
        if( srcX< 0 || srcX>= inCount )
          value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                          inStride, srcX );
        else
          value+= weight * inLine[srcX*inStride];
      }
      
      *outLine= value / weightSum;
    }
  }
  
  
  /** rational resampling with linear (hat) filter */
  template<class T>
  void
  LinearResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
					 unsigned long inStride, T *outLine,
					 unsigned long outCount,
					 unsigned long outStride,
					 double *samples,
					 unsigned long sampleStride )
  {
    double srcXContinuous, srcXContinuousOld, suppRadius;
    T value;
    double weight, weightSum;
    long srcX, dstX;
    bool haveAverage= false;
    T average= 0.0;
    
    srcXContinuousOld= samples[1]; // right-sided difference at left of array
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      suppRadius= fabs( srcXContinuous - srcXContinuousOld );
      // adjust support radius in magnification regions
      if( suppRadius< 1.0 )
        suppRadius= 1.0;
      
      // if support radius exceeds width of line, we just use average or
      // background color
      if( suppRadius> inCount ||
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX )
      {
        if( Resampler<T>::boundary== Background )
          *outLine= Resampler<T>::background;
        else
        {
          // compute average over whole line if not already done
          if( !haveAverage )
          {
            for( unsigned long i= 0 ; i< inCount ; i++ )
              average+= inLine[i*inStride];
            average/= inStride;
            haveAverage= true;
          }
          *outLine= average;
        }
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        // weighted average over pixels in support
        for( value= weight= weightSum= 0.0 ;
            srcX<= srcXContinuous + suppRadius ;
            srcX++ )
        {
          weight= 1.0 - fabs( srcXContinuous-srcX ) / suppRadius;
          weightSum+= weight;
          if( srcX< 0 || srcX>= inCount )
            value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                            inStride, srcX );
          else
            value+= weight * inLine[srcX*inStride];
        }
        *outLine= value / weightSum;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  
  //
  // CubicResampler methods
  //
  
  
  /** linear resampling with cubic filter */
  template<class T>
  void
  CubicResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                    unsigned long inStride, T *outLine,
                                    unsigned long outCount,
                                    unsigned long outStride,
                                    double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    double dist, weight, weightSum;
    double scale= fabs( a )< 1.0 ? 1.0 : fabs( a );
    double suppRadius= 2.0 * scale;
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      srcX= (long)ceil( srcXContinuous - suppRadius );
      // weighted average over support
      for( value= weightSum= 0.0 ; srcX<= srcXContinuous + suppRadius ; srcX++ )
      {
        dist= fabs( srcXContinuous-srcX ) / scale;
        if( dist<= 1.0 )
          weight= (dist - 2.0) * dist*dist + 1.0;
        else
          weight= ((5.0 - dist) * dist - 8.0) * dist + 4.0;
        weightSum+= weight;
        if( srcX< 0 || srcX>= inCount )
          value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                          inStride, srcX );
        else
          value+= weight * inLine[srcX*inStride];
      }
      
      *outLine= value / weightSum;
    }
  }
  
  
  /** rational resampling with cubic filter */
  template<class T>
  void
  CubicResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
					unsigned long inStride, T *outLine,
					unsigned long outCount,
					unsigned long outStride,
					double *samples,
					unsigned long sampleStride )
  {
    double srcXContinuous, srcXContinuousOld, suppRadius, scale;
    T value;
    double weight, weightSum, dist;
    long srcX, dstX;
    bool haveAverage= false;
    T average= 0.0;
    
    srcXContinuousOld= samples[1]; // right-sided difference for first sample
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      scale= fabs( srcXContinuous - srcXContinuousOld );
      // adjust scale in magnification regions
      if( scale< 1.0 )
        scale= 1.0;
      suppRadius= 2.0 * scale;
      
      // if scale exceeds width of line, we just use average or
      // background color
      if( scale> inCount ||
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX )
      {
        if( Resampler<T>::boundary== Background )
          *outLine= Resampler<T>::background;
        else
        {
          // compute average over whole line if not already done
          if( !haveAverage )
          {
            for( unsigned long i= 0 ; i< inCount ; i++ )
              average+= inLine[i*inStride];
            average/= inStride;
            haveAverage= true;
          }
          *outLine= average;
        }
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        // weighted average over pixels in support
        for( value= weight= weightSum= 0.0 ;
            srcX<= srcXContinuous + suppRadius ;
            srcX++ )
        {
          dist= fabs( srcXContinuous-srcX ) / scale;
          if( dist<= 1.0 )
            weight= (dist - 2.0) * dist*dist + 1.0;
          else
            weight= ((5.0 - dist) * dist - 8.0) * dist + 4.0;
          weightSum+= weight;
          if( srcX< 0 || srcX>= inCount )
            value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                            inStride, srcX );
          else
            value+= weight * inLine[srcX*inStride];
        }
        *outLine= value / weightSum;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  
  //
  // GaussResampler methods
  //
  
  /** linear resampling Gauss with filter */
  template<class T>
  void
  GaussResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                    unsigned long inStride, T *outLine,
                                    unsigned long outCount,
                                    unsigned long outStride,
                                    double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    double dist, weight, weightSum;
    double scale= fabs( a )< 1.0 ? 1.0 : fabs( a );
    double suppRadius= 3.0 * sigma * scale;
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      srcX= (long)ceil( srcXContinuous - suppRadius );
      // average over support
      for( value= weightSum= 0.0 ; srcX<= srcXContinuous + suppRadius ; srcX++ )
      {
        dist= (srcXContinuous - srcX) / scale;
        weight= exp( multiplier * dist*dist );
        weightSum+= weight;
        if( srcX< 0 || srcX>= inCount )
          value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                          inStride, srcX );
        else
          value+= weight * inLine[srcX*inStride];
      }
      
      *outLine= value / weightSum;
    }
  }
  
  
  /** rational resampling with filter in the minification case */
  template<class T>
  void
  GaussResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
					unsigned long inStride, T *outLine,
					unsigned long outCount,
					unsigned long outStride,
					double *samples,
					unsigned long sampleStride )
  {
    double srcXContinuous, srcXContinuousOld, suppRadius, scale, dist;
    T value;
    double weight, weightSum;
    long srcX, dstX;
    bool haveAverage= false;
    T average= 0.0;
    
    srcXContinuousOld= samples[1]; // right sided difference for first sample
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      scale= fabs( srcXContinuous - srcXContinuousOld );
      // adjust scale in magnification regions
      if( scale< 1.0 )
        scale= 1.0;
      suppRadius= 3.0 * sigma * scale;
      
      // if support radius exceeds width of line, we just use average /
      // background color
      if( scale> inCount ||
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX )
      {
        if( Resampler<T>::boundary== Background )
          *outLine= Resampler<T>::background;
        else
        {
          // compute average over whole line if not already done
          if( !haveAverage )
          {
            for( unsigned long i= 0 ; i< inCount ; i++ )
              average+= inLine[i*inStride];
            average/= inStride;
            haveAverage= true;
          }
          *outLine= average;
        }
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        // average over pixels in support with linear (i.e. hat) weight function
        for( value= weight= weightSum= 0.0 ;
            srcX<= srcXContinuous + suppRadius ;
            srcX++ )
        {
          dist= (srcXContinuous - srcX) / scale;
          weight= exp( multiplier * dist*dist );
          weightSum+= weight;
          if( srcX< 0 || srcX>= inCount )
            value+= weight * Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                            inStride, srcX );
          else
            value+= weight * inLine[srcX*inStride];
        }
        *outLine= value / weightSum;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  //
  // MinResampler methods
  //
  
  /** linear resampling with a minumum filter */
  template<class T>
  void
  MinResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                  unsigned long inStride, T *outLine,
                                  unsigned long outCount,
                                  unsigned long outStride,
                                  double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    double suppRadius= fabs( a )<= 1.001 ? 0.5005 : fabs( a ) / 2.0;
    
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      // assign value of first spurce pixel within support
      srcX= (long)ceil( srcXContinuous-suppRadius );
      if( srcX< 0 || srcX>= inCount )
        *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                 inStride, srcX );
      else
        *outLine= inLine[srcX*inStride];
      
      // find min over source pixels within support
      for( srcX++ ; srcX<= srcXContinuous + suppRadius ; srcX++ )
      {
        if( srcX< 0 || srcX>= inCount )
          value= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                inStride, srcX );
        else
          value= inLine[srcX*inStride];
        if( value< *outLine )
          *outLine= value;
      }
    }
  }
  
  
  /** rational resampling with box filter */
  template<class T>
  void
  MinResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
				      unsigned long inStride, T *outLine,
				      unsigned long outCount,
				      unsigned long outStride,
				      double *samples,
				      unsigned long sampleStride )
  {
    double srcXContinuous, srcXContinuousOld, suppRadius;
    T value;
    long srcX, dstX;
    bool haveLineMin= false;
    T lineMin= 0.0;
    
    srcXContinuousOld= samples[1]; // right sided difference for first sample
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      suppRadius= fabs( srcXContinuous - srcXContinuousOld ) / 2.0;
      // adjust support radius for magnification regions
      if( suppRadius< 0.5005 )
        suppRadius= 0.5005;
      
      // beyond a certain minification level, we just give up and use
      // the background color or min value, depending on boundary mode...
      if( suppRadius > inCount/2 ||
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX )
      {
        if( !haveLineMin )
        {
          if( Resampler<T>::boundary== Background )
            lineMin= Resampler<T>::background;
          else
            lineMin= inLine[0];
          
          // compute min over whole line if not already done
          for( unsigned long i= 0 ; i< inCount ; i++ )
            if( lineMin> inLine[i*inStride] )
              lineMin= inLine[i*inStride];
          haveLineMin= true;
        }
        *outLine= lineMin;
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        if( srcX< 0 || srcX>= inCount )
          *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                   inStride, srcX );
        else
          *outLine= inLine[srcX*inStride];
        
        // min over pixels in support
        for( srcX++ ; srcX<= srcXContinuous + suppRadius ; srcX++ )
          if( srcX< 0 || srcX>= inCount )
            value= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                  inStride, srcX );
          else
            value= inLine[srcX*inStride];
        if( value< *outLine )
          *outLine= value;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  //
  // MaxResampler methods
  //
  
  /** linear resampling with a minumum filter */
  template<class T>
  void
  MaxResampler<T>::resampleLinear( T *inLine, unsigned long inCount,
                                  unsigned long inStride, T *outLine,
                                  unsigned long outCount,
                                  unsigned long outStride,
                                  double a, double b )
  {
    double srcXContinuous;
    T value;
    long srcX, dstX;
    double suppRadius= fabs( a )<= 1.001 ? 0.5005 : fabs( a ) / 2.0;
    
    
    for( dstX= 0, srcXContinuous= b ;
        dstX< outCount ;
        dstX++, srcXContinuous+= a, outLine+= outStride )
    {
      // assign value of first spurce pixel within support
      srcX= (long)ceil( srcXContinuous-suppRadius );
      if( srcX< 0 || srcX>= inCount )
        *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                 inStride, srcX );
      else
        *outLine= inLine[srcX*inStride];
      
      // find min over source pixels within support
      for( srcX++ ; srcX<= srcXContinuous + suppRadius ; srcX++ )
      {
        if( srcX< 0 || srcX>= inCount )
          value= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                inStride, srcX );
        else
          value= inLine[srcX*inStride];
        if( value> *outLine )
          *outLine= value;
      }
    }
  }
  
  
  /** rational resampling with box filter */
  template<class T>
  void
  MaxResampler<T>::resampleIrregular( T *inLine, unsigned long inCount,
				      unsigned long inStride, T *outLine,
				      unsigned long outCount,
				      unsigned long outStride,
				      double *samples,
				      unsigned long sampleStride )
  {
    double srcXContinuous, srcXContinuousOld, suppRadius;
    T value;
    long srcX, dstX;
    bool haveLineMax= false;
    T lineMax= 0.0;
    
    srcXContinuousOld= samples[1]; // right sided difference for first sample
    for( dstX= 0 ; dstX< outCount ; dstX++, outLine+= outStride )
    {
      srcXContinuous= samples[dstX*sampleStride];
      suppRadius= fabs( srcXContinuous - srcXContinuousOld ) / 2.0;
      // adjust support radius for magnification regions
      if( suppRadius< 0.5005 )
        suppRadius= 0.5005;
      
      // beyond a certain minification level, we just give up and use
      // the background color or min value, depending on boundary mode...
      if( suppRadius > inCount/2 ||
         srcXContinuous == DBL_MAX || srcXContinuousOld == DBL_MAX )
      {
        if( !haveLineMax )
        {
          if( Resampler<T>::boundary== Background )
            lineMax= Resampler<T>::background;
          else
            lineMax= inLine[0];
          
          // compute min over whole line if not already done
          for( unsigned long i= 0 ; i< inCount ; i++ )
            if( lineMax> inLine[i*inStride] )
              lineMax= inLine[i*inStride];
          haveLineMax= true;
        }
        *outLine= lineMax;
      }
      else
      {
        srcX= (long)ceil( srcXContinuous - suppRadius );
        if( srcX< 0 || srcX>= inCount )
          *outLine= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                   inStride, srcX );
        else
          *outLine= inLine[srcX*inStride];
        
        // min over pixels in support
        for( srcX++ ; srcX<= srcXContinuous + suppRadius ; srcX++ )
          if( srcX< 0 || srcX>= inCount )
            value= Resampler<T>::getBoundaryPixel( inLine, inCount,
                                                  inStride, srcX );
          else
            value= inLine[srcX*inStride];
        if( value> *outLine )
          *outLine= value;
      }
      
      srcXContinuousOld= srcXContinuous;
    }
  }
  
  
  
  
  
  //
  // template instantiation code
  //
  template class Resampler<unsigned char>;
  template class NearestNeighborResampler<unsigned char>;
  template class BoxResampler<unsigned char>;
  template class LinearResampler<unsigned char>;
  template class CubicResampler<unsigned char>;
  template class GaussResampler<unsigned char>;
  template class MinResampler<unsigned char>;
  template class MaxResampler<unsigned char>;
  
  
  template class Resampler<float>;
  template class Resampler<double>;
  template class NearestNeighborResampler<float>;
  template class NearestNeighborResampler<double>;
  template class BoxResampler<float>;
  template class BoxResampler<double>;
  template class LinearResampler<float>;
  template class LinearResampler<double>;
  template class CubicResampler<float>;
  template class CubicResampler<double>;
  template class GaussResampler<float>;
  template class GaussResampler<double>;
  template class MinResampler<float>;
  template class MinResampler<double>;
  template class MaxResampler<float>;
  template class MaxResampler<double>;
  
  
} /* namespace */

#endif /* RESAMPLING_RESAMPLER_C */
