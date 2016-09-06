// ==========================================================================
// $Id: Linear1DFilter.C 638 2010-03-08 19:34:24Z heidrich $
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

#ifndef FILTERS_LINEAR1DFILTER_C
#define FILTERS_LINEAR1DFILTER_C

#include <math.h>

#include "MDA/Config.hh"
#include "MDA/Base/Errors.hh"
#include "MDA/Array/Boundary.hh"

#include "Linear1DFilter.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


//
// Linear1DFilter Code
//

/** apply filter to one channel of an array */
template<class T>
bool
Linear1DFilter<T>::apply( Array<T> &array,
			  BoundaryMethod boundary,
			  unsigned int axis,
			  unsigned int inChannel,
			  unsigned int outChannel )
{
  // check boundary mode
  if( boundary== Renormalize )
  {
    double integral= 0.0;
    for( unsigned i= 2*radius+1 ; i> 0 ; )
      integral+= filter[--i];
    if( !warnCond( fabs( 1.0-integral )> FILTER_NORM_THRESHOLD,
		   "  renormalize mode requires a filter integrating to 1!\n"
		   "  (falling back to clamp mode)" ) )
      boundary= Clamp;
  }
  
  // call superclass method for actually applying the filter
  return Filter1D<T>::apply( array, boundary, axis, inChannel, outChannel );
}
  
/** provides an estimate of the computational affort involved in
    processing a certain number of elements (i.e #elem*(radius*2+1)) */
template <class T>
double
Linear1DFilter<T>::getLineCost( unsigned long numElements )
{
  return numElements * (radius*2+1);
}

/** apply filter to a single line in the array */
template <class T>
void
Linear1DFilter<T>::apply( T *startPos, unsigned long incr,
			  unsigned long numElements, BoundaryMethod boundary,
			  T background, T *startPosOut )
{
  unsigned long j, k, l;
  
  T* lineBuf= new T[numElements+2*radius];
  
  // fetch line into temp buffer
  fetchLine<T>( lineBuf, startPos, incr, numElements, radius,
		boundary, background );
  
  // start and stop are positions along the line - set these to
  // cover full line (overruled by renormalize mode)
  unsigned long start= 0;
  unsigned long stop= numElements;
  
  // in renormalize mode, we compute the boundary pixels separately
  if( boundary== Renormalize )
  {
    for( j= 0 ; j< radius ; j++ )
    {
      double accum, weight;
      
      // pixel at front of line
      accum= weight= 0.0;
      for( k= radius, l= radius-j ; k<= j+2*radius ; k++, l++ )
      {
	accum+= lineBuf[k] * filter[l];
	weight+= filter[l];
      }
      startPosOut[j*incr]= accum/weight;
      
      // pixel at back of line
      accum= weight= 0.0;
      for( k= numElements-1-j, l= 0 ; k<= numElements-1+radius ; k++, l++ )
      {
	accum+= lineBuf[k] * filter[l];
	weight+= filter[l];
      }
      startPosOut[(numElements-1-j)*incr]= accum/weight;
    }
    
    // finally, start and stop need to be adjusted for renormalize
    // mode to exclude the pixels we just processed
    start+= radius;
    stop-= radius;
  }
    
  // apply filter to center pixels in renormalize mode, all pixels
  // in all other modes:
  for( j= start ; j< stop ; j++ )
  {
    double h= 0.0;
    for( k= j, l= 0 ; k<= j+2*radius ; k++, l++ )
      h+= filter[l] * lineBuf[k];
    startPosOut[j*incr]= h;
  }

  // clean up
  delete [] lineBuf;
}
    


//
// Individual filter class constructors
//


/** BoxFilter1D constructor */
template<class T>
BoxFilter1D<T>::BoxFilter1D( unsigned int _radius, double multiplier )
  : Linear1DFilter<T>( _radius )
{
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;
  // normalization constant
  double val= multiplier / (double)(2*rad+1);
  
  for( int i= 0 ; i< 2*rad+1 ; i++ )
    f[i]= val;
}


/** HatFilter1D constructor */
template<class T>
HatFilter1D<T>::HatFilter1D( unsigned int _radius, double multiplier )
  : Linear1DFilter<T>( _radius )
{
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;
  // this normalizes the integral
  multiplier/= (1.0 + Linear1DFilter<T>::radius);
  
  for( int i= 0 ; i< rad ; i++ )
    f[i]= f[2*rad-i]= multiplier * (double)(i+1)/(double)(rad+1);
  f[rad]= multiplier;
}


/** Gaussian1D constructor */
template<class T>
Gaussian1D<T>::Gaussian1D( double sigma, int _radius,
			   double multiplier )
  : Linear1DFilter<T>( _radius<= 0 ? (int)(2.0*sigma+.5) : _radius )
{
  int i;
  double integral= 0.0;
  
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;
  
  for( i= 0 ; i< rad ; i++ )
    integral+= f[i]= f[2*rad-i]=
      exp( (double)((rad-i)*(rad-i))/(-2.0*sigma*sigma) );
  f[rad]= 1.0;
  // normalize
  integral= 2.0*integral + 1.0;
  for( i= 0 ; i<= 2*rad ; i++ )
    f[i]*= multiplier / integral;
}


/** constructor */
template<class T>
FastGaussian1D<T>::FastGaussian1D( double _sigma, double multiplier )
  : Linear1DFilter<T>( 2 ), sigma( _sigma )
{
  double q;
  
  // choose pseudo- standard deviation
  // constants are from pg. 144 of Young & van Vliet (eq. 11b)
  if( sigma>= 2.5 )
    q= 0.98711*sigma - 0.96330;
  else
    q= 3.97156 - 4.14554*sqrt( 1.0 - 0.26891*sigma );
  
  // simplify notation
  double *f= Linear1DFilter<T>::filter;
  double norm= multiplier / (1.57825 + q*(2.44413 + q*(1.4281 + q*0.422205)));
  f[1]= q*(2.44413 + q*(2.85619 + q*1.26661)) * norm;
  f[2]= -q*q*(1.4281 + q*1.26661) * norm;
  f[3]= q*q*q*0.422205 * norm;
  f[0]= multiplier - f[1] - f[2] - f[3];
}

    
/** apply filter to one channel of an array */
template<class T>
bool
FastGaussian1D<T>::apply( Array<T> &array,
			  BoundaryMethod boundary,
			  unsigned int axis,
			  unsigned int inChannel,
			  unsigned int outChannel )
{
  // check boundary mode
  if( !warnCond( boundary!= Renormalize,
		 "  FastGaussian does not support renormalize mode\n"
		 "  (falling back to clamp mode)" ) )
    boundary= Clamp;
  
  // then call super class method to del with the application of the filter
  return Filter1D<T>::apply( array, boundary, axis, inChannel, outChannel );
}

/** apply filter to a single line in the array */
template <class T>
void
FastGaussian1D<T>::apply( T *startPos, unsigned long incr,
			  unsigned long numElements, BoundaryMethod boundary,
			  T background, T *startPosOut )
{
  unsigned long j, k, l;
  
  // width of boundary region
  // experimentally, a factor of 3 is necessary to avoid artifacts.
  unsigned bWidth= (unsigned)(3.0*sigma);
  // allocate temporary memory for one line
  T* lineBuf= new T[numElements+2*bWidth];
  
  // notational convenience
  double* f= Linear1DFilter<T>::filter;
  
  // fetch line into temp buffer
  fetchLine<T>( lineBuf, startPos, incr, numElements, bWidth,
		boundary, background );
  
  // forward pass
  for( j= 3 ; j< numElements+2*bWidth ; j++ )
    lineBuf[j]=
      f[3] * lineBuf[j-3] +
      f[2] * lineBuf[j-2] +
      f[1] * lineBuf[j-1] +
      f[0] * lineBuf[j];
  
    // backward pass
  for( j= numElements+2*bWidth-3 ; j> 0 ; )
  {
    j--;
    lineBuf[j]=
      f[3] * lineBuf[j+3] +
      f[2] * lineBuf[j+2] +
      f[1] * lineBuf[j+1] +
      f[0] * lineBuf[j];
  }
    
  // write out result
  for( j= 0, k= bWidth, l= 0 ; j< numElements ; j++, k++, l+= incr )
    startPosOut[l]= lineBuf[k];
  
  // clean up
  delete [] lineBuf;
}
   


/** LaplacianOfGaussian1D constructor */
template<class T>
LaplacianOfGaussian1D<T>::LaplacianOfGaussian1D( double sigma,
						 int _radius,
						 double multiplier )
  : Linear1DFilter<T>( _radius<= 0 ? (int)(4.0*sigma+1.5) : _radius )
{
  int i;
  
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;
  double *tmp= new double[rad+1];
  
  // compute Gaussian first
  for( i= 0 ; i<= rad ; i++ )
    tmp[i]= exp( (double)((rad-i)*(rad-i))/(-2.0*sigma*sigma) );
  
  // now convolve w/ derivative filter
  f[0]= f[2*rad]= 0.0;
  f[rad]= 2.0 * (tmp[rad-1] - tmp[rad]);
  double integral= 0.0;
  for( i= 1 ; i< rad ; i++ )
    integral+= f[i]= f[2*rad-i]= tmp[i-1] + tmp[i+1] - 2.0*tmp[i];
  
  // finally, renormalize each half of the filter to 1
  for( i= 0 ; i<= 2*rad ; i++ )
    f[i]/= integral;
  
  delete [] tmp;
}


/** FirstDerivative1D constructor */
template<class T>
FirstDerivative1D<T>::FirstDerivative1D( unsigned int baseline )
  : Linear1DFilter<T>( baseline )
{
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;

  f[0]= -1.0;
  f[2*rad]= 1.0;
  
  for( int i= 1 ; i< 2*rad ; i++ )
    f[i]= 0.0;
}


/** SecondDerivative1D constructor */
template<class T>
SecondDerivative1D<T>::SecondDerivative1D( unsigned int baseline )
  : Linear1DFilter<T>( baseline )
{
  // simplify notation
  int rad= Linear1DFilter<T>::radius;
  double *f= Linear1DFilter<T>::filter;

  for( int i= 1 ; i< 2*rad ; i++ )
    f[i]= 0.0;
  f[0]= -1.0;
  f[rad]= 2.0;
  f[2*rad]= -1.0;
}



/* we use explicit instantiation for classes in this file */
template class Linear1DFilter<float>;
template class Linear1DFilter<double>;
template class BoxFilter1D<float>;
template class BoxFilter1D<double>;
template class HatFilter1D<float>;
template class HatFilter1D<double>;
template class Gaussian1D<float>;
template class Gaussian1D<double>;
template class FastGaussian1D<float>;
template class FastGaussian1D<double>;
template class LaplacianOfGaussian1D<float>;
template class LaplacianOfGaussian1D<double>;
template class FirstDerivative1D<float>;
template class FirstDerivative1D<double>;
template class SecondDerivative1D<float>;
template class SecondDerivative1D<double>;

} /* namespace */

#endif /* FILTERS_LINEAR1DFILTER_C */

