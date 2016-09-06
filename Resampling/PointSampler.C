// ==========================================================================
// $Id: PointSampler.C 735 2010-08-19 02:00:52Z heidrich $
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

#ifndef RESAMPLING_POINTSAMPLER_C
#define RESAMPLING_POINTSAMPLER_C

#include "PointSampler.hh"

#if defined (_WIN32) || defined (_WIN64)
#define _USE_MATH_DEFINES
#include <math.h>
#endif



namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


  //
  // PointSampler methods
  //


/** default constructor */
template<class T>
PointSampler<T>::PointSampler( CoordinateVector &arrayDimension,
			       unsigned long supportWidth )
  : dim( arrayDimension ), supWidth( supportWidth ),
    currentPos( arrayDimension.vec.size() )
{
  unsigned long i;
  
  // initialize dimensions and position vector
  dimension= dim.vec.size();
  
  // calculate # pixels in support
  supSize= 1;
  for( i= 0 ; i< dimension ; i++ )
    supSize*= supWidth;
  
  // initialize support offsets and weights (to be overridden whenever
  // a sample is placed)
  supOffsets= new long[supSize];
  supWeights= new T[supSize];
  for( i= 0 ; i< supSize ; i++ )
  {
    supOffsets[i]= -1;
    supWeights[i]= 0.0;
  }
}


/** default constructor */
template<class T>
PointSampler<T>::~PointSampler()
{
  delete [] supWeights;
  delete [] supOffsets;
}


/** calculate the array offset for a given integer position
    (including boundary effects) */
template<class T>
unsigned long
PointSampler<T>::calculateOffset( CoordinateVector &pos,
				  BoundaryMethod boundary )
{
  unsigned long offset= 0;
  long tmp;
  
  for( int i= dimension ; i> 0 ; )
  {
    offset*= dim.vec[--i];

    if( pos.vec[i]>= 0 && pos.vec[i]< dim.vec[i] )
      // interior of the array; straightforward
      offset+= pos.vec[i];
    else
    {
      // nasty boundary stuff
      switch( boundary )
      {
      case Background:
	// return index for boundary color
	return -1;
      case Renormalize:
      case Clamp:
	if( pos.vec[i]>= dim.vec[i] )
	  // "right" boundary: clip to max value
	  offset+= dim.vec[i]-1;
	// else clip to 0 (no change in offset)
	// i.e. offset+= 0;
	break;
      case Cyclic:
	if( pos.vec[i]< 0 )
	  // "left" boundary (the extra addition is because % is negative)
	  offset+= pos.vec[i] % dim.vec[i] + dim.vec[i];
	else
	  // "right" boundary
	  offset+= pos.vec[i] % dim.vec[i];
	break;
      case Mirror:
	tmp= pos.vec[i];
	if( tmp< 0 )
	  tmp= tmp%(dim.vec[i]*2) + dim.vec[i]*2;
	tmp= tmp%(dim.vec[i]*2);
	if( tmp>= dim.vec[i] )
	  // "reverse" traversal direction
	  offset+= (2*dim.vec[i]-1 - tmp);
	else
	  // "forward" traversal direction
	  offset+= tmp;
      }
    }
  }
  
  return offset;
}


/** calculate the array offset for a given integer position
    (including boundary effects) */
template<class T>
void
PointSampler<T>::setSampleLocation( Vector pos,
				    BoundaryMethod boundary )
{
  unsigned i, j;
  double w1, w2;
  unsigned supRad= supWidth/2;
  
  // "upper left" corner
  for( i= 0 ; i< PointSampler<T>::dimension ; i++ )
    PointSampler<T>::currentPos.vec[i]= (long)(pos[i]) - supRad+1;
  // compute offsets
  for( i= 0 ; i< PointSampler<T>::supSize ; i++ )
  {
    PointSampler<T>::supOffsets[i]=
      PointSampler<T>::calculateOffset( PointSampler<T>::currentPos,
					boundary );
    // update position 
    for( j= 0 ; j< PointSampler<T>::dimension ; j++ )
    {
      PointSampler<T>::currentPos.vec[j]++;
      if( PointSampler<T>::currentPos.vec[j]> (long)(pos[j]) + supRad )
	PointSampler<T>::currentPos.vec[j]= (long)(pos[j]) - supRad+1;
      else
	break;
    }
  }
  
  // now compute weights using interp method from subclasses
  supWeights[0]= 1.0;
  unsigned n= 1;
  for( i= dimension ; i> 0 ; )
  {
    --i;
    n= this->interp( n, pos[i]-(long)(pos[i]) );
  }
}

  //
  // NearestNeighborPointSampler methods
  //


/** setting a sample in NN is a bit different from the other
    reconstuction filters... */
template<class T>
void
NearestNeighborPointSampler<T>::setSampleLocation( Vector pos,
						   BoundaryMethod boundary )
{
  for( unsigned i= 0 ; i< PointSampler<T>::dimension ; i++ )
    PointSampler<T>::currentPos.vec[i]= (int)(pos[i] + .5);
  PointSampler<T>::supOffsets[0]=
    PointSampler<T>::calculateOffset( PointSampler<T>::currentPos, boundary );
}



  //
  // LinearPointSampler methods
  //

template<class T>
unsigned
LinearPointSampler<T>::interp( unsigned numSamples, double x )
{
  double w= 1.0-x;
  for( unsigned i= numSamples ; i> 0 ; )
  {
    --i;
    PointSampler<T>::supWeights[2*i+1]= x * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[2*i]=   w * PointSampler<T>::supWeights[i];
  }
  
  return 2*numSamples;
}  



  //
  // CubicPointSampler methods
  //

template<class T>
unsigned
CubicPointSampler<T>::interp( unsigned numSamples, double x )
{
  double w0= ((-x + 2.0)*x - 1.0)*x;	// f(-1), f centered at 0 <= x < 1
  double w1= (x - 2.0)*x*x + 1.0;	// f( 0), f centered at 0 <= x < 1
  double w2= ((-x + 1.0)*x + 1.0)*x;	// f( 1), f centered at 0 <= x < 1
  double w3= (x - 1.0)*x*x;		// f( 2), f centered at 0 <= x < 1
  
  for( unsigned i= numSamples ; i> 0 ; )
  {
    --i;
    PointSampler<T>::supWeights[4*i+3]= w3 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i+2]= w2 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i+1]= w1 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i]=   w0 * PointSampler<T>::supWeights[i];
  }
  
  return 4*numSamples;
}



  //
  // GaussPointSampler methods
  //

/** The scale factor for sigma that best approximates (in a
    least-squares sense) a sinc with a Gaussian:
	sigma= sqrt( 3.0 * log( 2.0 ) ) / M_PI;
	factor= -1.0 / (2.0 * sigma * sigma);
    Essentially, sigma minimizes
	integral( (sinc(x)-G(x,sigma)^2, x=-infinity..infinity )
*/
template<class T>
const double GaussPointSampler<T>::scaleFactor= -M_PI*M_PI / (6.0 * log( 2.0 ));

template<class T>
unsigned
GaussPointSampler<T>::interp( unsigned numSamples, double x )
{
  double w0= exp( scaleFactor*(1+x)*(1+x) ); // G(-1), G centered at 0 <= x < 1
  double w1= exp( scaleFactor*x*x );         // G( 0), G centered at 0 <= x < 1
  double w2= exp( scaleFactor*(1-x)*(1-x) ); // G( 1), G centered at 0 <= x < 1
  double w3= exp( scaleFactor*(2-x)*(2-x) ); // G( 2), G centered at 0 <= x < 1
  double wSum= w0+w1+w2+w3;
  w0/= wSum; w1/= wSum; w2/= wSum; w3/= wSum;
  
  for( unsigned i= numSamples ; i> 0 ; )
  {
    --i;
    PointSampler<T>::supWeights[4*i+3]= w3 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i+2]= w2 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i+1]= w1 * PointSampler<T>::supWeights[i];
    PointSampler<T>::supWeights[4*i]=   w0 * PointSampler<T>::supWeights[i];
  }
  
  return 4*numSamples;
}  


// explicit template instation code

template class PointSampler<float>;
template class PointSampler<double>;
template class NearestNeighborPointSampler<float>;
template class NearestNeighborPointSampler<double>;
template class LinearPointSampler<float>;
template class LinearPointSampler<double>;
template class CubicPointSampler<float>;
template class CubicPointSampler<double>;
template class GaussPointSampler<float>;
template class GaussPointSampler<double>;


} /* namespace */

#endif /* RESAMPLING_POINTSAMPLER_C */

