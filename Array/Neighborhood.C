// ==========================================================================
// $Id: Neighborhood.C 931 2012-01-29 00:15:09Z heidrich $
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

#ifndef ARRAY_NEIGHBORHOOD_C
#define ARRAY_NEIGHBORHOOD_C

#include <limits.h>
#include <math.h>

#include <MDA/Config.hh>
#include "Neighborhood.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

//
// local helper function for computing the p-norm of a
// CoordinateVector this should really be done using Vector norms, but
// I don't want the basic Array to be dependent on the LinearAlgebra
// subsystem
//
double
norm( CoordinateVector x, double p )
{
  unsigned long dimension= x.vec.size();
  unsigned long i;
  double dist;
  
  // really small p -> Lp is the minimum (semi)norm
  if( p< LP_NORM_THRESHOLD )
  {
    dist= DBL_MAX;
    for( i= 0 ; i< dimension ; i++ )
      dist= dist< fabs( (float)x.vec[i] ) ? dist : fabs( (float)x.vec[i] );
    return dist;
  }
  // really large p -> Lp is the maximum (Linfinity) norm
  if( p> 1.0/LP_NORM_THRESHOLD )
  {
    dist= 0.0;
    for( i= 0 ; i< dimension ; i++ )
      dist= dist> fabs( (float)x.vec[i] ) ? dist : fabs( (float)x.vec[i] );
    return dist;
  }
  
  // everything else - use the actual Lp formula
  dist= 0.0;
  for( i= 0 ; i< dimension ; i++ )
    dist+= pow( fabs( (float)x.vec[i] ), p );
  return pow( dist, 1/p );
}


//
// Neighborhood methods
//

/** hyercube constructor using radius */
void
Neighborhood::init( CoordinateVector &dim, unsigned long padding,
		    unsigned long radius, double norm )
{
  // set up a coordinate iterator from -radius...+radius in each direction
#if defined(_WIN32) || defined(_WIN64)
  LongRange range(-(static_cast<long>(radius)), radius );
#else
  LongRange range(-radius, radius);
#endif
  LongRangeList box;
  unsigned long dimension= dim.vec.size();
  for( unsigned long i= 0 ; i< dimension ; i++ )
    box.vec.push_back( range );
  
  init( dim, padding, box, norm );
}


/** hyercube constructor using radius and list of axes */
void
Neighborhood::init( CoordinateVector &dim, unsigned long padding,
		    AxisList &axes, unsigned long radius, double norm )
{
  // set up a coordinate iterator from -radius...+radius in each
  // direction contained in axes
#if defined(_WIN32) || defined(_WIN64)
	LongRange range(-(static_cast<long>(radius)), radius);
#else
	LongRange range(-radius, radius);
#endif
  LongRange empty( 0, 0 );
  LongRangeList box;
  unsigned long dimension= dim.vec.size();
  for( unsigned long i= 0 ; i< dimension ; i++ )
  {
    bool found= false;
    
    for( unsigned long j= axes.vec.size() ; j> 0 ; )
      if( axes.vec[--j]== i )
      {
	box.vec.push_back( range );
	found= true;
	break;
      }
    
    if( !found )
      box.vec.push_back( empty );
  }
  
  init( dim, padding, box, norm );
}


/** actual initialization method used by constructors */
void
Neighborhood::init( CoordinateVector &dim, unsigned long padding,
		    LongRangeList &box, double p )
{
  unsigned long i;
  
  // compute size of neighborhood, and allocate space for full box
  unsigned long dimension= dim.vec.size();
  numPixels= 1;
  for( i= 0 ; i< dimension ; i++ )
    numPixels*= box.vec[i].val.second - box.vec[i].val.first + 1;
  pixels= new Pixel[numPixels];
  
  // traverse the box
  CoordinateVectorIter iter( box );
  for( numPixels= 0, iter.begin() ; !iter.isAtEnd() ; ++iter )
  {
    CoordinateVector &pos= iter.getPos();
    
    // skip this position if it is rejected by the masking function
    if( !this->mask( pos ) )
      continue;
    
    // compute offsets for current position
    pixels[numPixels].uOff= pos.vec[dimension-1];
    pixels[numPixels].pOff= pos.vec[dimension-1];
    for( i= dimension-1 ; i> 0 ; )
    {
      i--;
      pixels[numPixels].uOff=
	pos.vec[i] + dim.vec[i]*pixels[numPixels].uOff;
      pixels[numPixels].pOff=
	pos.vec[i] + (dim.vec[i]+2*padding)*pixels[numPixels].pOff;
    }
    
    // compute distance for current position
    pixels[numPixels].dist= norm( pos, p );
    
    // update count of pixels within the mask
    numPixels++;
  }
}

/** mask out certain pixels in the box during neighborhood
    construction. True means the pixel is part of the
    neighborhood, false that it gets removed.  This method should
    be overloaded in subclasses for specific behaviors. */
bool
Neighborhood::mask( CoordinateVector &pos ) const
{
  return true;
}


//
// BallNeighborhood methods
//

/** mask out regions outside the ball */
bool
BallNeighborhood::mask( CoordinateVector &pos ) const
{
  return norm( pos, pNorm )<= ballRadius;
}


//
// StarNeighborhood methods
//

/** mask out regions outside the ball */
bool
StarNeighborhood::mask( CoordinateVector &pos ) const
{
  unsigned long count= 0;
  for( unsigned long i= pos.vec.size() ; i> 0 ; )
    if( pos.vec[--i]!= 0 && ++count> 1 )
      return false;
  return true;
}


} /* namespace */

#endif /* ARRAY_NEIGHBORHOOD_C */

