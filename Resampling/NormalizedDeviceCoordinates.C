// ==========================================================================
// $Id: NormalizedDeviceCoordinates.C 676 2010-04-04 00:03:47Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef WARPS_NORMALIZEDDEVICECOORDINATES_C
#define WARPS_NORMALIZEDDEVICECOORDINATES_C

#include "NormalizedDeviceCoordinates.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** constructor */
  NormalizedDeviceCoordinates::NormalizedDeviceCoordinates( const CoordinateVector &d, int aspect, bool yFlip )
  : dim( d ), dimension( dim.vec.size() )
{
  //
  // create matrices and vectors
  //

	  // NDC to pixel
  ndc2pixelH= Matrix( dimension+1, dimension+1 );
  ndc2pixelH.identity();
  ndc2pixelScale= ndc2pixelH.getDiagonal();
  ndc2pixelOff= ndc2pixelH.getColumnVector( dimension );

  // pixel to NDC
  pixel2ndcH= Matrix( dimension+1, dimension+1 );
  pixel2ndcH.identity();
  pixel2ndcScale= pixel2ndcH.getDiagonal();
  pixel2ndcOff= pixel2ndcH.getColumnVector( dimension );
  
  //
  // fill in the parameters
  //
  for( unsigned i= 0 ; i< dimension ; i++ )
  {
    // axis with respect to which we normalize
    unsigned normAxis= aspect< 0 ? i : aspect;
    
    pixel2ndcScale[i]= 2.0 / dim.vec[normAxis];
    pixel2ndcOff[i]= -(double)(dim.vec[i]-1)/(double)dim.vec[normAxis];
      
    ndc2pixelScale[i]= dim.vec[normAxis] / 2.0;
    ndc2pixelOff[i]= (dim.vec[i]-1) / 2.0;
  }

  // flip y if desired
  if( yFlip && dimension> 1 )
  {
    pixel2ndcScale[1]*=	-1.0;
    pixel2ndcOff[1]*=	-1.0;
    ndc2pixelScale[1]*=	-1.0;
  }

}



} /* namespace */

#endif /* WARPS_NORMALIZEDDEVICECOORDINATES_C */

