// ==========================================================================
// $Id: NormalizedDeviceCoordinates.hh 676 2010-04-04 00:03:47Z heidrich $
// mapping between NDC and pixel indices
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

#ifndef WARPS_NORMALIZEDDEVICECOORDINATES_H
#define WARPS_NORMALIZEDDEVICECOORDINATES_H

/*! \file  NormalizedDeviceCoordinates.hh
    \brief mapping between NDC and pixel indices
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/LinearAlgebra/LinAlg.hh>

namespace MDA {

  /** \class NormalizedDeviceCoordinates NormalizedDeviceCoordinates.hh
      mapping between NDC and pixel indices.
  
      The mapping makes sure that images tile seamlessly, so the
      max/min coordinates are usually NOT exactly -1 and 1.
  */
  
  class NormalizedDeviceCoordinates {
    
  public:
    
    /** constructor
	\param d: image dimensions
	\param aspect: mode for dealing w/ aspect ratios. Negative
	numbers: ignore aspect ratios (true NDC mapping, -1..1 for all
	aex). Non-negative numbers represent the axis with respect to
	which the others are scaled to give "square" pixels
	(NDC+viewport mapping)
	\param yFlip: whether to flip along axis 1 to map from typical
	image coordinates (origin in top left) to normal RHS
	coordinates (origin bottom left)*/
    NormalizedDeviceCoordinates( const CoordinateVector &d,
				 int aspect= -1, bool yFlip= true );
    
    //
    // NDC to pixel
    //
    
    /** return the (homogeneous) NDC-to-pixel matrix */
    inline Matrix &getNDCToPixelMatrix()
    {
      return ndc2pixelH;
    }

    /** map point from NDC to pixels (in-place) */
    inline Vector &pointFromNDCToPixel( Vector &pt )
    {
      return pointFromNDCToPixel( pt, pt );
    }
    
    /** map point from NDC to pixels (separate result vector) */
    inline Vector &pointFromNDCToPixel( Vector &pt, Vector &res )
    {
      for( unsigned i= 0 ; i< dimension ; i++ )
	res[i]= pt[i]*ndc2pixelScale[i] + ndc2pixelOff[i];
      return res;
    }
    
    /** map single point component from NDC to pixels */
    inline double pointCompFromNDCToPixel( double comp, unsigned which )
    {
      return comp*ndc2pixelScale[which] + ndc2pixelOff[which];
    }
    
    /** map vector from NDC to pixels (in place) */
    inline Vector &vecFromNDCToPixel( Vector &vec )
    {
      return vecFromNDCToPixel( vec, vec );
    }
    
    /** map vector from NDC to pixels (separate result vector) */
    inline Vector &vecFromNDCToPixel( Vector &vec, Vector &res )
    {
      for( unsigned i= 0 ; i< dimension ; i++ )
	res[i]= vec[i]*ndc2pixelScale[i];
      return res;
    }
    
    /** map single vector component from NDC to pixels */
    inline double vecCompFromNDCToPixel( double comp, unsigned which )
    {
      return comp*ndc2pixelScale[which];
    }
    
    //
    // pixel to NDC
    //
    
    /** return the (homogeneous) pixel-to-NDC matrix */
    inline Matrix &getPixelToNDCMatrix()
    {
      return pixel2ndcH;
    }
    
    /** map point from pixels to NDC (in place) */
    inline Vector &pointFromPixelToNDC( Vector &pt )
    {
      return pointFromPixelToNDC( pt, pt );
    }
    
    /** map point from pixels to NDC (separate result vector) */
    inline Vector &pointFromPixelToNDC( Vector &pt, Vector &res )
    {
      for( unsigned i= 0 ; i< dimension ; i++ )
	res[i]= pt[i]*pixel2ndcScale[i] + pixel2ndcOff[i];
      return res;
    }
    
    /** map single point component from pixels to NDC */
    inline double pointCompFromPixelToNDC( double comp, unsigned which )
    {
      return comp*pixel2ndcScale[which] + pixel2ndcOff[which];
    }
    
    /** map vector from pixels to NDC (in place) */
    inline Vector &vecFromPixelToNDC( Vector &vec )
    {
      return vecFromPixelToNDC( vec, vec );
    }
    
    /** map vector from pixels to NDC (separate result vector) */
    inline Vector &vecFromPixelToNDC( Vector &vec, Vector &res )
    {
      for( unsigned i= 0 ; i< dimension ; i++ )
	res[i]= vec[i]*pixel2ndcScale[i];
      return res;
    }
    
    /** map single vector component from pixels to NDC */
    inline double vecCompFromPixelToNDC( double comp, unsigned which )
    {
      return comp*pixel2ndcScale[which];
    }
    
  protected:
    
    // array size
    CoordinateVector dim;
    
    // dimensionality of array
    unsigned dimension;
    
    // mapping from NDC to pixels (full homogeneous version)
    Matrix ndc2pixelH;
    
    // scale from NDC to pixels (diagonal of ndc2pixelH)
    Vector ndc2pixelScale;

    // mapping from NDC to pixels (translation part)
    Vector ndc2pixelOff;
    
    // mapping from pixels to NDC (full homogeneous version)
    Matrix pixel2ndcH;
    
    // mapping from pixels to NDC (diagonal of pixel2ndcH)
    Vector pixel2ndcScale;
    
    // mapping from pixels to NDC (translation part)
    Vector pixel2ndcOff;
  };


} /* namespace */

#endif /* WARPS_NORMALIZEDDEVICECOORDINATES_H */

