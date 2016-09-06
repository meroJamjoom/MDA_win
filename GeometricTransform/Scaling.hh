// ==========================================================================
// $Id: Scaling.hh 642 2010-03-09 20:10:31Z heidrich $
// Simple non-uniform scaling
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

#ifndef GEOMETRICTRANSFORM_SCALING_H
#define GEOMETRICTRANSFORM_SCALING_H

/*! \file  Scaling.hh
    \brief Simple non-uniform scaling
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/LinearAlgebra/LinAlg.hh"
#include "GeometricTransformation.hh"

namespace MDA {

  /** \class Scaling Scaling.hh
      Simple non-uniform scaling */
  template<class T>
  class Scaling: public GeometricTransformation<T> {

  public:

    /** default constructor */
    inline Scaling( Resampler<T> *_resampler= NULL )
      : GeometricTransformation<T>( _resampler )
    {}
    
    /** set the scaling parameters */
    inline void setScaling( const Vector &scale )
    {
      unsigned dimension= scale.getSize();
      scaling= Matrix( dimension+1, dimension+1 );
      scaling.identity();
      for( int i= 0 ; i< dimension ; i++ )
	// one over, since our matrix represents inverse mapping
	scaling[i][i]= 1.0/scale[i];
    }
    
    /** apply the transform to a single channel in one array */
    virtual bool apply( Array<T> &srcArray, int srcChannel,
			Array<T> &dstArray, int dstChannel );
    
  protected:
    
    /** homogeneous scaling matrix */
    Matrix scaling;
    
  };


} /* namespace */


#endif /* GEOMETRICTRANSFORM_SCALING_H */

