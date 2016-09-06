// ==========================================================================
// $Id: Camera.hh $
// camera information like intrinsic and extrinsic parameters
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: bradleyd ()
// Email:   bradleyd@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIB_CAMERA_H
#define CAMERACALIB_CAMERA_H

/*! \file  Camera.hh
    \brief camera information like intrinsic and extrinsic parameters
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/Base/MetaData.hh>
#include <MDA/LinearAlgebra/LinAlg.hh>

namespace MDA {
  
  /** \class Camera Camera.hh
      camera information like intrinsic and extrinsic parameters */
  
  class Camera {

  public:

    /** empty default constructor */
    inline Camera()
      : K(3,3),R(3,3),T(3),D(2)
    {}
    
    /** constructor using actual values */
    inline Camera( const Matrix &Kin, const Matrix &Rin, const Vector &Tin, const Vector &Din )
    {
      K.copy( Kin );
      R.copy( Rin );
      T.copy( Tin );
      D.copy( Din );
    }

    /** access method for principal point */
    inline Vector getPrincipalPoint() { return Vector(K[0][2], K[1][2]); }

    /** access method for focal length */
    inline Vector getFocalLength() { return Vector(K[0][0], K[1][1]); }
    
    /** Intrinsic parameters in NDC (3x3 matrix)
	
        Focal length in x and y is (K[0][0], K[1][1]).  
	Principal point in x and y is (K[0][2], K[1][2]).
     */
    Matrix K;
    
    /** Extrinsic parameters 
	
        R is the rotation component (3x3 matrix).  
	T is the translation component (3x1 vector).
     */
    Matrix R;
    Vector T;

    /** Radial Distortion co-efficients in NDC (2x1 vector)
	
        D[0] is coefficient k1.  
	D[1] is coefficient k2.
     */
    Vector D;
    
  };
  
  //
  // setting and getting from MetaData
  //
  
  /** get a Camera tag */
  bool get( MetaData &m, const char *path, Camera &c );
  
  /** set a Camera tag */
  void set( MetaData &m, const char *path, const Camera &c );
  
  
} /* namespace */

#endif /* CAMERACALIB_CAMERA_H */

