// ==========================================================================
// $Id: PointCorrespondence.hh 744 2010-08-31 23:19:59Z bradleyd $
// point correspondence between 2D (image) and 3D (world) coordinates
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

#ifndef CAMERACALIB_POINTCORRESPONDENCE_H
#define CAMERACALIB_POINTCORRESPONDENCE_H

/*! \file  PointCorrespondence.hh
    \brief point correspondence between 2D (image) and 3D (world) coordinates
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/Base/MetaData.hh>
#include <MDA/LinearAlgebra/LinAlg.hh>

namespace MDA {
  
  /** \class PointCorrespondence PointCorrespondence.hh
      point correspondence between 2D (image) and 3D (world) coordinates */
  
  class PointCorrespondence {

  public:

    /** empty default constructor */
    inline PointCorrespondence()
      : iPt(2),wPt(3)
    {}
    
    /** constructor using actual points */
    inline PointCorrespondence( const Vector &imagePt, const Vector &worldPt )
    {
      iPt.copy( imagePt );
      wPt.copy( worldPt );
    }
    
    /** 2D point in image space */
    Vector iPt;
    
    /** 3D point in world space */
    Vector wPt;
    
  };
  
  //
  // setting and getting from MetaData
  //
  
  /** get a single PointCorrespondence tag */
  bool get( MetaData &m, const char *path, PointCorrespondence &c );
  
  /** set a single PointCorrespondence tag */
  void set( MetaData &m, const char *path, const PointCorrespondence &c );
  
  /** get a vector of PointCorrespondence tags that are all siblings
      under a common <corr> tag */
  bool get( MetaData &m, const char *path, vector<PointCorrespondence> &c );
  
  /** set a vector of PointCorrespondence tags that are all siblings
      under a common <corr> tag */
  void set( MetaData &m, const char *path,
	    const vector<PointCorrespondence> &c );
  
} /* namespace */

#endif /* CAMERACALIB_POINTCORRESPONDENCE_H */

