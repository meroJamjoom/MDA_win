// ==========================================================================
// $Id: Camera.C  bradleyd $
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

#ifndef CAMERACALIB_CAMERA_C
#define CAMERACALIB_CAMERA_C

#include <sstream>

#include "Camera.hh"
#include <MDA/Base/MetaData.hh>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;
  
  
  /** get a Camera tag */
   bool get( MetaData &m, const char *path, Camera &c )
   {
     string p( path );
    
     if (!warnCond(c.K.getNumRows() == 3 && c.K.getNumColumns() == 3 &&
		   c.R.getNumRows() == 3 && c.R.getNumColumns() == 3 &&
		   c.T.getSize() == 3 && c.D.getSize() == 3, "dimensional mismatch in Camera"))
       return false;

     bool success = true;
    
     success&= get( m, (p+".intrinsic.K").c_str(), c.K);
     success&= get( m, (p+".intrinsic.D").c_str(), c.D);
     success&= get( m, (p+".extrinsic.R").c_str(), c.R);
     success&= get( m, (p+".extrinsic.T").c_str(), c.T);

     return success;
   }
  
  /** set a Camera tag */
  void set( MetaData &m, const char *path, const Camera &c )
  {
      string p( path );
      
      // for some reason, the compiler makes me copy these matrices
      // and vectors locally before I can call the set function. 
      // I think it has something to do with consts, but I can't figure it out.
      Matrix K = c.K;
      Matrix R = c.R;
      Vector T = c.T;
      Vector D = c.D;
      set (m, (p+".intrinsic.K").c_str(), K);
      set (m, (p+".intrinsic.D").c_str(), D);
      set (m, (p+".extrinsic.R").c_str(), R);
      set (m, (p+".extrinsic.T").c_str(), T);
  }
  

} /* namespace */

#endif /* CAMERACALIB_CAMERA_C */

