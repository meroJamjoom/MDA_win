// ==========================================================================
// $Id: JacobiRotation.hh 638 2010-03-08 19:34:24Z heidrich $
// Jacobi Rotation for solving the Eigenvalue problem for a symmetric matrix
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

#ifndef LINEARALGEBRA_JACOBIROTATION_H
#define LINEARALGEBRA_JACOBIROTATION_H

/*! \file  JacobiRotation.hh
    \brief Jacobi Rotation for solving the Eigenvalue problem for a symmetric matrix
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/Config.hh>
#include "EigenSolver.hh"

namespace MDA {
  
  /** \class JacobiRotation JacobiRotation.hh
      Jacobi Rotation for solving the Eigenvalue problem for a
      symmetric matrix */
  
  class JacobiRotation: public EigenSolver {

  public:

    /** default constructor */
    JacobiRotation()
      : threshold( NUM_ZERO_THRESHOLD*NUM_ZERO_THRESHOLD ),
	zeroThreshold( NUM_ZERO_THRESHOLD ),
	maxIter( 50 ), iter( 0 )
    {}
    
    /** compute vectors e_i, scalars d_i with A*e_i = d_i * e_i
	(Eigenvectors e_i are the ROWS of matrix E. They are sorted by
	Eigenvalue.) */
    virtual void solve( Matrix &A, Vector &d, Matrix &E );
    
    /** set convergence threshold */
    inline void setThreshold( double th )
    {
      threshold= th;
    }
    
    /** set maximum number of iterations */
    inline void setMaxIter( unsigned it )
    {
      maxIter= it;
    }
    
    /** report number of iterations in last solve */
    unsigned getNumIter() const
    {
      return iter;
    }
    
  protected:
    
    /** maximum number of iterations */
    unsigned maxIter;
    
    /** convergence threshold */
    double threshold;
    
    /** threshold for setting off-axis elements to zero */
    double zeroThreshold;
    
    /** number of iterations in last solve */
    unsigned iter;
    
  };


} /* namespace */

#endif /* LINEARALGEBRA_JACOBIROTATION_H */

