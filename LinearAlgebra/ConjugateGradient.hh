// ==========================================================================
// $Id: ConjugateGradient.hh 782 2011-02-09 11:11:01Z nasarouf $
// Conjugate Gradient solver for symmetric, positive-definit matrices.
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

#ifndef LINEARALGEBRA_CONJUGATEGRADIENT_H
#define LINEARALGEBRA_CONJUGATEGRADIENT_H

/*! \file  ConjugateGradient.hh
    \brief Conjugate Gradient solver for symmetric, positive-definite matrices.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/Config.hh>

#include "Preconditioner.hh"
#include "LinearSolver.hh"

namespace MDA {

  /** \class ConjugateGradient ConjugateGradient.hh
      Conjugate Gradient solver for symmetric, positive-definite matrices */
  class ConjugateGradient: public LinearSolver {
    
  public:
    
    /** constructor */
    inline ConjugateGradient()
      : threshold( NUM_ZERO_THRESHOLD ), iter( 0 ),
	precond( new Preconditioner ), maxiter(-1)
    {}
    
    /** constructor */
    inline ConjugateGradient( Preconditioner *pre )
      : threshold( NUM_ZERO_THRESHOLD ), iter( 0 ), 
	precond( pre ), maxiter(-1)
    {}
    
    /** destructor */
    inline ~ConjugateGradient()
    {
      delete precond;
    }
    
    /** solving the linear system Ax = b */
    virtual void solve( const LinearOperator &A,
			const Vector &b, Vector &x );
    
    /** solving the linear system A^T*Ax = A^T*b (without explicitly
	computing the matrix) */
    virtual void solveLeastSquares( const LinearOperator &A,
				    const Vector &b, Vector &x ); 
    
    /** set a new preconditioner */
    inline void setPreconditioner( Preconditioner *pre )
    {
      delete precond;
      precond= pre;
    }
    
    /** set convergence threshold */
    inline void setThreshold( double th )
    {
      threshold= th;
    }
    
    /** get convergence threshold */
    inline double getThreshold( )
    {
      return threshold;
    }
    
    /** get number of iterations (for last run) */
    unsigned getNumIter() const
    {
      return iter;
    }

    /** set max number of iterations (for next run) */
    inline void setMaxNumIter(unsigned iter)
    {
      maxiter=iter;
    }
    
    /** set max number of iterations) */
    inline unsigned getMaxNumIter()
    {
      return maxiter;
    }
    
  protected:

    /** max number of iterations allowed, defaults to inf */
    unsigned maxiter;

    /** preconditioner */
    Preconditioner *precond;
    
    /** termination threshold */
    double threshold;
    
    /** number of iterations from last solve */
    unsigned iter;
    
  };

} /* namespace */

#endif /* LINEARALGEBRA_CONJUGATEGRADIENT_H */

