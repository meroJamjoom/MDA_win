// ==========================================================================
// $Id: LinearSolver.hh 292 2009-03-20 03:08:52Z heidrich $
// abstract base class for linear solvers
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

#ifndef LINEARALGEBRA_LINEARSOLVER_H
#define LINEARALGEBRA_LINEARSOLVER_H

/*! \file  LinearSolver.hh
    \brief abstract base class for linear solvers
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "LinearOperator.hh"

namespace MDA {

  /** \class LinearSolver LinearSolver.hh
      abstract base class for linear solvers */
  
  class LinearSolver {

  public:

    /** solving the linear system Ax = b */
    virtual void solve( const LinearOperator &A,
			const Vector &b, Vector &x ) = 0;
    
    /** solving the linear system A^T*Ax = A^T*b */
    virtual void solveLeastSquares( const LinearOperator &A,
				    const Vector &b, Vector &x ) = 0;
    
  };

} /* namespace */

#endif /* LINEARALGEBRA_LINEARSOLVER_H */

