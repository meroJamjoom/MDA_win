// ==========================================================================
// $Id: JacobiPreconditioner.hh 292 2009-03-20 03:08:52Z heidrich $
// a Jacobi preconditioner for linear solvers
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

#ifndef LINEARALGEBRA_JACOBIPRECONDITIONER_H
#define LINEARALGEBRA_JACOBIPRECONDITIONER_H

/*! \file  JacobiPreconditioner.hh
    \brief a Jacobi preconditioner for linear solvers
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Preconditioner.hh"

namespace MDA {

  /** \class JacobiPreconditioner JacobiPreconditioner.hh
      a Jacobi preconditioner for linear solvers */
  class JacobiPreconditioner: public Preconditioner {
    
  public:
    
    /** default constructor */
    inline JacobiPreconditioner()
      : initialized( false ), diag()
    {}
    
    /** constructor
     *  \param m: the linear system
     *  \param leastSquares: whether the actual matrix is of the form m^T * m
     */
    inline JacobiPreconditioner( const LinearOperator &m,
				 bool leastSquares= false )
    {
      if( leastSquares )
	initLeastSquares( m );
      else
	init( m );
    }
    
    /** initialize the preconditioner for a specific linear system */
    virtual void init( const LinearOperator &m );
    
    /** initialize the preconditioner for a specific least squares
	system m^T * m */
    virtual void initLeastSquares( const LinearOperator &m );
    
    /** right multiply the preconditioner matrix P^{-1} with a vector */
    virtual void rightMultiply( const Vector &v, Vector &result );
    
  protected:
    
    /** whether or not the preconditioner is initialized */
    bool initialized;
    
    /** diagonal of the inverse preconditioning matrix P^{-1} */
    Vector diag;
    
    /** size of the system */
    unsigned numElem;
    
  };

} /* namespace */


#endif /* LINEARALGEBRA_JACOBIPRECONDITIONER_H */

