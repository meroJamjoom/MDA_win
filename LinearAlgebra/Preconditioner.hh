// ==========================================================================
// $Id: Preconditioner.hh 292 2009-03-20 03:08:52Z heidrich $
// a baseclass for preconditioners of linear solvers (this class produces an identity matrix "preconditioner", that is a null operation. See subclasses for actual preconditioners)
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

#ifndef LINEARALGEBRA_PRECONDITIONER_H
#define LINEARALGEBRA_PRECONDITIONER_H

/*! \file  Preconditioner.hh
    \brief a baseclass for preconditioners of linear solvers. This
    class produces an identity matrix "preconditioner". See subclasses
    for actual preconditioners.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "LinearOperator.hh"

namespace MDA {

  /** \class Preconditioner Preconditioner.hh
      A baseclass for preconditioners of linear solvers. This class
      produces an identity matrix "preconditioner". See subclasses for
      actual preconditioners. */
  class Preconditioner {

  public:

    /** initialize the preconditioner for a specific linear system */
    inline virtual void init( const LinearOperator &m )
    {}
    
    /** initialize the preconditioner for a specific least squares system */
    inline virtual void initLeastSquares( const LinearOperator &m )
    {}
    
    /** right multiply the preconditioner matrix P^{-1} with a vector */
    virtual void rightMultiply( const Vector &v, Vector &result );
    
    /** left multiply the preconditioner matrix P^{-1} with a vector
     *  (the implementation here applies to all symmetric preconditioners) */
    inline virtual void leftMultiply( const Vector &v, Vector &result )
    {
      rightMultiply( v, result );
    }
    
  };

} /* namespace */

#endif /* LINEARALGEBRA_PRECONDITIONER_H */

