// ==========================================================================
// $Id: MultigridPreconditioner.hh 782 2011-02-09 11:11:01Z nasarouf $
// a Multigrid preconditioner for linear solvers
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: nasarouf (Mushfiqur Rouf)
// Email:   nasarouf@cs.ubc.ca
// ==========================================================================

#ifndef MULTIGRIDPRECONDITIONER_H
#define MULTIGRIDPRECONDITIONER_H

/*! \file  MultigridPreconditioner.hh
    \brief a Multigrid preconditioner for linear solvers
 */

#ifndef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/LinearAlgebra/Preconditioner.hh"
#include "MDA/LinearAlgebra/PoissonSystem.hh"
#include "MDA/LinearAlgebra/ConjugateGradient.hh"

namespace MDA {

  /** \class MultigridPreconditioner MultigridPreconditioner.hh
      a Multigrid preconditioner for linear solvers */
  class MultigridPreconditioner: public Preconditioner {

    
  protected:

    int firstiter;

    int veryfirst;

    /** counter used to keep track of how many times it has been
	called. useful for skipping going deeper (check gamma and
	maxRecursionCount) */
    unsigned rightMultiplyCount;

    /** once in every gamma times, go deeper once */
    unsigned gamma;

    /** maximum number of recursions. default=(unsigned)(-1) */
    unsigned maxRecursionCount;
    
    /** whether or not the preconditioner is initialized */
    bool initialized;
    
    /** the system */
    MultigridSystem &system;

    /** the half-system preconditioner */
    MultigridPreconditioner* halfMultigridPreconditioner;

    /** the half solver */
    ConjugateGradient halfSolver;

    
  public:

    void setMaxRecursionCount(unsigned _count){
      maxRecursionCount=_count;
    }

    unsigned getMaxRecursionCount(){
      return maxRecursionCount;
    }
    
    /** constructor
     *  \param _system: the image space system
     *  \param _gamma: once in every gamma times, go deeper. default:2
     */
    inline MultigridPreconditioner( MultigridSystem& _system, 
				    unsigned _gamma=2) 
      : initialized( true ), rightMultiplyCount(0), gamma(_gamma),
	firstiter(0),veryfirst(1),system(_system),maxRecursionCount(-1)
    {
      if (system.getHalfSystem()!=NULL && 
	  system.getHalfSystem()->getHalfSystem()!=NULL){
	halfMultigridPreconditioner=
          new MultigridPreconditioner(*system.getHalfSystem());
	halfSolver.setPreconditioner(halfMultigridPreconditioner);
      }
    }
    
    /** initialize the preconditioner for a specific linear system. 
	We ignore the parameter here. m is ignored here. */
    inline virtual void init( const LinearOperator &m);
    
    /** initialize the preconditioner for a specific least squares
	system m^T * m. m is ignored here. */
    inline virtual void initLeastSquares( const LinearOperator &m ){ init(m); }
    
    /** right multiply the preconditioner matrix P^{-1} with a vector */
    virtual void rightMultiply( const Vector &v, Vector &result );
        
  };

} /* namespace */


#endif /* MULTIGRIDPRECONDITIONER_H */

