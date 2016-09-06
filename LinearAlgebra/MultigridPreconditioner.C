// ==========================================================================
// $Id: MultigridPreconditioner.C 782 2011-02-09 11:11:01Z nasarouf $
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

#ifndef IMAGESPACESYSTEMS_MULTIGRIDPRECONDITIONER_C
#define IMAGESPACESYSTEMS_MULTIGRIDPRECONDITIONER_C

#include <float.h>

#include "MDA/Base/Errors.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "MultigridPreconditioner.hh"
#include "ImageSpaceSystem.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  void MultigridPreconditioner::init( const LinearOperator &m)
  {
    rightMultiplyCount=0;
    firstiter=1;
    veryfirst=0;
  } 

  void MultigridPreconditioner::rightMultiply( const Vector &v, Vector &result ){    
    int numDimensions = system.getDimensions().vec.size();
    int systemSize = 1;
    for(int i = 0; i < numDimensions; ++i)
      systemSize *= system.getDimensions().vec[i];
    
    errorCond( initialized,
	       "  preconditioner has to be initialized first" );
    errorCond( systemSize== result.getSize() && systemSize== v.getSize(),
	       "  incompatible matrix/vector dimensions" );

    result.assign(v);
      
    if (system.getHalfSystem()!=NULL){
      unsigned step=rightMultiplyCount/gamma;
      if (step<maxRecursionCount 
	  && (rightMultiplyCount++)%gamma==0 ) {
        int numDimensionsHalfSys = system.getDimensions().vec.size();
        int systemSizeHalfSys = 1;
        for(int i = 0; i < numDimensionsHalfSys; ++i)
          systemSizeHalfSys *= system.getHalfSystem()->getDimensions().vec[i];
	Vector half(systemSizeHalfSys,0.0);
	Vector solution(systemSizeHalfSys,0.0);
	Vector rhs(systemSizeHalfSys,0.0);
	system.restrict(v,half);

	if (system.getHalfSystem()->getHalfSystem()==NULL){
	  halfSolver.setMaxNumIter(-1);
	}else if (halfMultigridPreconditioner->getMaxRecursionCount()>3){  
	  halfSolver.setMaxNumIter(-1);
	  halfMultigridPreconditioner->setMaxRecursionCount(3);
	}else if (halfMultigridPreconditioner->getMaxRecursionCount()==3){
	  halfSolver.setMaxNumIter(2*gamma-1);
	  halfMultigridPreconditioner->setMaxRecursionCount(1);
	}
	halfSolver.solve( *system.getHalfSystem(), half, solution);
	int it=halfSolver.getNumIter();
        system.prolongate(solution,result);
      }
    }
  }
} /* namespace */

#endif /* MULTIGRIDPRECONDITIONER_C */

