// ==========================================================================
// $Id: EigenSolver.hh 447 2009-10-23 06:58:37Z heidrich $
// pure baseclass for algorithms that solve an eigenvalue problem
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

#ifndef LINEARALGEBRA_EIGENSOLVER_H
#define LINEARALGEBRA_EIGENSOLVER_H

/*! \file  EigenSolver.hh
    \brief pure baseclass for algorithms that solve an eigenvalue problem
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Vector.hh"
#include "Matrix.hh"

namespace MDA {

  /** \class EigenSolver EigenSolver.hh
      pure baseclass for algorithms that solve an eigenvalue problem */
  
  class EigenSolver {

  public:

    /** compute vectors e_i, scalars d_i with A*e_i = d_i * e_i
	(Eigenvectors e_i are the ROWS of matrix E) */
    virtual void solve( Matrix &A, Vector &d, Matrix &E )= 0;
    
    
    /** \class EigenValue EigenSolver.hh
	a pair of eigen value and index of eigen vector */
    class EigenValue {
    public:
      /** actual eigen value */
      double value;
      /** index for eigen vector, corresponding to a row in the result
	  matrix */
      unsigned long index;
    };
    
    /** an array of EigenValue */ 
    typedef vector<EigenValue> EigenValues;
    
    /** return EigenValues sorted by magnitude */
    inline const EigenValues &getSortedEigenValues() const
    {
      return sortedEigenValues;
    }
    
  protected:
    
    /** sorted eigen values */
    EigenValues sortedEigenValues;
    
    /** create the sorted list of eigen values */
    void sort( Vector &evs );
    
  };


} /* namespace */

#endif /* LINEARALGEBRA_EIGENSOLVER_H */

