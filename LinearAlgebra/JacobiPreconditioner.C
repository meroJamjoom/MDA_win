// ==========================================================================
// $Id: JacobiPreconditioner.C 326 2009-06-04 01:25:51Z curtd $
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

#ifndef LINEARALGEBRA_JACOBIPRECONDITIONER_C
#define LINEARALGEBRA_JACOBIPRECONDITIONER_C

#include <float.h>

#include "MDA/Base/Errors.hh"
#include "Vector.hh"
#include "Matrix.hh"
#include "SparseMatrix.hh"
#include "JacobiPreconditioner.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;
  
/** initialize the preconditioner for a specific linear system */
void
JacobiPreconditioner::init( const LinearOperator &m )
{
  unsigned long i;
  numElem= m.getNumColumns();
  errorCond( numElem== m.getNumRows(), "  require square matrix..." );
  
  // obtain the diagonal of the matrix
  const Matrix *mm= dynamic_cast<const Matrix *>( &m );
  const SparseMatrix *ms= dynamic_cast<const SparseMatrix *>( &m );
  
  if( mm!= NULL )
    // m is a Matrix
    diag= mm->getDiagonal();
  else if( ms!= NULL )
  {
    // m is a SparseMatrix
    diag= Vector( numElem );
    for( i= 0 ; i< numElem ; i++ )
      diag[i]= ms->get( i, i );
  }
  else
  {
    // m is of an unsupported type
    warning( "  unsupported matrix type - using identity preconditioner\n" );
    diag= Vector( numElem, 1.0 );
  }
  
  // invert the vector
  for( i= 0 ; i< numElem ; i++ )
  {
    if( fabs( diag[i] )> DBL_EPSILON )
      diag[i]= 1.0 / diag[i];
    else
      diag[i]= 1.0;
  }
  initialized= true;
}


/** initialize the preconditioner for a specific least squares
    system m^T * m */
void
JacobiPreconditioner::initLeastSquares( const LinearOperator &m )
{
  unsigned long i;
  numElem= m.getNumColumns();
  diag= Vector( numElem );
  
  // obtain the diagonal of m^T * m
  const Matrix *mm= dynamic_cast<const Matrix *>( &m );
  if( mm!= NULL )
  {
    // m is of type Matrix
    for( i= 0 ; i< numElem ; i++ )
    {
      Vector col= mm->getColumnVector( i );
      double d= dot( col, col );
      if( fabs( d )> DBL_EPSILON )
	diag[i]= 1.0 / d;
      else
	diag[i]= 1.0;
    }
  }
  else
  {
    // m is not of a supported matrix type
    warning( "  unsupported matrix type - using identity preconditioner\n" );
    diag= Vector( numElem );
    for( i= 0 ; i< numElem ; i++ )
      diag[i]= 1.0;
  }
  
  initialized= true;
}
    

    
/** right multiply the preconditioner matrix P^{-1} with a vector */
void
JacobiPreconditioner::rightMultiply( const Vector &v, Vector &result )
{
  errorCond( initialized,
	     "  preconditioner hasn to been initialized to a specific matrix" );
  errorCond( numElem== result.getSize() && numElem== v.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  // multipy vector with diagonal matrix
  for( unsigned long i= 0 ; i< numElem ; i++ )
    result[i]= v[i]*diag[i];
}


} /* namespace */

#endif /* LINEARALGEBRA_JACOBIPRECONDITIONER_C */

