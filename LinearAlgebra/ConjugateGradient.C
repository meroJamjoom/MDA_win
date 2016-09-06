// ==========================================================================
// $Id: ConjugateGradient.C 831 2011-08-18 23:21:57Z joelaf $
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

#ifndef LINEARALGEBRA_CONJUGATEGRADIENT_C
#define LINEARALGEBRA_CONJUGATEGRADIENT_C

#include "MDA/Base/Errors.hh"

#include "Preconditioner.hh"
#include "JacobiPreconditioner.hh"
#include "ConjugateGradient.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/* solving the linear system Ax = b */
void
ConjugateGradient::solve( const LinearOperator &A, const Vector &b, Vector &x )
{ 
  unsigned size= A.getNumRows();
  errorCond( size== A.getNumColumns() && size== b.getSize() &&
	     size== x.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  // initialize preconditioner
  precond->init( A );
  
  Vector h( size );
  double a, g, zDotR, newZDotNewR;
  
  iter= 0;
  
  Vector r( size );
  A.rightMultiply( x, r );
  (r-= b).negate();			// r= b - Ax
  

  Vector z( size );
  precond->rightMultiply( r, z );	// z= Preconditioner * r
  Vector p( size );
  p.copy( z );				// p= z
  
  zDotR= dot( z, r );  
  
  while( iter<maxiter && zDotR> threshold )
  { 
    A.rightMultiply( p, h );		// h= Ap
    a= zDotR / dot( p, h );		// a= dot(z,r)/dot(p,h)
    x.addScalarTimesVector( a, p );	// x= x + ap
    r.addScalarTimesVector( -a, h );	// r= r - ah
    precond->rightMultiply( r, z );	// z= Preconditioner * r
    newZDotNewR= dot( z, r );
    g= newZDotNewR / zDotR;		// g= dot(new_z,new_r)/dot(old_z,old_r)
    p*= g;
    p+= z;				// p= z + g*p
    zDotR= newZDotNewR;
    iter++;
  }
}
  
/** solving the linear system A^T*Ax = A^T*b (without explicitly
    computing the matrix) */
void
ConjugateGradient::solveLeastSquares( const LinearOperator &A,
				      const Vector &b, Vector &x )
{
  unsigned size= A.getNumColumns();
  unsigned hhSize= A.getNumRows();
  errorCond( hhSize== b.getSize() && size== x.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  // initialize preconditioner
  precond->initLeastSquares( A );
  
  // right hand side: ab= A^T * b
  Vector ab( size );
  A.leftMultiply( b, ab );
  
  Vector h( size );
  Vector hh( hhSize ); 
  double a, g, zDotR, newZDotNewR;
  
  iter= 0;
  
  Vector r( size );
  A.rightMultiply( x, hh );		// r= A^T * A * x
  A.leftMultiply( hh, r );
  (r-= ab).negate();			// r= A^T * b - A^T * Ax
  Vector z( size );
  precond->rightMultiply( r, z );	// z= Preconditioner * r
  Vector p( size );
  p.copy( z );				// p= z
  
  zDotR= dot( z, r );  
  do
  {
    A.rightMultiply( p, hh );		// h= A^T & A * p
    A.leftMultiply( hh, h );
    a= zDotR / dot( p, h );		// a= dot(z,r)/dot(p,h)
    x.addScalarTimesVector( a, p );	// x= x + ap
    r.addScalarTimesVector( -a, h );	// r= r - ah
    precond->rightMultiply( r, z );	// z= Preconditioner * r
    newZDotNewR= dot( z, r );
    g= newZDotNewR / zDotR;		// g= dot(new_z,new_r)/dot(old_z,old_r)
    p*= g;
    p+= z;				// p= r + g*p
    zDotR= newZDotNewR;
    
    iter++;
  }
  while( zDotR> threshold );
}

} /* namespace */

#endif /* LINEARALGEBRA_CONJUGATEGRADIENT_C */

