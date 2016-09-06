// ==========================================================================
// $Id: JacobiRotation.C 447 2009-10-23 06:58:37Z heidrich $
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

#ifndef LINEARALGEBRA_JACOBIROTATION_C
#define LINEARALGEBRA_JACOBIROTATION_C

#include "MDA/Base/Errors.hh"

#include "JacobiRotation.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** compute vectors e_i, scalars d_i with A*e_i = d_i * e_i
    (Eigenvectors e_i are the ROWS of matrix E. They are sorted by
    Eigenvalue.) */
void
JacobiRotation::solve( Matrix &A, Vector &d, Matrix &E )
{
  unsigned size= d.getSize();
  errorCond( size== A.getNumRows() && size== A.getNumColumns() &&
	     size== E.getNumRows() && size== E.getNumColumns(),
	     "  incompatible matrix/vector dimensions" );
  
  unsigned i, j, k, l;
  Vector b( size );
  Vector z( size, 0.0 );
  double c, s, g, h, t, theta;
  
  E.identity();
  // b, d are initially the diagonal of A
  for( i= 0 ; i< size ; i++ )
    b[i]= d[i]= A[i][i];
  
  // iterate
  iter= 0;
  for( k= 0 ; k< maxIter ; k++ )
  {
    // sum of off-axis elements
    double sum= 0.0;
    for( i= 0 ; i< size ; i++ )
      for( j= i+1 ; j< size ; j++ )
	sum+= fabs( A[i][j] );
    
    // terminate if the upper triangle is close to zero
    if( sum<= threshold )
      break;
    
    // first 4 iterations vs the rest
    double thres= (k< 4) ? 0.2*sum / (size*size) : 0.0;
    
    // traverse upper triangle
    for( i= 0 ; i< size ; i++ )
      for( j= i+1 ; j< size ; j++ )
      {
	g= fabs( A[i][j] );
	// skip rotations for small entries in later iterations
	if( k> 4 &&
	    fabs( d[i] )*zeroThreshold > g &&
	    fabs( d[j] )*zeroThreshold > g )
	  A[i][j]= 0.0;
	else if( g> thres )
	{
	  h= d[j]-d[i];
	  if( fabs( h )*zeroThreshold> g )
	    t= A[i][j] / h;
	    else
	  {
	    theta= 0.5 * h / A[i][j];
	    t= (theta< 0.0 ? -1.0 : 1.0) /
	      (fabs( theta ) + sqrt( 1.0 + theta*theta ));
	  }
	  c= 1.0 / sqrt( 1.0 + t*t );
	  s= t*c;
	  h= t * A[i][j];
	  z[i]-= h;
	  z[j]+= h;
	  d[i]-= h;
	  d[j]+= h;
	  A[i][j]= 0.0;
	  for( l= 0 ; l< i ; l++ )
	  {
	    g= A[l][i]; h= A[l][j];
	    A[l][i]= c*g - s*h;
	    A[l][j]= s*g + c*h;
	  }
	  for( l= i+1 ; l< j ; l++ )
	  {
	    g= A[i][l]; h= A[l][j];
	    A[i][l]= c*g - s*h;
	    A[l][j]= s*g + c*h;
	  }
	  for( l= j+1 ; l< size ; l++ )
	  {
	    g= A[i][l]; h= A[j][l];
	    A[i][l]= c*g - s*h;
	    A[j][l]= s*g + c*h;
	  }
	  for( l= 0 ; l< size ; l++ )
	  {
	    g= E[i][l]; h= E[j][l];
	    E[i][l]=  c*g - s*h;
	    E[j][l]=  s*g + c*h;
	  }
	}
      }
    
    b+= z;
    d.copy( b );
    z.zero();
    iter++;
  }
  
  sort( d );
}


} /* namespace */

#endif /* LINEARALGEBRA_JACOBIROTATION_C */

