// ==========================================================================
// $Id:$
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

#ifndef LAPACK_LAPACKBRIDGE_C
#define LAPACK_LAPACKBRIDGE_C

#include <MDA/Config.hh>
#include "LAPACKBridge.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// declaration of Fortran LAPACK functions
extern "C" {
  // singular value decomposition for larger matrices
  void dgesdd_( char *JOBZ, long *M, long *N, double *A, long *LDA,
		double *S, double *U, long *LDU, double *VT, long *LDVT,
		double *WORK, long *LWORK, long *IWORK, long *INFO );
  
  // singular value decomposition for smaller matrices
  void dgesvd_( char *JOBU, char *JOBVT, long *M, long *N,
		double *A, long *LDA, double *S, double *U, long *LDU,
		double *VT, long *LDVT, double *WORK, long *LWORK,
		long *info );
  
  // LU - factorization
  void dgetrf_( long *R, long *C, double *A, long *LDA,
		long *IPIV, long *INFO );
  
  // matrix inverse, given an LU factorization
  void dgetri_( long *N, double *A, long *LDA,
		long *IPIV, double *WORK, long *LWORK, long *INFO ); 
}

/** inverse of a square matrix -- the original matrix could be
    destroyed unless copyFirst is specified */
Matrix
inverse( Matrix &m, bool preserveMatrix )
{
  long info= 0l;
  long rows= m.getNumRows();
  errorCond( rows== m.getNumColumns(), "  not a square matrix!" );
  
  // use original matrix if we can, create Fortran-layout copy if not
  Matrix mat;
  if( preserveMatrix || !m.hasFortranLayout() )
    mat.copy( m, true );
  else
    mat= m;
  
  // first, do LU factorization
  long *pivot= new long[rows];
  dgetrf_( &rows, &rows, &(mat[0][0]), &rows, pivot, &info );
  errorCond( info>= 0, "  dgetrf - illegal value\n" );
  warnCond( info<= 0, "  singular matrix\n" );
  
  // before the final solve, query the optimal worksize first
  long workSize;
  double ws;
  workSize= -1;
  dgetri_( &rows, &(mat[0][0]), &rows, pivot, &ws, &workSize, &info );
  workSize= (long)ws;
  //!! do something about really large space requirements
  warnMem( workSize*sizeof(double), MEGA_BYTES( 8 ) );
  double *work= new double[workSize];
  
  // then, solve the system using the factorization
  dgetri_( &rows, &(mat[0][0]), &rows, pivot, work, &workSize, &info );
  errorCond( info>= 0, "  dgetri - illegal value\n" );
  warnCond( info<= 0, "  singular matrix\n" );
  
  delete [] work;
  delete [] pivot;
  return mat;
}


/** singular value decomposition
 *  we use LAPACK's dgesvd for small matrixes (m,n <
 *  SVD_SMALL_MATRIX_THRESHOLD), and dgesdd for large ones
 */
Vector
svd( Matrix &m, Matrix &U, Matrix &Vt )
{
  //
  // set up paramters
  //
  char JOBU= 'A';
  char JOBVT= 'A';
  long M= m.getNumRows();
  long N= m.getNumColumns();
  bool large= M> SVD_SMALL_MATRIX_THRESHOLD || N> SVD_SMALL_MATRIX_THRESHOLD;
  
  // re-allocate input matrix Fortran-style if necessary
  Matrix A;
  if( !m.hasFortranLayout() )
    A.copy( m, true );
  else
    A= m;
  
  // allocate output vector and matrices (Fortran style)
  Vector S( M< N ? M : N );
  U= *(new Matrix( M, M, false, true )); // Fortran-style MxM matrix, no init
  Vt= *(new Matrix( N, N, false, true )); // Fortran-style NxN matrix, no init
  
  long INFO= 0l;
  
  // query best work area size & allocate work area
  long LWORK= -1l; // ask for correct size
  long *IWORK= new long[M<N ? 8*M : 8*N];
  double ws;
  if( large )
    dgesdd_( &JOBU, &M, &N, &(A[0][0]), &M, &(S[0]), &(U[0][0]), &M,
	     &(Vt[0][0]), &N, &ws, &LWORK, IWORK, &INFO );
  else
    dgesvd_( &JOBU, &JOBVT, &M, &N, &(A[0][0]), &M, &(S[0]),
	     &(U[0][0]), &M, &(Vt[0][0]), &N, &ws, &LWORK, &INFO );
  LWORK= (long)ws;
  //!! do something about really large space requirements
  warnMem( LWORK*sizeof(double), MEGA_BYTES( 8 ) );
  double *WORK= new double[LWORK];
  
  //
  // actually compute the SVD
  //
  if( large )
    dgesdd_( &JOBU, &M, &N, &(A[0][0]), &M, &(S[0]), &(U[0][0]), &M,
	     &(Vt[0][0]), &N, WORK, &LWORK, IWORK, &INFO );
  else
    dgesvd_( &JOBU, &JOBVT, &M, &N, &(A[0][0]), &M, &(S[0]),
	     &(U[0][0]), &M, &(Vt[0][0]), &N, WORK, &LWORK, &INFO );
  errorCond( INFO>= 0, "  illegal value\n" );
  warnCond( INFO<= 0, "  did not converge\n" );
  
  // clean up
  delete [] WORK;
  delete [] IWORK;
  return S;
}


} /* namespace */

#endif /* LAPACK_LAPACKBRIDGE_C */

