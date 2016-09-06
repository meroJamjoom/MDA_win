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

#ifndef LINEARALGEBRA_MATRIX_C
#define LINEARALGEBRA_MATRIX_C

#include <math.h>

#include "MDA/Threading/SMPJobManager.hh"

#include <MDA/Config.hh>
#include "Matrix.hh"
#include "Vector.hh"
#include "LinAlgThreading.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** copy a matrix (into a new memory object memory) */
void
Matrix::copy( const Matrix &orig, bool fortranLayout )
{
  // de-reference current matrix, create new one
  m->deref();
  m= new Mat( orig.m->rows, orig.m->columns, fortranLayout );
  
  // then use assignment method for actually moving the data
  assign( orig );
}


/** assign the values of another matrix to the memory currently held
    by this matrix (dimensions must match) */
Matrix &
Matrix::assign( const Matrix &orig )
{
  unsigned long i, j;
  
  // check for matching dimensions
  errorCond( m->rows== orig.m->rows && m->columns== orig.m->columns,
	     "  matrix dimensions do not match!" );
  
  // now copy the data from the original
  double *d= m->data + m->offset;
  double *s= orig.m->data + orig.m->offset;
  if( m->columnStride== 1ul && orig.m->columnStride== 1ul &&
      m->rowStride== m->columns && orig.m->rowStride== m->columns )
    // simple, successive memory layout on both sides -> use memcpy
    memcpy( d, s, m->rows * m->columns * sizeof( double ) );
  else
    // else, copy by element
    for( i= 0 ; i< m->rows ; i++ )
    {
      double *dr= d + i*m->rowStride;
      double *sr= s + i*orig.m->rowStride;
      for( j= m->columns ; j> 0 ;
	   j--, dr+= m->columnStride, sr+= orig.m->columnStride )
	(*dr)= (*sr);
    }
  
  return *this;
}


/** Frobenius norm (sum of sqares of all entries) */
double
Matrix::frobeniusNorm()
{
  unsigned long i;
  double result= 0.0;
  if( m->rows< LINALG_SMP_THRESHOLD*LINALG_SMP_THRESHOLD )
  {
    // serial computation for small matrices
    for( i= 0 ; i< m->rows ; i++ )
      result+= dot( m->rowVec[i], m->rowVec[i] );
  }
  else
  {
    // threading for large matrixes
    SMPJobList jobs;
    for( i= 0 ; i< m->rows ; i++ )
      jobs.push_back( new DotProductJob( m->rowVec[i], m->rowVec[i],
					 result, true ) );
    SMPJobManager::getJobManager()->batch( jobs );
  }
  return sqrt( result );
}

/** matrix-vector product */
Vector &
Matrix::rightMultiply( const Vector &v, Vector &result ) const
{
  unsigned long i;
  errorCond( m->rows== result.getSize() && m->columns== v.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  if( m->rows< LINALG_SMP_THRESHOLD )
    // for small matrices, we use a single thread
    for( i= 0 ; i< m->rows ; i++ )
      result[i]= dot( m->rowVec[i], v );
  else
  {
    // for large matrices, we use threading (if enabled)
    SMPJobList dpJobs;
    for( i= 0 ; i< m->rows ; i++ )
      dpJobs.push_back( new DotProductJob( m->rowVec[i], v, result[i] ) );
    SMPJobManager::getJobManager()->batch( dpJobs );
  }
  
  return result;
}

/** vector-matrix product (i.e. left multiply of a row vector) for
    non-sparse matrices */
Vector &
Matrix::leftMultiply( const Vector &v, Vector &result ) const
{
  unsigned long i;
  errorCond( m->columns== result.getSize() && m->rows== v.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  result.zero();
  if( m->columns< LINALG_SMP_THRESHOLD )
    // for small matrices, we use a single thread
    for( i= 0 ; i< m->columns ; i++ )
      result[i]= dot( v, getColumnVector( i ) );
  else
  {
    // for large matrices, we use threading (if enabled)
    SMPJobList dpJobs;
    for( i= 0 ; i< m->columns ; i++ )
      dpJobs.push_back( new DotProductJob( v, getColumnVector( i ),
					   result[i]) );
    SMPJobManager::getJobManager()->batch( dpJobs );
  }
  
  return result;
}


//
// Methods
//

/** matrix-matrix product */
Matrix &
multMatrixMatrix( const Matrix &m1, const Matrix &m2, Matrix &result )
{
  unsigned long i, j;
  
  unsigned long rows1= m1.getNumRows();
  unsigned long cols1= m1.getNumColumns();
  unsigned long rows2= m2.getNumRows();
  unsigned long cols2= m2.getNumColumns();
  errorCond( rows1== result.getNumRows() && cols2== result.getNumColumns() &&
	     cols1== rows2, "  incompatible matrix dimensions" );
  
  if( cols1 * cols2 < LINALG_SMP_THRESHOLD )
  {
    // if little work is to be done for each new row, just do things
    // serially
    for( i= 0 ; i< rows1 ; i++ )
      for( j= 0 ; j< cols2 ; j++ )
	result.set( i, j, dot( m1[i], m2.getColumnVector( j ) ) );
  }
  else
  {
    SMPJobManager *scheduler= SMPJobManager::getJobManager();
    SMPJobList jobs;
    
    for( i= 0 ; i< rows1 ; i++ )
      for( j= 0 ; j< cols2 ; j++ )
	jobs.push_back( new DotProductJob( m1[i], m2.getColumnVector( j ),
					   result[i][j] ) );
    scheduler->batch( jobs );
  }
  
  return result;
}


//
// MetaData set/get
//

/** specialized get for Matrix objects */
bool
get( MetaData &m, char const *path, Matrix &val )
{
  string p( path );
  
  // see whether the vector is inlined or stored in a file by reading
  // the "type" attribute. Default is "inline"
  string type;
  if( !m.get( (p+":type").c_str(), type ) )
    type= "inline";
  
  if( type== "inline" )
  {
    // read an inlined array of doubles
    
    // first, get the number of rows and columns
    unsigned rows, cols;
    if( !warnCond( get( m, (p+":rows").c_str(), rows ),
		   "missing row attribute in inlined Matrix") ||
	!warnCond( get( m, (p+":cols").c_str(), cols ),
		   "missing columns (cols) attribute in inlined Matrix") )
      return false;
    
    // create matrix
    val= Matrix( rows, cols );
    // then read elements as raw array of doubles
    return get( m, path, &val[0][0], rows*cols );
  }
  else if( type== "file" ) 
  {
    // read from a file
    
    // get file name
    string file;
    if( !m.get( path, file ) )
      return false;
    
    // Get the desired channel index
    unsigned channel;
    if( !get( m, (p+":channel").c_str(), channel ) )
      channel= 0u;
    
    // read array, 
    Array<double> array;
    if( !warnCond( array.read( file.c_str() ), "cannot read file" ) )
      return false;
    
    if( !warnCond( channel< array.getNumChannels(), "channel out of range" ) )
      return false;
    
    val= Matrix( array, channel );
    return true;
  }
  
  // if we reach this, things have gone wrong
  warning( "unknown storage type" );
  return false;  
}


/** specialized set for Matrix objects */
void
set( MetaData &m, char const *path, Matrix &val )
{
  string p( path );
  
  // set the row and column attributes and then just use the raw array write
  unsigned rows= val.getNumRows();
  unsigned cols= val.getNumColumns();
  set( m, (p+":rows").c_str(), rows );
  set( m, (p+":cols").c_str(), cols );
  set( m, path, &val[0][0], rows*cols );
}


} /* namespace */

#endif /* LINEARALGEBRA_MATRIX_C */

