// ==========================================================================
// $Id: SparseMatrix.C 319 2009-05-27 21:17:01Z heidrich $
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

#ifndef LINEARALGEBRA_SPARSEMATRIX_C
#define LINEARALGEBRA_SPARSEMATRIX_C

#include "SparseMatrix.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** set an element of the vector */
void
SparseMatrix::Vec::set( unsigned long index, double value )
{
  unsigned i;
  
  if( numEntries== 0 )
  {
    data.push_back( Entry( index, value ) );
    firstEntry= lastEntry= index;
    numEntries= 1;
  }
  else if( index> lastEntry )
  {
    data.push_back( Entry( index, value ) );
    lastEntry= index;
    numEntries++;
  }
  else if( index< firstEntry )
  {
    data.push_back( Entry( index, value ) );
    firstEntry= index;
    numEntries++;
  }
  else
  {
    Entry *d= &(data[0]);
    for( i= 0 ; i< numEntries ; i++ )
      if( d[i].index== index )
      {
	d[i].val= value;
	break;
      }
    if( i== numEntries )
    {
      data.push_back( Entry( index, value ) );
      numEntries++;
    }
  }
}

/** get value of an element */
double
SparseMatrix::Vec::get( unsigned long index ) const
{
  const Entry *d= &(data[0]);
  
  if( index>= firstEntry || index<= lastEntry )
    for( unsigned i= 0 ; i< numEntries ; i++ )
      if( d[i].index== index )
	return d[i].val;
  
  return 0.0;
}

/** helper function for matrix-vector product */
double
SparseMatrix::rowDot( unsigned long row, const Vector &v ) const
{
  const Entry *d= &(m->data[row].data[0]);
  
  double result= 0.0;
  for( unsigned i= m->data[row].numEntries ; i> 0 ; i--, d++ )
    result+= d->val * v[d->index];
  
  return result;
}

/** helper function for vector-matrix product */
void
SparseMatrix::addRowMult( unsigned long row, double mult, Vector &accum ) const
{
  const Entry *d= &(m->data[row].data[0]);
  
  for( unsigned i= m->data[row].numEntries ; i> 0 ; i--, d++ )
    accum[d->index]+= mult * d->val;
}

/** sparse matrix-vector product */
Vector &
SparseMatrix::rightMultiply( const Vector &v, Vector &result ) const
{
  unsigned long i;
  errorCond( m->rows== result.getSize() && m->columns== v.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  //!! multithreading!
  for( i= 0 ; i< m->rows ; i++ )
    result[i]= rowDot( i, v );
  
  return result;
}

/** vector-sparse matrix product */
Vector &
SparseMatrix::leftMultiply( const Vector &v, Vector &result ) const
{
  unsigned long i;
  errorCond( m->rows== v.getSize() && m->columns== result.getSize(),
	     "  incompatible matrix/vector dimensions" );
  
  //!! multithreading
  for( i= 0 ; i< m->rows ; i++ )
    addRowMult( i, v[i], result );
  
  return result;
}
  



} /* namespace */

#endif /* LINEARALGEBRA_SPARSEMATRIX_C */

