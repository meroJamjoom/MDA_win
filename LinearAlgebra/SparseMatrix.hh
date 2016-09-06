// ==========================================================================
// $Id: SparseMatrix.hh 753 2010-09-09 00:54:03Z bradleyd $
// a sparse matrix representation
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

#ifndef LINEARALGEBRA_SPARSEMATRIX_H
#define LINEARALGEBRA_SPARSEMATRIX_H

/*! \file  SparseMatrix.hh
    \brief a sparse matrix representation
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>

#include "LinearOperator.hh"

namespace MDA {
  
  using namespace std;
  
  /** \class SparseMatrix SparseMatrix.hh
      a sparse matrix representation */
  
  class SparseMatrix: public LinearOperator {

  public:

    /** default constructor */
    SparseMatrix();
    
    /** constructor with given dimensions */
    SparseMatrix( unsigned long rows, unsigned long columns )
    {
      m= new Mat( rows, columns );
    }
    
    /** destructor */
    ~SparseMatrix()
    {
      m->deref();
    }

    //
    // implementations of LinearOperator methods
    //
    
    /** get number of rows */
    inline virtual unsigned long getNumRows() const
    {
      return m->rows;
    }
    
    /** get number of columns */
    inline virtual unsigned long getNumColumns() const
    {
      return m->columns;
    }
    
    /** right-multiply a vector */
    virtual Vector &rightMultiply( const Vector &v, Vector &result ) const;
    
    /** left-multiply a vector */
    virtual Vector &leftMultiply( const Vector &v, Vector &result ) const;
    
    //
    // initialization
    //
    
    /** set matrix to zero */
    inline void zero()
    {
      for( unsigned long i= 0 ; i< m->rows ; i++ )
      {
	m->data[i].data.clear();
	m->data[i].numEntries= m->data[i].firstEntry= m->data[i].lastEntry= 0ul;
      }
    }
    
    /** set matrix to identity */
    inline void identity()
    {
      zero();
      
      unsigned long numElem= m->rows< m->columns ? m->rows : m->columns;
      for( unsigned long i= 0 ; i< numElem ; i++ )
	m->data[i].set( i, 1.0 );
    }
    
    // access operators
    
    /** set an element to a new value */
    inline void set( unsigned long r, unsigned long c, double val )
    {
      errorCond( r< m->rows && c< m->columns, "  index out of range!" );
      
      m->data[r].set( c, val );
    }
    
    /** get the current value of an element */
    inline double get( unsigned long r, unsigned long c ) const
    {
      errorCond( r< m->rows && c< m->columns, "  index out of range!" );
      
      return m->data[r].get( c );
    }
    
    /** helper function for matrix-vector product - dot product of a
	row and a vector */
    double rowDot( unsigned long row, const Vector &v ) const;
    
    /** helper function for vector-matrix product - add a multiple of
	a row to a vector */
    void addRowMult( unsigned long row, double mult, Vector &accum ) const;
    
  protected:
    
    /** \class Entry SparseMatrix.hh
	an entry with the associated (column) index */
    class Entry {
    public:
      /** constructor */
      inline Entry( unsigned long _index, double _value )
	: index( _index ), val( _value )
      {}
      /** element index */
      unsigned long index;
      /** value */
      double val;
    };
    
    /** \class Vec SparseMatrix.hh
	a (sparse) row vector */
    class Vec {
    public:
      /** constructor */
      inline Vec()
	: numEntries( 0ul ), firstEntry( 0ul ), lastEntry( 0ul )
      {}
      
      /** set an element of the vector */
      void set( unsigned long index, double value );
      
      /** get value of an element */
      double get( unsigned long index ) const;
      
      /** index of the first non-zero entry */
      unsigned long firstEntry;
      /** index of the last non-zero entry */
      unsigned long lastEntry;
      /** number of non-zero entries */
      unsigned numEntries;
      /** (unsorted) array of entries */
      vector<Entry> data;
    };
    
    /** class Mat SparseMatrix.hh
	a memory object storing sparse matrix data */
    class Mat: public MemoryObject {
    public:
      /** constructor */
      inline Mat( unsigned long _rows= 0ul, unsigned long _cols= 0ul )
	: MemoryObject(), rows( _rows ), columns( _cols )
      {
	data.resize( rows );
      }
      
      /** sparse matrix data, by row (hides the superclass data, which
	  is of type double, and not used in sparse matrices) */
      vector<Vec> data;
      
      /** whether or not the matrix has been transposed */
      bool transposed;
      
      /** number of rows */
      unsigned long rows;
      
      /** number of rows */
      unsigned long columns;
      
    protected:
      
      /** destructor - private (use deref to delete references) */
      virtual ~Mat()
      {}
      
    };
    
    
    /** reference to the actual sparse matrix data */
    Mat *m;
    
  };
  
} /* namespace */

#endif /* LINEARALGEBRA_SPARSEMATRIX_H */

