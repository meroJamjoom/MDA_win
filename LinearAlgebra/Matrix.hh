// ==========================================================================
// $Id:$
// a n x m matrix
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

#ifndef LINEARALGEBRA_MATRIX_H
#define LINEARALGEBRA_MATRIX_H

/*! \file  Matrix.hh
    \brief a n x m matrix
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>

#include "MDA/Base/Errors.hh"
#include "MDA/Base/MetaData.hh"
#include "MDA/Array/Array.hh"
#include "LinearOperator.hh"

namespace MDA {

  using namespace std;
  
  /** \class Matrix Matrix.hh
      a n x m matrix */
  class Matrix: public LinearOperator {

  public:
    
    /** default constructor - yields uninitialized matrix */
    inline Matrix()
    {
      m= new Mat();
    }
    
    /** constructor for matrix (zero matrix if initialized) */
    inline Matrix( unsigned long rows, unsigned long cols,
		   bool initialize= false, bool fortranLayout= false )
    {
      m= new Mat( rows, cols, fortranLayout );
      if( initialize )
	zero();
    }
    
    /** constructor for an outer product matrix */
    inline Matrix( const Vector &v1, const Vector &v2 )
    {
      unsigned long i, j;
      unsigned long size1= v1.getSize();
      unsigned long size2= v2.getSize();
      
      m= new Mat( size1, size2 );
      
      double *d= m->data;
      for( i= 0 ; i< size1 ; i++ )
	for( j= 0 ; j< size2 ; j++, d++ )
	  (*d)= v1[i] * v2[j];
    }
    
    /** contructor from a channel in an array */
    inline Matrix( Array<double> &a, unsigned ch )
    {
      CoordinateVector dim= a.getDimension();
      errorCond( dim.vec.size()== 2, "  only 2D arrays supported" );
      
      m= new Mat( dim.vec[1], dim.vec[0], a[ch], 0, dim.vec[0], 1 );
    }
    
    /** contructor from external raw data (this data is NOT copied,
	and changes to the vector will overwrite it!) */
    inline Matrix( unsigned long rows, unsigned long cols,
		   double *rawData, unsigned long offset= 0ul,
		   unsigned long rStride= 0ul, unsigned long cStride= 1 )
    {
      m= new Mat( rows, cols, rawData, offset, rStride, cStride );
    }
    
    /** constructor for matrix sharing data with another memory object */
    inline Matrix( unsigned long rows, unsigned long cols,
		   MemoryObject *owner, unsigned long offset= 0ul,
		   unsigned long rStride= 0ul, unsigned long cStride= 1 )
    {
      m= new Mat( rows, cols, owner, offset, rStride, cStride );
    }
    
    /** copy constructor */
    inline Matrix( const Matrix &mat )
    {
      m= (Mat *)(mat.m->ref());
    }
    
    /** assignment operator (adds a reference to an existing matrix -
	for a physical copy use the copy() or assign() methods) */
    inline const Matrix &operator=( const Matrix &mat )
    {
      // de-reference current vector
      m->deref();
      // and reference the target instead
      m= (Mat *)(mat.m->ref());
      return *this;
    }
    
    /** copy a matrix into a newly created memory object, which
	becomes the property of this matrix */
    void copy( const Matrix &orig, bool fortranLayout= false );
    
    /** assign the values of another matrix to the memory currently
	held by this matrix (dimensions must match) */
    Matrix &assign( const Matrix &orig );
    
    /** destructor */
    inline ~Matrix()
    {
      m->deref();
    }
    
    /** check if the matrix is stored in Fortran layout (i.e. if it is
	ready to be used by LINPACK routines) */
    inline bool hasFortranLayout() const
    {
      return (m->rowStride== 1) & (m->columnStride== m->rows);
    }
    
    /** create a submatrix that shares memory with this matrix */
    inline Matrix getSubMatrix( unsigned long firstRow, unsigned long lastRow,
				unsigned long firstCol, unsigned long lastCol,
				unsigned long rowStride= 1ul,
				unsigned long colStride= 1ul )
    {
      return Matrix( (lastRow-firstRow)/rowStride+1,
		     (lastCol-firstCol)/colStride+1,
		     m,
		     m->offset +firstRow*m->rowStride +firstCol*m->columnStride,
		     rowStride*m->rowStride, colStride*m->columnStride );
    }
    
    //
    // implementations of LinearOperator methods
    //
    
    /** get number of rows */
    inline unsigned long getNumRows() const
    {
      return m->rows;
    }
    
    /** get number of columns */
    inline unsigned long getNumColumns() const
    {
      return m->columns;
    }
    
    /** right-multiply a vector */
    virtual Vector &rightMultiply( const Vector &v, Vector &result ) const;
    
    /** left-multiply a vector */
    virtual Vector &leftMultiply( const Vector &v, Vector &result ) const;
    
    //
    // initializations
    //
    
    /** set matrix to zero */
    inline void zero()
    {
      errorCond( m->data!= NULL, "  uninitialized matrix!" );
      
      unsigned long i;
      double *d= m->data + m->offset;
      
      // if the matrix is a contiguous memory region, use bzero(),
      // otherwise clear each row vector individually
      if( m->columnStride== 1ul && m->rowStride== m->columns )
	bzero( d, m->rows * m->columns * sizeof( double ) );
      else
	for( i= 0 ; i< m->rows ; i++ )
	  m->rowVec[i].zero();
    }
    
    /** set matrix to identity */
    inline void identity()
    {
      double *d= m->data + m->offset;
      unsigned long stride= m->rowStride + m->columnStride; // diagonal stride
      unsigned long i;
      
      zero();
      for( i= (m->rows< m->columns ? m->rows : m->columns) ;
	   i> 0 ; i--, d+= stride )
	(*d)= 1.0;
    }
    
    /** transpose THIS matrix (without copying elements) */
    inline Matrix &transpose()
    {
      unsigned long h;
      
      // swap dimensions, strides; then re-generate row vectors
      h= m->rows; m->rows= m->columns; m->columns= h;
      h= m->rowStride; m->rowStride= m->columnStride; m->columnStride= h;
      m->makeRowVecs();
      
      return *this;
    }
    
    /** create a NEW matrix with the transpose elements (sharing the
	original data) */
    inline Matrix getTranspose() const
    {
      return Matrix( m->columns, m->rows, m,
		     m->offset, m->columnStride, m->rowStride );
    }
    
    /** swap two rows */
    inline Matrix &swapRows( unsigned long r1, unsigned long r2 )
    {
      errorCond( r1< m->rows && r2< m->rows, "invalid row number" );
      double *d1=  m->data + m->offset + r1*m->rowStride;
      double *d2=  m->data + m->offset + r2*m->rowStride;
      double tmp;
      for( unsigned long i= 0 ; i< m->columns ;
	   i++, d1+= m->columnStride, d2+= m->columnStride )
      {
	tmp= *d1; *d1= *d2; *d2= tmp;
      }
      return *this;
    }
    
    /** swap two columns */
    inline Matrix &swapColumns( unsigned long r1, unsigned long r2 )
    {
      errorCond( r1< m->columns && r2< m->columns, "invalid column number" );
      double *d1=  m->data + m->offset + r1*m->columnStride;
      double *d2=  m->data + m->offset + r2*m->columnStride;
      double tmp;
      for( unsigned long i= 0 ; i< m->rows ;
	   i++, d1+= m->rowStride, d2+= m->rowStride )
      {
	tmp= *d1; *d1= *d2; *d2= tmp;
      }
      return *this;
    }
    
    // norms
    
    /** Frobenius norm (sum of squares of all entries) */
    double frobeniusNorm();
    
    
    // access operators
    
    /** set an element to a new value */
    inline void set( unsigned long r, unsigned long c, double val )
    {
      errorCond( m->data!= NULL, "  uninitialized matrix!" );
      
      m->data[r*m->rowStride + c*m->columnStride + m->offset]= val;
    }
    
    /** get the current value of an element */
    inline double get( unsigned long r, unsigned long c ) const
    {
      errorCond( m->data!= NULL, "  uninitialized matrix!" );
      
      return m->data[r*m->rowStride + c*m->columnStride + m->offset];
    }
    
    /** read/write access to row vectors */
    inline Vector& operator[]( unsigned long i )
    {
      errorCond( m->rowVec!= NULL, "  uninitialized matrix!" );
      
      return m->rowVec[i];
    }
    
    /** read-only access to row vectors */
    inline const Vector& operator[]( unsigned long i ) const
    {
      errorCond( m->rowVec!= NULL, "  uninitialized matrix!" );
      
      return m->rowVec[i];
    }
    
    /** row vector (precomputed; same as operator[]) */
    inline Vector &getRowVector( unsigned long i )
    {
      errorCond( m->rowVec!= NULL, "  uninitialized matrix!" );
      
      return m->rowVec[i];
    }
    
    /** column vector (creates a new vector that shares memory with
	the matrix) */
    inline Vector getColumnVector( unsigned long i ) const
    {
      errorCond( m->data!= NULL, "  uninitialized matrix!" );
      
      return Vector( m->rows, m,
		     i*m->columnStride + m->offset, m->rowStride );
    }
    
    /** diagonal vector (shares mamory with matrix) */
    inline Vector getDiagonal() const
    {
      errorCond( m->data!= NULL, "  uninitialized matrix!" );
      
      return Vector( m->rows< m->columns ? m->rows : m->columns, m,
		     m->offset, m->rowStride + m->columnStride );
    }
      
    // assignment-style operators
    
    /** multiplication with scalar (not threaded) */
    inline Matrix& operator*=( double s )
    {
      for( unsigned long i= 0 ; i< m->rows ; i++ )
	m->rowVec[i]*= s;
      return *this;
    }
    
    /** division by scalar (not threaded) */
    inline Matrix& operator/=( double s )
    {
      for( unsigned long i= 0 ; i< m->rows ; i++ )
	m->rowVec[i]/= s;
      return *this;
    }
    
    /** sum of two matrices (not threaded) */
    inline Matrix &operator+=( const Matrix &other )
    {
      errorCond( m->rows== other.m->rows, "  number of rows doesn't match" );
      
      for( unsigned long i= 0 ; i< m->rows ; i++ )
	m->rowVec[i]+= other.m->rowVec[i];
      return *this;
    }
    
    /** difference of two matrices (not threaded) */
    inline Matrix &operator-=( const Matrix &other )
    {
      errorCond( m->rows== other.m->rows, "  number of rows doesn't match" );
      
      for( unsigned long i= 0 ; i< m->rows ; i++ )
	m->rowVec[i]-= other.m->rowVec[i];
      return *this;
    }
    
    /** add the outer product of a vector with itself */
    inline Matrix &addOuterProduct( const Vector &v )
    {
      errorCond( m->rows== m->columns && m->rows== v.getSize(),
		 "  matrix and vector dimensions don't match" );
      
      double h;
      double *d= m->data + m->offset;
      for( unsigned long i= 0 ; i< m->rows ; i++ )
      {
	d[i * (m->rowStride + m->columnStride)]+= v[i] * v[i];
	for( unsigned long j= i+1; j< m->rows ; j++ )
	{
	  h= v[i] * v[j];
	  d[i*m->rowStride + j*m->columnStride]+= h;
	  d[j*m->rowStride + i*m->columnStride]+= h;
	}
      }
      return *this;
    }
    
  protected:
  
    /** \class Mat Matrix.hh
	actual matrix data with reference counter */
    class Mat: public MemoryObject {
    public:
      /** constructor */
      inline Mat( unsigned long _rows= 0ul, unsigned long _cols= 0ul,
		  bool fortranLayout= false )
	: MemoryObject(),
	  rowVec( NULL ), rows( _rows ), columns( _cols ), offset( 0ul ),
	  rowStride( fortranLayout ? 1ul : _cols ),
	  columnStride( fortranLayout ? _rows : 1ul )
      {
	allocate( rows*columns );
	makeRowVecs();
      }
      
      /** constructor using external memory object */
      inline Mat( unsigned long _rows, unsigned long _cols,
		  MemoryObject *owner, unsigned long _offset= 0ul,
		  unsigned long rStride= 0ul, unsigned long cStride= 1ul )
	: MemoryObject( owner->data, owner ),
	  rowVec( NULL ), rows( _rows ), columns( _cols ), offset( _offset ),
	  rowStride( rStride== 0 ? _cols : rStride ), columnStride( cStride )
      {
	// create row vectors
	makeRowVecs();
      }
      
      /** constructor using external data */
      inline Mat( unsigned long _rows, unsigned long _cols,
		  double *_data, unsigned long _offset= 0ul,
		  unsigned long rStride= 0ul, unsigned long cStride= 1ul )
	: MemoryObject( _data ),
	  rowVec( NULL ), rows( _rows ), columns( _cols ), offset( _offset ),
	  rowStride( rStride== 0 ? _cols : rStride ), columnStride( cStride )
      {
	// create row vectors
	makeRowVecs();
      }
      
      /** regenerate row vectors */
      inline void makeRowVecs()
      {
	if( rowVec!= NULL )
	  delete [] rowVec;
 	rowVec= new Vector[rows];
	for( unsigned long i= 0 ; i< rows ; i++ )
	  rowVec[i]= Vector( columns, this, i*rowStride+offset, columnStride );
      }
      
      /** overloaded deref operator to destroy row vectors if necessary */
      inline void deref()
      {
	// if the only references are from myself plus my row vectors, 
	// then I really am deleting this Mat
	if( refCounter==(rows+1) && rowVec!= NULL )
	  {
	    delete [] rowVec;
	    rowVec = NULL;
	  }
	MemoryObject::deref();
      }

      /** row vectors referencing into the data array (for performance
	  and convenience) */
      Vector *rowVec;
      
      /** offset from data to the beginning of the first vector element
       *  (i.e. v[0] is at data[offset])
       */
      unsigned long offset;
      
      /** stride between row elements */
      unsigned long rowStride;
      
      /** stride between column elements */
      unsigned long columnStride;
      
      /** number of rows */
      unsigned long rows;
      
      /** number of rows */
      unsigned long columns;
      
    protected:
      
      /** destructor -- private (use deref to delete references) */
      virtual ~Mat()
      {
	// rowVec should already be cleaned up now, but just in case I'll leave this code here
 	if( rowVec!= NULL )
 	  delete [] rowVec;
      }
      
    };
    
    /** reference to the actual matrix data */
    Mat *m;
    
  };

  
  //
  // methods involving Matrix
  //
  
  /** matrix-matrix product (with target matrix) */
  Matrix &multMatrixMatrix( const Matrix &m1, const Matrix &m2,
			    Matrix &result );
  
  
  //
  // specialized MetaData set/get functions
  //
  
  /** MetaData get for Vector objects */
  bool get( MetaData &m, char const *path, Matrix &val );
  
  /** MetaData set for Vector objects (currently only inline) */
  void set( MetaData &m, char const *path, Matrix &val );

  
} /* namespace */

#endif /* LINEARALGEBRA_MATRIX_H */

