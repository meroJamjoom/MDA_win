// ==========================================================================
// $Id:$
// a n-dimensional vector
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

#ifndef LINEARALGEBRA_VECTOR_H
#define LINEARALGEBRA_VECTOR_H

/*! \file  Vector.hh
    \brief a n-dimensional vector
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <string.h>
#include <math.h>
#include <float.h>

#include "MDA/Base/Errors.hh"
#include "MDA/Base/MetaData.hh"
#include "MDA/Array/Array.hh"

#if defined(_WIN32) || defined(_WIN64)
#define bzero(p, l) memset(p, 0, l)
#endif

namespace MDA {
  
  using namespace std;
  
  /** \class Vector Vector.hh
      a n-dimensional vector */
  class Vector {
    
  public:
    
    /** constructor (no initialization but with memory allocation if
	the number of elements is >0) */
    inline Vector( unsigned long n= 0ul )
    {
      v= new Vec( n );
    }
    
    /** constructor (with initialization) */
    inline Vector( unsigned long n, double val )
    {
      v= new Vec( n );
      if( val== 0.0 )
	zero();
      else
	for( unsigned long i= 0 ; i< n ; i++ )
	  v->data[i]= val;
    }
    
    /** contructor for vector sharing data with another memory object */
    inline Vector( unsigned long n, MemoryObject *owner,
		   unsigned long offset= 0ul, unsigned long stride= 1ul )
    {
      v= new Vec( n, owner, offset, stride );
    }
    
    /** contructor from a channel in an array */
    inline Vector( Array<double> &a, unsigned ch )
    {
      CoordinateVector dim= a.getDimension();
      unsigned long size= 1;
      for( int i= dim.vec.size()-1 ; i>= 0 ; i-- )
	size*= dim.vec[i];
      
      v= new Vec( size, a[ch], 0ul, 1ul );
    }
    
    /** contructor from external raw data (this data is NOT copied,
	and changes to the vector will overwrite it!) */
    inline Vector( unsigned long n, double *rawData,
		   unsigned long offset= 0ul, unsigned long stride= 1ul )
    {
      v= new Vec( n, rawData, offset, stride );
    }
    
    /** copy constructor */
    inline Vector( const Vector &vec )
    {
      v= (Vec *)(vec.v->ref());
    }
    
    /** assignment operator (adds a reference to an existing vector -
	for a new physical copy use the copy method) */
    inline const Vector &operator=( const Vector &vec )
    {
      // de-reference current vector
      v->deref();
      // and reference the target instead
      v= (Vec *)(vec.v->ref());
      return *this;
    }
    
    /** copy a vector into a newly created memory object, which
	becomes the property of this vector */
    void copy( const Vector &orig );
    
    /** assign the values of another vector to the memory currently
	held by this vector (dimensions must match) */
    Vector &assign( const Vector &orig );
    
    /** destructor */
    inline ~Vector()
    {
      v->deref();
    }
    
    /** report size of vector */
    inline unsigned long getSize() const
    {
      return v->size;
    }
    
    /** report estimate of non-zero elements - same as vector size for
	a non-sparse Vector */
    inline unsigned long getNumNonZero() const
    {
      return v->size;
    }
    
    /** set all vector components to zero */
    inline void zero()
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      unsigned long i;
      double *d= v->data + v->offset;
      
      if( v->stride== 1ul )
	bzero( d, v->size*sizeof( double ) );
      else
	for( i= v->size ; i> 0 ; i--, d+= v->stride )
	  (*d)= 0.0;
	}
    
    /** randomize each component in a range of 0..1 */
	inline void randomize()
	{
		errorCond(v->data != NULL, "  uninitialized vector!");

		double *d = v->data + v->offset;
		for (unsigned long i = v->size; i > 0; i--, d += v->stride)

#if defined (_WIN32) || defined(_WIN64)
#define drand48() (rand()/(RAND_MAX+1.0))
#else
			*d= drand48();
#endif 

		;
	}
    // norms
    
    /** L2 norm squared */
    inline double normSq() const
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      double result= 0.0;
      double *d= v->data + v->offset;
      for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
	result+= (*d) * (*d);
      return result;
    }
    
    /** L2 norm */
    inline double norm() const
    {
      return sqrt( normSq() );
    }
    
    /** Lp norm for arbitrary p.
	p>100 corrensponds to max norm (Linfinity), p<1/100 to min norm */
    double norm( double p ) const;
    
    /** normalize nonzero vector */
    inline void normalize()
    {
      double length= norm();
      if( length> DBL_EPSILON )
	*this/= length;
    }
    
    // access operators
    
    /** get the value of an entry */
    inline double get( unsigned long pos ) const
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      return v->data[pos * v->stride + v->offset];
    }
    
    /** set an entry */
    inline void set( unsigned long pos, double value )
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      v->data[pos * v->stride + v->offset]= value;
    }
    
    /** read/write element access */
    inline double& operator[]( unsigned long i )
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      return v->data[i * v->stride + v->offset];
    }
    
    /** read-only element access */
    inline double operator[]( unsigned long i ) const
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      return v->data[i * v->stride + v->offset];
    }
    
    /** create a subvector that shares memory */
    inline Vector getSubVector( unsigned long first, unsigned long last,
				unsigned long stride= 1 )
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      return Vector( (last-first)/stride+1, v,
		     v->offset + first * v->stride,
		     stride * v->stride );
    }
    
    // assignment-style operators
    
    /** multiplication with scalar */
    inline Vector& operator*=( double s )
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      double *d= v->data + v->offset;
      for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
	(*d)*= s;
      return *this;
    }
    
    /** division by scalar */
    inline Vector& operator/=( double s )
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      double *d= v->data + v->offset;
      for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
	(*d)/= s;
      return *this;
    }
    
    /** negate (same as vector*= -1.0) */
    inline Vector& negate()
    {
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      double *d= v->data + v->offset;
      for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
	(*d)= -(*d);
      return *this;
    }
      
    /** sum of two vectors */
    inline Vector &operator+=( const Vector &other )
    {
      errorCond( v->size== other.v->size, "  vector dimensions don't match" );
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      unsigned long i;
      double *d= v->data + v->offset;
      double *s= other.v->data + other.v->offset;
      for( i= v->size ; i> 0 ; i--, d+= v->stride, s+= other.v->stride )
	(*d)+= (*s);
      return *this;
    }
    
    /** difference of two vectors */
    inline Vector &operator-=( const Vector &other )
    {
      errorCond( v->size== other.v->size, "  vector dimensions don't match" );
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      unsigned long i;
      double *d= v->data + v->offset;
      double *s= other.v->data + other.v->offset;
      for( i= v->size ; i> 0 ; i--, d+= v->stride, s+= other.v->stride )
	(*d)-= (*s);
      return *this;
    }
    
    /** (*this) += sc * vec */
    inline Vector &addScalarTimesVector( double sc, const Vector &other )
    {
      errorCond( v->size== other.v->size, "  vector dimensions don't match" );
      errorCond( v->data!= NULL, "  uninitialized vector!" );
      
      unsigned long i;
      double *d= v->data + v->offset;
      double *s= other.v->data + other.v->offset;
      for( i= v->size ; i> 0 ; i--, d+= v->stride, s+= other.v->stride )
	(*d)+= sc * (*s);
      return *this;
    }
    
  protected:
    
    /** \class Vec Vector.hh
	actual vector data with reference counter */
    class Vec: public MemoryObject {
    public:
      /** constructor - new vector of certain length */
      inline Vec( unsigned long _size= 0ul )
	: MemoryObject(),
	  size( _size ), offset( 0ul ), stride( 1ul )
      {
	allocate( size );
      }

      /** constructor using data from another linear algebra object */
      inline Vec( unsigned long _size, MemoryObject *other,
		  unsigned long _offset= 0ul, unsigned long _stride= 1ul )
	: MemoryObject( other->data, other ),
	  size( _size ), offset( _offset ), stride( _stride )
      {}
      
      /** constructor using completely external data */
      inline Vec( unsigned long _size, double *_data,
		  unsigned long _offset= 0ul, unsigned long _stride= 1ul )
	: MemoryObject( _data ),
	  size( _size ), offset( _offset ), stride( _stride )
      {}
      
      /** offset from data to the beginning of the first vector element
       *  (i.e. v[0] is at data[offset])
       */
      unsigned long offset;
      
      /** stride between data elements
       *  (i.e. if v[i] is at data[j], then v[i+1] is at data[j+stride])
       */
      unsigned long stride;
      
      /** vector length */
      unsigned long size;
      
    protected:
      
      /** destructor -- private (use deref to delete references) */
      virtual ~Vec()
      {}
      
    };
    
    /** reference to the actual vector data */
    Vec *v;
    
  };
  
  
  //
  // methods involving vectors
  //
  
  /** dot product of two vectors */
  double dot( const Vector &v1, const Vector &v2 );
  
  /** Euclidean distance squared */
  double distSq( const Vector &v1, const Vector &v2 );
  
  /** Euclidean distance */
  inline double dist( const Vector &v1, const Vector &v2 )
  {
    return sqrt( distSq( v1, v2 ) );
  }
  
  /** Lp distance */
  double dist( const Vector &v1, const Vector &v2, double p );
  
  //
  // specialized MetaData set/get functions
  //
  
  /** MetaData get for Vector objects */
  bool get( MetaData &m, char const *path, Vector &val );
  
  /** MetaData set for Vector objects (currently only inline) */
  void set( MetaData &m, char const *path, Vector &val );
  
  
} /* Namespace */

#endif /* LINEARALGEBRA_VECTOR_H */

