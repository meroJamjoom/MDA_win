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

#ifndef LINEARALGEBRA_VECTOR_C
#define LINEARALGEBRA_VECTOR_C

#include <MDA/Config.hh>
#include "Vector.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** copy a vector into a newly created memory object, which becomes
    the property of this vector */
void
Vector::copy( const Vector &orig )
{
  // create new Vec object (dereference old)
  v->deref();
  v= new Vec( orig.v->size );
  
  // use the assign method to actually move the data
  assign( orig );
}
  
/** assign the values of another vector to the memory currently
    held by this vector (dimensions must match) */
Vector &
Vector::assign( const Vector & orig )
{
  unsigned long i;
  
  // check for matching dimensions
  errorCond( v->size== orig.v->size, "  vector dimensions do not match!" );
  
  // then copy the data
  double *d= v->data + v->offset;
  double *s= orig.v->data + orig.v->offset;
  if( v->stride== 1ul && orig.v->stride== 1ul )
    // memcopy if both vectors have simple, contiguous memory layout
    memcpy( d, s, v->size*sizeof( double ) );
  else
    // element-by-element copy, otherwise
    for( i= v->size ; i> 0ul ; i--, d+= v->stride, s+= orig.v->stride )
      (*d)= (*s);
  
  return *this;
}

/** Lp norm for arbitrary p (p<= 0.01 is minimum semi-norm, p>100 corrensponds to Linfinity)*/
double
Vector::norm( double p ) const
{
  unsigned long i;
  
  errorCond( v->data!= NULL, "  uninitialized vector!" );
  
  double result= 0.0;
  double *d= v->data + v->offset;
  
  // for very small p, the Lp norm converges to the minimum of the
  // components, while for large p, it converges to the maximum
  // (i.e. the Linfinity norm).
  //
  // Since the standard Lp formulas get instable for exteme values of
  // p, we use special cases for p<LP_NORM_THRESHOLD and p>1/LP_NORM_THRESHOLD
  if( p< LP_NORM_THRESHOLD )
  {
    // L0 = min of all components - this is a seminorm (many vectors
    // have zero distance!)
    result= DBL_MAX;
    for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
      result= result< fabs(*d) ? result : fabs(*d);
    return result;
  }
  else if( p> 1.0/LP_NORM_THRESHOLD )
  {
    // Linfinity = max of all components
    for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
      result= result> fabs(*d) ? result : fabs(*d);
    return result;
  }
  
  // everything else is a regular Lp norm
  for( unsigned long i= v->size ; i> 0 ; i--, d+= v->stride )
    result+= pow( fabs( *d ), p );
  return pow( result, 1.0/p );
}


//
// methods
//


/** dot product of two vectors (generic version) */
double
dot( const Vector &v1, const Vector &v2 )
{
  unsigned long size= v1.getSize();
  errorCond( v2.getSize()== size, "  vector dimensions must match" );

  double result= 0.0;
  for( unsigned long i= 0 ; i< size ; i++ )
    result+= v1[i] * v2[i];

  return result;
}

/** Euclidean distance squared */
double
distSq( const Vector &v1, const Vector &v2 )
{
  unsigned long size= v1.getSize();
  errorCond( size== v2.getSize(), "  vector dimensions must match" );
  
  double h, result= 0.0;
  for( unsigned long i= 0 ; i< size ; i++ )
  {
    h= v1[i] - v2[i];
    result+= h*h;
  }
  return result;
}

/** Lp distance */
double
dist( const Vector &v1, const Vector &v2, double p )
{
  unsigned long size= v1.getSize();
  errorCond( size== v2.getSize(), "  vector dimensions must match" );
  
  double result= 0.0;
  for( unsigned long i= 0 ; i< size ; i++ )
    result+= pow( fabs( v1[i] - v2[i] ), p );
  return pow( result, 1.0/p );
}


//
// MetaData set/get
//

/** specialized get for Vector objects */
bool
get( MetaData &m, char const *path, Vector &val )
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
    
    // first, get the size of the vector
    unsigned length;
    if( !warnCond( get( m, (p+":length").c_str(), length ),
		   "missing length attribute in inlined Vector") )
      return false;
    
    // create vector
    val= Vector( length );
    // then read elements as raw array of doubles
    return get( m, path, &val[0], length );
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
    
    val= Vector( array, channel );
    return true;
  }
  
  // if we reach this, things have gone wrong
  warning( "unknown storage type" );
  return false;
}


/** specialized set for Vector objects */
void
set( MetaData &m, char const *path, Vector &val )
{
  string p( path );
  
  // set the length attribute and then just use the raw array write
  unsigned length= val.getSize();
  set( m, (p+":length").c_str(), length );
  set( m, path, &val[0], length );
}


} /* namespace */

#endif /* LINEARALGEBRA_VECTOR_C */

