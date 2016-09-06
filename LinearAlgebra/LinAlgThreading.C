// ==========================================================================
// $Id: LinAlgThreading.C 292 2009-03-20 03:08:52Z heidrich $
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

#ifndef LINEARALGEBRA_LINALGTHREADING_C
#define LINEARALGEBRA_LINALGTHREADING_C

#include "LinAlgThreading.hh"
#include "Vector.hh"
#include "Matrix.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** constructor
 * \param v1: vector 1
 * \param v2: vector 2
 * \param result: result
 * \param addToResult: whether the result of the dot product is
 * added to the current value of the result, or just overwrites
 * it. Adding is implemented as a reduction operator.
 */
DotProductJob::DotProductJob( const Vector &_v1, const Vector &_v2,
			      double &_result, bool addToResult )
  : SMPJob( _v1.getNumNonZero()*_v2.getNumNonZero() ),
    v1( _v1 ), v2( _v2 ), result( _result )
{
  SMPJob::applyReduction= addToResult;
}

/** execute a dot product */
void
DotProductJob::execute( int jobID )
{
  tmpResult= dot( v1, v2 );
  if( !applyReduction )
    result= tmpResult;
}

/** reduction operator adds dot product result to existing value */
void
DotProductJob::reduce( int jobID )
{
  result+= tmpResult;
}




} /* namespace */

#endif /* LINEARALGEBRA_LINALGTHREADING_C */

