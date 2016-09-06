// ==========================================================================
// $Id: Preconditioner.C 292 2009-03-20 03:08:52Z heidrich $
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

#ifndef LINEARALGEBRA_PRECONDITIONER_C
#define LINEARALGEBRA_PRECONDITIONER_C

#include "Preconditioner.hh"
#include "Matrix.hh"
#include "Vector.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** right multiply the preconditioner matrix P^{-1} with a vector */
void
Preconditioner::rightMultiply( const Vector &v, Vector &result )
{
  result.copy( v );
}


} /* namespace */

#endif /* LINEARALGEBRA_PRECONDITIONER_C */

