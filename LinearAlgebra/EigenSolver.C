// ==========================================================================
// $Id: EigenSolver.C 499 2009-11-24 23:43:58Z bradleyd $
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

#ifndef LINEARALGEBRA_EIGENSOLVER_C
#define LINEARALGEBRA_EIGENSOLVER_C

#include "EigenSolver.hh"
#include <algorithm>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// local comparison method for EigenValues
// (this is the ">" function, so large EVs are always first)
bool
compare( const EigenSolver::EigenValue &ev1,
	 const EigenSolver::EigenValue &ev2 )
{
  return ev1.value> ev2.value;
}


/** create the sorted list of eigen values */
void
EigenSolver::sort( Vector &evs )
{
  unsigned long i;
  
  // create the unsorted array
  unsigned long numEntries= evs.getSize();
  sortedEigenValues.resize( numEntries );
  for( i= 0 ; i< numEntries ; i++ )
  {
    sortedEigenValues[i].value= evs[i];
    sortedEigenValues[i].index= i;
  }
  
  // then use the vector sort
  std::sort( sortedEigenValues.begin(), sortedEigenValues.end(), compare );
}


} /* namespace */

#endif /* LINEARALGEBRA_EIGENSOLVER_C */

