// ==========================================================================
// $Id: LinAlgThreading.hh 292 2009-03-20 03:08:52Z heidrich $
// symmetric multiprocessing jobs for common linear algebra operations
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

#ifndef LINEARALGEBRA_LINALGTHREADING_H
#define LINEARALGEBRA_LINALGTHREADING_H

/*! \file  LinAlgThreading.hh
    \brief symmetric multiprocessing jobs for common linear algebra operations
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Threading/SMPJob.hh"

#include "Vector.hh"

namespace MDA {

  class Vector;
  class Matrix;
  
  /** \class DotProductJob LinAlgThreading.hh
      an SMP job for performing a dot product of two vectors */
  class DotProductJob: public SMPJob {

  public:

    /** constructor
     * \param v1: vector 1
     * \param v2: vector 2
     * \param result: result
     * \param addToResult: whether the result of the dot product is
     * added to the current value of the result, or just overwrites
     * it. Adding is implemented as a reduction operator.
     */
    DotProductJob( const Vector &_v1, const Vector &_v2,
		   double &_result, bool addToResult= false );
    
    /** execute a dot product */
    virtual void execute( int jobID );
    
    /** reduction is addition */
    virtual void reduce( int jobID );
    
  protected:
    
    /** vector 1 */
    Vector v1;
    
    /** vector 2 */
    Vector v2;
    
    /** result */
    double &result;
    
    /** tmp result (before reduction) */
    double tmpResult;
    
  };

} /* namespace */

#endif /* LINEARALGEBRA_LINALGTHREADING_H */

