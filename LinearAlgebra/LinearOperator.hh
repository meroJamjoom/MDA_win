// ==========================================================================
// $Id: LinearOperator.hh 319 2009-05-27 21:17:01Z heidrich $
// abstract base class for various matrix representations
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

#ifndef LINEARALGEBRA_LINEAROPERATOR_H
#define LINEARALGEBRA_LINEAROPERATOR_H

/*! \file  LinearOperator.hh
    \brief abstract base class for various matrix representations
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Vector.hh"

namespace MDA {

  /** \class LinearOperator LinearOperator.hh
      abstract base class for various matrix representations      

      A linear system must support the following operations 
      \item report number of rows and columns
      \item multiply a Vector from the right
      \item multiply a Vector from the left
  */
  
  class LinearOperator {

  public:

    /** get number of rows */
    virtual unsigned long getNumRows() const= 0;
    
    /** get number of columns */
    virtual unsigned long getNumColumns() const= 0;
    
    /** right-multiply a vector */
    virtual Vector &rightMultiply( const Vector &v, Vector &result ) const= 0;
    
    /** left-multiply a vector */
    virtual Vector &leftMultiply( const Vector &v, Vector &result ) const= 0;
    
    /** convenience function for right multiply with new result vector */
    inline Vector rightMult( const Vector &v ) const
    {
      Vector result( this->getNumRows() );
      return this->rightMultiply( v, result );
    }
    
    /** convenience function for left multiply with new result vector */
    inline Vector leftMult( const Vector &v ) const
    {
      Vector result( this->getNumColumns() );
      return this->leftMultiply( v, result );
    }
    
  };


} /* namespace */

#endif /* LINEARALGEBRA_LINEAROPERATOR_H */

