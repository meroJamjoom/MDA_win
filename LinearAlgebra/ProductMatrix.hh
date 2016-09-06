// ==========================================================================
// $Id: ProductMatrix.hh 319 2009-05-27 21:17:01Z heidrich $
// product of two matrices (without explicitly computing&storing the product)
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

#ifndef LINEARALGEBRA_PRODUCTMATRIX_H
#define LINEARALGEBRA_PRODUCTMATRIX_H

/*! \file  ProductMatrix.hh
    \brief product of two matrices (without explicitly
    computing&storing the product)
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/Errors.hh"

#include "LinearOperator.hh"

namespace MDA {

  /** \class ProductMatrix ProductMatrix.hh
      product of two matrices (without explicitly
      computing&storing the product)
  
      recursive product matrices are possible, but the tmeplate needs
      to be instantiated elsewhere...
  */
  class ProductMatrix: public LinearOperator {
    
  public:
    
    /** default constructor */
    inline ProductMatrix()
    {}
    
    /** constructor */
    inline ProductMatrix( LinearOperator *left, LinearOperator *right )
      : lmat( left ), rmat( right )
    {
      errorCond( lmat->getNumColumns()== rmat->getNumRows(),
		 "  matrix dimensions do not match!" );
    }
    
    /** get number of rows */
    inline virtual unsigned long getNumRows() const
    {
      return lmat->getNumRows();
    }
    
    /** get number of columns */
    inline virtual unsigned long getNumColumns() const
    {
      return rmat->getNumColumns();
    }
    
    /** right-multiply a vector */
    inline virtual Vector &rightMultiply( const Vector &v, Vector &result )
    {
      Vector tmp( lmat->getNumColumns() );
      
      rmat->rightMultiply( v, tmp );
      lmat->rightMultiply( tmp, result );
  
      return result;
    }
    
    /** left-multiply a vector */
    inline virtual Vector &leftMultiply( const Vector &v, Vector &result )
    {
      Vector tmp( lmat->getNumColumns() );
      
      lmat->leftMultiply( v, tmp );
      rmat->leftMultiply( tmp, result );
  
      return result;
    }
    
  protected:
    
    /** left matrix */
    LinearOperator *lmat;
    
    /** right matrix */
    LinearOperator *rmat;
    
  };

} /* namespace */

#endif /* LINEARALGEBRA_PRODUCTMATRIX_H */

