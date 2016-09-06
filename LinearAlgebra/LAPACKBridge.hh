// ==========================================================================
// $Id:$
// bridge to LAPACK routines
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

#ifndef LAPACK_LAPACKBRIDGE_H
#define LAPACK_LAPACKBRIDGE_H

/*! \file  LAPACKBridge.hh
    \brief bridge to LAPACK routines
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/Errors.hh"
#include "MDA/LinearAlgebra/LinAlg.hh"

namespace MDA {
  
  /** inverse of a square matrix
      \param m: square matrix to be inverted
      \param preserveMatrix: whether the contents of m must be preserved
      \returns the inverse of m
  */
  Matrix inverse( Matrix &m, bool preserveMatrix= false );
  
  /** singular value decomposition
      \param m: input matrix
      \param U: matrix left singular vectors (output)
      \param Vt: (transpose of) matrix for right singular vectors (output)
      \returns vector of singular values
  */
  Vector svd( Matrix &m, Matrix &U, Matrix &Vt );
  
} /* namespace */

#endif /* LAPACK_LAPACKBRIDGE_H */

