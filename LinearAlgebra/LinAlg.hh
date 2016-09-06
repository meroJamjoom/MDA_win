
/** \file LinAlg.hh
    \brief meta-include for basic linear algebra classes
*/

// basic types
#include "Vector.hh"
#include "Matrix.hh"
#include "SparseMatrix.hh"
#include "ProductMatrix.hh"

// LAPACK
#include "LAPACKBridge.hh"

// linear solvers & preconditioners
#include "JacobiPreconditioner.hh"
#include "ConjugateGradient.hh"

// Eigenvector computation
#include "JacobiRotation.hh"
