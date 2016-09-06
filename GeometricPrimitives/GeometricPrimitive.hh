// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: Felix Heide
// Email:   heide@informatik.uni-siegen.de
// ==========================================================================

#ifndef GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVE_H
#define GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVE_H

/*! \file  GeometricPrimitive.hh
    \brief An abstract base class for geometric primitives in n dimensions.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>
#include "MDA/LinearAlgebra/Vector.hh"

namespace MDA {

  using namespace std;

  /** \class GeometricPrimitive GeometricPrimitive.hh
      An abstract base class for geometric primitives in n dimensions. */
  class GeometricPrimitive {

  public:

    /** returns true if the point is inside the primitive, else false */
    virtual bool inside( Vector &point )= 0;

     /** destructor */
    virtual ~GeometricPrimitive() {}
  };


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVE_H */

