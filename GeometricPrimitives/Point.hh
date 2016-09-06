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

#ifndef GEOMETRICPRIMITIVES_POINT_H
#define GEOMETRICPRIMITIVES_POINT_H

/*! \file  Point.hh
    \brief A class for a Point in n dimensions
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/LinearAlgebra/Vector.hh"
#include "GeometricPrimitive.hh"

namespace MDA {

  using namespace std;

  /** \class Point Point.hh
      A class for a Point in n dimensions*/

  class Point : public GeometricPrimitive {

  public:

    /** default constructor, point is initialized with empty position */
    Point();

    /** point is initialized as coordinate system origin in parameterized dimensions */
    Point( unsigned long dimensions );

    /** explicit point constructor */
    Point( const Vector &point);

    /** default destructor */
    ~Point();

    /** returns always false, as a point is infinitively thin */
    virtual bool inside( Vector &point );

    /** Position of the point */
    Vector pos;

  };

  //
  // methods involving Points
  //
  /** distance between two points */
  double dist( const Point &p1, const Vector &p2);

} /* namespace */

#endif /* GEOMETRICPRIMITIVES_POINT_H */

