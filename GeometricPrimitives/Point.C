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

#ifndef GEOMETRICPRIMITIVES_POINT_C
#define GEOMETRICPRIMITIVES_POINT_C

#include "Point.hh"
#include <math.h>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  /** default constructor, point is initialized with empty position */
  Point::Point(){}

  /** point is initialized as coordinate system origin in parameterized dimensions */
  Point::Point( unsigned long dimensions ) : pos( dimensions, 0.0 ){}

  /** explicit point constructor */
  Point::Point( const Vector &point) { pos.copy( point ); }

  /** default destructor */
  Point::~Point(){}

  /** returns always false, as a point is infinitively thin */
  bool Point::inside( Vector &point )
  {
      return false;
  }
    
  //
  // methods involving Points
  //
  /** distance between two points */
  double dist( const Point &p1, const Vector &p2)
  {
      return dist( p1.pos, p2 );
  }

} /* namespace */

#endif /* GEOMETRICPRIMITIVES_POINT_C */

