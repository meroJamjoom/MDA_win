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

#ifndef GEOMETRICPRIMITIVES_LINE_H
#define GEOMETRICPRIMITIVES_LINE_H

/*! \file  Line.hh
    \brief A simple class for a line in n dimensions.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "GeometricPrimitive.hh"

namespace MDA {

  using namespace std;

  /** \class Line Line.hh
      A simple class for a line in n dimensions. */

  class Line : public GeometricPrimitive  {

  public:

    /** default constructor, line is empty*/
    Line();

    /** line is first axis of coordinate system with parameterized dimensions*/
    Line( unsigned long dimensions);

    /** a line defined by two points on it*/
    Line( Vector &p1, Vector &p2 );

    /** default destructor*/
    ~Line();

    /** returns always false, as a line is infinitively thin*/
    virtual bool inside( Vector &point );

    /** The line is defined by two points on it */
    Vector p1;
    Vector p2;

  };

  //
  // methods involving Lines
  //
  /** distance between a line and a point */
  double dist( const Line &line, const Vector &point);


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_LINE_H */

