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

#ifndef GEOMETRICPRIMITIVES_LINE_C
#define GEOMETRICPRIMITIVES_LINE_C

#include "Line.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  /** default constructor, line is empty*/
  Line::Line(){}

  /** line is first axis of coordinate system with parameterized dimensions*/
  Line::Line( unsigned long dimensions): p1(dimensions, 0.0) , p2(dimensions, 0.0)
  {
      p2[0] = 1.0;
      return;
  }

  /** a line defined by two points on it*/
  Line::Line( Vector &p1, Vector &p2 )
  {
      //Exit, if wrong size
      if( p1.getSize() != p2.getSize() )
      {
	 cerr << "Vector dimensions must match! " << std::endl;
	 return;
      }

      (this->p1).copy( p1 );
      (this->p2).copy( p2 );

      return;
  }

  /** default destructor */
  Line::~Line(){}

  /** returns allways false, as a line is infinitively thin*/
  bool Line::inside( Vector &point )
  {
      return false;
  }

  //
  // methods involving Lines
  //
  /** distance between a line and a point */
  double dist( const Line &line, const Vector &point)
  {
      //Express the nearest point p on the line as p = line.p1 + u*(line.p2-line.p1)
      //then it has to be (point - p)dot(line.p2 - line.p1) = 0, and consequently:
      Vector lineVec( 2 );
      lineVec.copy( line.p2 );
      lineVec -= line.p1;
      Vector pointVec( 2 );
      pointVec.copy(  point );
      pointVec -= line.p1;
      
      double u = dot( pointVec, lineVec ) / dot( lineVec, lineVec );
      Vector p( 2 );
      p.copy( line.p1 );
      lineVec *= u;
      p += lineVec;
      return dist(point, p);
  }

} /* namespace */

#endif /* GEOMETRICPRIMITIVES_LINE_C */

