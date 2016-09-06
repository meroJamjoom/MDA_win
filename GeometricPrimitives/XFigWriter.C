// ==========================================================================
// $Id: XFigWriter.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef GEOMETRICPRIMITIVES_XFIGWRITER_C
#define GEOMETRICPRIMITIVES_XFIGWRITER_C

#if defined (_WIN32) || defined (_WIN64)
#define _USE_MATH_DEFINES
#include <math.h>
#else 
#include <math.h>
#endif

#include "XFigWriter.hh"




namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


const XFigWriter::Style XFigWriter::defaultStyle= {
  50,	// depth
  0,	// line style = solid
  1,	// thickness = 1/80"
  -1,	// default pen color
  -1,	// default fill color
  -1,	// area fill = off
  1.0,	// space between dashes in 1/80"
  1,	// round line joints
  1,	// round line caps
  10,	// arc box radius
  0,	// forward arrow off
  0,	// backward arrow off
  1,	// triangular arrow heads
  1,	// filled arrow heads
  1.0,	// arrow thickness in 1/80"
  8.0,	// arrow width in effective DPI units
  8.0,	// arrow height in effective DPI units
};


/** open output stream */
bool
XFigWriter::open( ostream &_os )
{
  os= &_os;
  os->clear();
  ownOstream= false;
  if( os->good() )
    return writeHeader();
  return false;
}
  
/** open output file */
bool
XFigWriter::open( char *fileName )
{
  myOstream.open( fileName, ofstream::out & ofstream::binary );
  ownOstream= myOstream.good();
  if( ownOstream )
  {
    os= &myOstream;
    return writeHeader();
  }
  
  return false;
}

/** write header */
bool
XFigWriter::writeHeader()
{
  (*os) << "#FIG 3.2\n"
	<< "Portrait\n"
	<< "Center\n"
	<< "Inches\n"
	<< "Letter\n"
	<< "1\n"
	<< "Single\n"
	<< "-2\n"
	<< XFIG_INTERNAL_DPI << " 2\n";
  return os->good();
}

/** close stream */
bool
XFigWriter::close()
{
  if( ownOstream )
  {
    myOstream.close();
    ownOstream= false;
  }
  os= NULL;
  return true;
}


  //
  // primitives
  //
  
/** a single line segment */
void
XFigWriter::drawLine( int p1x, int p1y, int p2x, int p2y )
{
  errorCond( os!= NULL, "No file open" );
  
  // base options for polylines
  (*os) << 2 << ' '			// polyline object
	<< 1 << ' '			// polyline subobject
	<< style.lineStyle << ' '
	<< style.thickness << ' '
	<< style.penColor << ' '
	<< style.fillColor << ' '
	<< style.depth << ' '
	<< 0 << ' '			// unused: pen style
	<< style.areaFill << ' '
	<< style.dashStyle << ' '
	<< style.jointStyle << ' '
	<< style.capStyle << ' '
	<< style.arcBoxRadius << ' '
	<< style.forwardArrow << ' '
	<< style.backwardArrow << ' '
	<< 2 << endl;			// number of points on polyline
  
  // defnitition of forward arrow head
  if( style.forwardArrow )
    (*os) << '\t' << style.arrowStyle
	  << ' ' << style.arrowFill
	  << ' ' << style.arrowThickness
	  << ' ' << style.arrowWidth/dpiMult
	  << ' ' << style.arrowHeight/dpiMult << endl;
  
  // definition of backward arrow head
  if( style.backwardArrow )
    (*os) << '\t' << style.arrowStyle
	  << ' ' << style.arrowFill
	  << ' ' << style.arrowThickness
	  << ' ' << style.arrowWidth/dpiMult
	  << ' ' << style.arrowHeight/dpiMult << endl;
  
  // endpoints
  (*os) << '\t' << p1x << ' ' << p1y
	<< ' ' << p2x << ' ' << p2y
	<< endl;
}


/** an ellipse */
void
XFigWriter::drawEllipse( int type, int cx, int cy, int r1, int r2, double angle)
{
  errorCond( os!= NULL, "No file open" );
  
  (*os) << 1 << ' '			// ellipse type
	<< type << ' '			// use specified sub-type
	<< style.lineStyle << ' '
	<< style.thickness << ' '
	<< style.penColor << ' '
	<< style.fillColor << ' '
	<< style.depth << ' '
	<< 0 << ' '			// unused: pen style
	<< style.areaFill << ' '
	<< style.dashStyle << ' '
	<< 1 << ' '			// unused: direction
	<< angle << ' '
	<< cx << ' '			// center x
	<< cy << ' '			// center y
	<< r1 << ' '			// radius along first axis
	<< r2 << ' '			// radius along second axis
	<< cx << ' '			// first point entered, x
	<< cy << ' ';			// first point entered, y
  
  // last point semantics are different for ellipses and points
  if( type== 1 )
    // for ellipses, it is the corner of a bounding box
    (*os) << lround( cx+r1*cos(angle)+r2*cos(angle+M_PI_2) ) << ' '
	  << lround( cy-r1*sin(angle)-r2*sin(angle+M_PI_2) ) << endl;
  else
    // for circles, it is an actual point on the circle
    (*os) << cx + r1 << ' ' << cy << endl;
}


/** an image (from an external image file) */
void
XFigWriter::drawImage( const char *fileName, int width, int height,
		       int xOff, int yOff )
{
  errorCond( os!= NULL, "No file open" );
  
  // basic style options for image primitive
  (*os) << 2 << ' '			// polyline object
	<< 5 << ' '			// image subobject
	<< style.lineStyle << ' '
	<< style.thickness << ' '
	<< style.penColor << ' '
	<< style.fillColor << ' '
	<< style.depth << ' '
	<< 0 << ' '			// unused: pen style
	<< style.areaFill << ' '
	<< style.dashStyle << ' '
	<< style.jointStyle << ' '
	<< style.capStyle << ' '
	<< style.arcBoxRadius << ' '
	<< 0 << ' ' << 0 << ' '		// no arrow heads for images...
	<< 5 << endl;			// number of points (one corner twice)
  
  // image name and corner points
  (*os) << '\t' << 0			// don't flip image
	<< ' ' << fileName << endl
	<< '\t' << xOff << ' ' << yOff
	<< ' ' << xOff+width << ' ' << yOff
	<< ' ' << xOff+width << ' ' << yOff+height
	<< ' ' << xOff << ' ' << yOff+height
	<< ' ' << xOff << ' ' << yOff << endl;
}


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_XFIGWRITER_C */

