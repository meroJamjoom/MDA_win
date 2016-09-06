// ==========================================================================
// $Id: XFigWriter.hh 638 2010-03-08 19:34:24Z heidrich $
// Writer for creatign an XFig fiel from Ggeometric primitives
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

#ifndef GEOMETRICPRIMITIVES_XFIGWRITER_H
#define GEOMETRICPRIMITIVES_XFIGWRITER_H

/*! \file  XFigWriter.hh
    \brief Writer for creating an XFig file from Ggeometric primitives
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <fstream>

#include "MDA/Config.hh"
#include "MDA/LinearAlgebra/Vector.hh"

namespace MDA {

  using namespace std;
  
  /** \class XFigWriter XFigWriter.hh
      Writer for creating an XFig file from Ggeometric primitives */
  
  class XFigWriter {
    
  public:
    
    /** xfig style information */
    typedef struct {
    public:
      int depth;
      int lineStyle;
      int thickness;
      int penColor;
      int fillColor;
      int areaFill;
      float dashStyle;
      int jointStyle;
      int capStyle;
      int arcBoxRadius;
      int forwardArrow;
      int backwardArrow;
      int arrowStyle;
      int arrowFill;
      float arrowThickness;
      float arrowWidth;
      float arrowHeight;
    } Style;
    
    
    /** default constructor */
    inline XFigWriter()
      : ownOstream( false ), os( NULL ),
	dpiMult( XFIG_INTERNAL_DPI/XFIG_EFFECTIVE_DPI )
    {
      resetStyle();
    }
    
    /** open output stream */
    bool open( ostream &os );
    
    /** open output file */
    bool open( char *fileName );
    
    /** close stream */
    bool close();
    
    //
    // style access
    //
    
    /** resest style to default */
    inline void resetStyle()
    {
      style= defaultStyle;
    }
    
    /** read/write access to the style struct */
    inline Style &getStyle()
    {
      return style;
    }
    
    /** replace style with a new struct */
    inline void setStyle( Style &newStyle )
    {
      style= newStyle;
    }
    
    //
    // primitives...
    //
    
    /** a single line segment */
    inline void lineSegment( double p1x, double p1y, double p2x, double p2y )
    {
      drawLine( lround( dpiMult*p1x ), lround( dpiMult*p1y ),
		lround( dpiMult*p2x ), lround( dpiMult*p2y ) );
    }
    
    /** a single line segment (CoordinateVector interface) */
    inline void lineSegment( const CoordinateVector &p1,
			     const CoordinateVector &p2 )
    {
      errorCond( p1.vec.size()>= 2 && p2.vec.size()>= 2,
		 "1D output is not supported" );
      lineSegment( p1.vec[0], p1.vec[1], p2.vec[0], p2.vec[1] );
    }
    
    /** a single line segment (Vector interface) */
    inline void lineSegment( const Vector &p1, const Vector &p2 )
    {
      errorCond( p1.getSize()>= 2 && p2.getSize()>= 2,
		 "1D output is not supported" );
      lineSegment( p1[0], p1[1], p2[0], p2[1] );
    }
    
    
    /** an ellipse */
    inline void ellipse( double cx, double cy, double r1, double r2,
			 double angle= 0.0 )
    {
      drawEllipse( 1, lround( dpiMult*cx ), lround( dpiMult*cy ),
		   lround( dpiMult*r1 ), lround( dpiMult*r2 ), angle );
    }
    
    /** an ellipse (CoordinateVector interface) */
    inline void ellipse( const CoordinateVector &center,
			 const CoordinateVector &radii, double angle=0.0 )
    {
      errorCond( center.vec.size()>= 2 && radii.vec.size()>= 2,
		 "1D output is not supported" );
      ellipse( center.vec[0], center.vec[1], radii.vec[0], radii.vec[1], angle);
    }
    
    /** an ellipse (Vector interface) */
    inline void ellipse( const Vector &center, const Vector &radii,
			 double angle=0.0 )
    {
      errorCond( center.getSize()>= 2 && radii.getSize()>= 2,
		 "1D output is not supported" );
      ellipse( center[0], center[1], radii[0], radii[1], angle );
    }
    
    
    /** a circle */
    inline void circle( double cx, double cy, double radius )
    {
      drawEllipse( 3, lround( dpiMult*cx ), lround( dpiMult*cy ),
		   lround( dpiMult*radius ), lround( dpiMult*radius ), 0.0 );
    }
    
    /** a circle (CoordinateVector interface) */
    inline void circle( const CoordinateVector &center, double radius )
    {
      errorCond( center.vec.size()>= 2, "1D output is not supported" );
      circle( center.vec[0], center.vec[1], radius );
    }
    
    /** a circle (Vector interface) */
    inline void circle( const Vector &center, double radius )
    {
      errorCond( center.getSize()>= 2, "1D output is not supported" );
      circle( center[0], center[1], radius );
    }
    
    
    /** an image (from an external image file) */
    inline void image( const char *fileName, double width, double height,
		       double xOff= 0.0, double yOff= 0.0 )
    {
      drawImage( fileName, lround( dpiMult*width ), lround( dpiMult*height ),
		 lround( dpiMult*xOff ), lround( dpiMult*yOff ) );
    }
    
    /** an image (from an external image file, CoordinateVector interface) */
    inline void image( const char *fileName, const CoordinateVector &size,
		       const CoordinateVector &off )
    {
      errorCond( size.vec.size()>= 2 && off.vec.size()>= 2,
		 "1D output is not supported" );
      image( fileName, size.vec[0], size.vec[1], off.vec[0], off.vec[1] );
    }
    
    /** an image (from an external image file, CoordinateVector interface) */
    inline void image( const char *fileName, const Vector &size,
		       const Vector &off )
    {
      errorCond( size.getSize()>= 2 && off.getSize()>= 2,
		 "1D output is not supported" );
      image( fileName, size[0], size[1], off[0], off[1] );
    }
    
  protected:
    
    //
    // actual output methods for various xfig primitives
    //
    
    /** actually draw a single line */
    void drawLine( int p1x, int p1y, int p2x, int p2y );
    
    /** actually draw a single ellipse or circle */
    void drawEllipse( int type, int cx, int cy, int r1, int r2, double angle );
    
    /** actually draw image */
    void drawImage( const char *fileName, int width, int height,
		    int xOff, int yOff );
    
    //
    // fields
    //
    
    /** mutliplier to map from effective DPI to internal DPI */
    double dpiMult;
    
    /** output an XFig header */
    bool writeHeader();
    
    /** output stream */
    ostream *os;
    
    /** whether or not the ostream was created within this object */
    bool ownOstream;
    
    /** file-based ostream ceated within this object */
    ofstream myOstream;
    
    /** drawing style */
    Style style;
    
    /** default drawing style */
    static const Style defaultStyle;
    
  };


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_XFIGWRITER_H */

