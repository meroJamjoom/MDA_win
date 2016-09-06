// ==========================================================================
// $Id:$
// A whitepoint in x,y
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

#ifndef COLOR_WHITEPOINT_H
#define COLOR_WHITEPOINT_H

/*! \file  Whitepoint.hh
    \brief A whitepoint in x,y
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "ColorSpace.hh"

namespace MDA {

  using namespace std;
  
  /** \class Whitepoint Whitepoint.hh
      A whitepoint in XYZ */
  
  class Whitepoint {

  public:
    
    /** a whitepoint from the standards list (doubles as default constructor) */
    inline Whitepoint( const char *name= "D65", double scale= 100.0 )
      : XYZ( 3 )
    {
      setWhitepoint( name, scale );
    }
    
    /** constructor from a color in an arbitrary color space */
    inline Whitepoint( const Vector &color, const ColorSpace &space )
      : XYZ( 3 )
    {
      setWhitepoint( color, space );
    }
    
    /** constructor from x,y ccordinates */
    inline Whitepoint( double x, double y, double scale= 100.0 )
      : XYZ( 3 )
    {
      setWhitepoint( x, y, scale );
    }
    
    /** return XYZ coordinates of the whitepoint */
    inline const Vector &getXYZ() const
    {
      return XYZ;
    }
    
    /** set whitepoint to names strandard */
    void setWhitepoint( const char *name, double relativeScale= 100.0 );
    
    /** set whitepoint in x,y coordinates */
    inline void setWhitepoint( double x, double y, double relativeScale= 100.0 )
    {
      XYZ[1]= relativeScale;
      XYZ[0]= relativeScale * x / y;
      XYZ[2]= relativeScale * (1.0-x-y) / y;
    }
    
    /** set whitepoint from a color in an arbitrary color space */
    inline void setWhitepoint(  const Vector &color, const ColorSpace &space )
    {
      space.toXYZ( color, XYZ );
    }
    
    /** return a list of the names of known standard gammas
	(names are appended to end of existing list) */
    static void getStandardNames( list<const char *> &names );

  protected:
    
    /** the XYZ coordinates of the whitepoint (normalized to Y=100) */
    Vector XYZ;
  };
  
  
  /** write Gamma info to an ostream */
  ostream &operator<<( ostream &os, const Whitepoint &wp );
  
  /** read Gamma info from an istream */
  istream &operator>>( istream &is, Whitepoint &wp );
  

  /** \class WhitepointOption Gamma.hh
      commandline option for whitepoints
  */
  class WhitepointOption: public CommandlineOption {
    
  public:
    
    /** constructor */
    WhitepointOption( Whitepoint &_wp, const char *msg= "\twhitepoint\n",
		      const char *longText= "--whitepoint",
		      const char *shortText= "-wp" )
      : CommandlineOption( msg, longText, shortText ), wp( _wp )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the gamma object */
    Whitepoint &wp;
    
  };
  


} /* namespace */

#endif /* COLOR_WHITEPOINT_H */

