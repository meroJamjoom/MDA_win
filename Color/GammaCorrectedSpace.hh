// ==========================================================================
// $Id:$
// gamma corrected version of some other space
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

#ifndef COLOR_GAMMACORRECTEDSPACE_H
#define COLOR_GAMMACORRECTEDSPACE_H

/*! \file  GammaCorrectedSpace.hh
    \brief gamma corrected version of some other space
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "LinearTristimulusSpace.hh"
#include "Gamma.hh"

namespace MDA {

  /** \class GammaCorrectedSpace GammaCorrectedSpace.hh
      gamma corrected version of some other space */
  
  class GammaCorrectedSpace: public ColorSpace {

  public:

    /** constructor from a known standard */
    inline GammaCorrectedSpace( char *name )
    {
      space= new LinearTristimulusSpace( name );
      gamma= new Gamma( name );
    }
    
    /** constructor from a standard space & gamma */
    inline GammaCorrectedSpace( char *spaceName, char *gammaName )
    {
      space= new LinearTristimulusSpace( spaceName );
      gamma= new Gamma( gammaName );
    }
    
    /** constructor from separately generated space & gamma
     *  (the GammaCorrectedSpace takes ownership of these objects)
     */
    inline GammaCorrectedSpace( ColorSpace *s, ToneCurve *g )
      : space( s ), gamma( g )
    {}
    
    /** destructor */
    virtual ~GammaCorrectedSpace()
    {
      delete space;
      delete gamma;
    }
    
    /** report dimensionality of space */
    inline virtual unsigned getDimension() const
    {
      return space->getDimension();
    }
    
    /** convert from this space to XYZ */
    virtual void toXYZ( const Vector &tristimulus, Vector &XYZ ) const;
    
    /** convert from XYZ to this space */
    virtual void fromXYZ( const Vector &XYZ, Vector &tristimulus ) const;

  protected:

    /** the (typically linear) reference space */
    ColorSpace *space;
    
    /** the gamma curve */
    ToneCurve *gamma;
 
  };


} /* namespace */

#endif /* COLOR_GAMMACORRECTEDSPACE_H */

