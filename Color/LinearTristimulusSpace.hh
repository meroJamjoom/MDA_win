// ==========================================================================
// $Id:$
// a tristimulus space that can be generated from XYZ through a linear transform
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

#ifndef COLOR_LINEARTRISTIMULUSSPACE_H
#define COLOR_LINEARTRISTIMULUSSPACE_H

/*! \file  LinearTristimulusSpace.hh
    \brief a tristimulus space that can be generated from XYZ through a linear transform
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>

#include "MDA/LinearAlgebra/LinAlg.hh"
#include "ColorSpace.hh"

namespace MDA {

  using namespace std;
  
  /** \class LinearTristimulusSpace LinearTristimulusSpace.hh
      a tristimulus space that can be generated from XYZ through a
      linear transform */
  
  class LinearTristimulusSpace: public ColorSpace {

  public:

    /** constructor for named space - doubles as default constructor */
    LinearTristimulusSpace( const char *name= "sRGB" );
    
    /** constructor form custom matrices (forward and inverse) */
    inline LinearTristimulusSpace( Matrix &fromXYZ,
				   Matrix &toXYZ )
      : toXYZMatrix( toXYZ ), fromXYZMatrix( fromXYZ )
    {}
    
    /** copy constuctor */
    inline LinearTristimulusSpace( const LinearTristimulusSpace &other )
      : toXYZMatrix( other.toXYZMatrix ), fromXYZMatrix( other.fromXYZMatrix )
    {}
    
    /** report dimensionality of space */
    inline virtual unsigned getDimension() const
    {
      return 3;
    }
    
    /** convert from this space to XYZ */
    virtual void toXYZ( const Vector &tristimulus, Vector &XYZ ) const;
    
    /** convert from XYZ to this space */
    virtual void fromXYZ( const Vector &XYZ, Vector &tristimulus ) const;
    
    /** return a list of the names of known standard spaces
	(the strings are apppended to the end of the existing list)  */
    static void getStandardNames( list<const char *> &names );
    
  protected:

    /** matrix for conversion to XYZ */
    Matrix toXYZMatrix;
    
    /** matrix for conversion from XYZ */
    Matrix fromXYZMatrix;

  };


} /* namespace */

#endif /* COLOR_LINEARTRISTIMULUSSPACE_H */

