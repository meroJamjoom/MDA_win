// ==========================================================================
// $Id: ColorSpaceFactory.hh 279 2009-02-15 01:50:46Z heidrich $
// factory for creating various color spaces
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef COLOR_COLORSPACEFACTORY_H
#define COLOR_COLORSPACEFACTORY_H

/*! \file  ColorSpaceFactory.hh
    \brief factory for creating various color spaces
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CommandlineParser.hh"
#include "ColorSpace.hh"
#include "ToneCurveFactory.hh"

namespace MDA {

  /** \class ColorSpaceFactory ColorSpaceFactory.hh
      factory for creating various color spaces */
  
  class ColorSpaceFactory {

  public:

    /** default constructor */
    ColorSpaceFactory();

    /** register the relevant command line options */
    void registerOptions( CommandlineParser &parser );
    
    /** create a new ColorSpace object using the current parameters */
    ColorSpace *makeColorSpace();
    
  protected:
    
    /** list of standard linear tristimulus spaces */
    list<const char *> linearTriSpaces;
    
    /** list of spaces directly supported by the factory */
    list<const char *> directSupportOpts;
    
    /** names of spaces directly supported by the factory */
    static const char *directSupportNames[];
    
    /** name of the color space */
    const char *spaceName;
    
    /** a factory for tone curves */
    ToneCurveFactory toneCurveFactory;
    
    /** option for selecting a named standard colro space */
    StringSelectorOption colorSpaceOpt;

  };


} /* namespace */

#endif /* COLOR_COLORSPACEFACTORY_H */

