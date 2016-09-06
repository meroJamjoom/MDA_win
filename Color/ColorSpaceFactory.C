// ==========================================================================
// $Id: ColorSpaceFactory.C 279 2009-02-15 01:50:46Z heidrich $
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

#ifndef COLOR_COLORSPACEFACTORY_C
#define COLOR_COLORSPACEFACTORY_C

#include "MDA/Color/CIELABSpace.hh"
#include "MDA/Color/CIELUVSpace.hh"
#include "MDA/Color/YuvSpace.hh"
#include "MDA/Color/LinearTristimulusSpace.hh"
#include "MDA/Color/GammaCorrectedSpace.hh"
#include "ColorSpaceFactory.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** names of spaces directly supported by the factory */
const char *ColorSpaceFactory::directSupportNames[]= {
  "Lab",
  "Luv",
  "Yuv",
  ""
};


/** default constructor */

ColorSpaceFactory::ColorSpaceFactory()
  : spaceName( "XYZ" ),
    colorSpaceOpt( spaceName,
		   "\tselect a standard color space\n",
		   "--color-space", "-cs" )
{
  LinearTristimulusSpace::getStandardNames( linearTriSpaces );
  colorSpaceOpt.addOptions( linearTriSpaces );
  for( int i= 0 ; directSupportNames[i][0]!= '\0' ; i++ )
    directSupportOpts.push_back( directSupportNames[i] );
  colorSpaceOpt.addOptions( directSupportOpts );
}

/** register the relevant command line options */
void
ColorSpaceFactory::registerOptions( CommandlineParser &parser )
{
  toneCurveFactory.registerOptions( parser );
  parser.registerOption( &colorSpaceOpt );
}

/** create a new ColorSpace object using the current parameters */
ColorSpace *
ColorSpaceFactory::makeColorSpace()
{
  if( !strcasecmp( spaceName, directSupportNames[0] ) )
    return new CIELABSpace();    // Lab
  else if( !strcasecmp( spaceName, directSupportNames[1] ) )
    return new CIELUVSpace();    // Luv
  else if( !strcasecmp( spaceName, directSupportNames[2] ) )
    return new YuvSpace();       // Yuv
  else
  {
    // linear space, possibly with gamma or other tone curve
    LinearTristimulusSpace *linSpace= new LinearTristimulusSpace( spaceName );
    if( toneCurveFactory.isLinear() )
      return linSpace;
    return new GammaCorrectedSpace( linSpace,
				    toneCurveFactory.makeToneCurve() );
  }
  
  // never reached
  return NULL;
}



} /* namespace */

#endif /* COLOR_COLORSPACEFACTORY_C */

