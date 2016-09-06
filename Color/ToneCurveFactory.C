// ==========================================================================
// $Id: ToneCurveFactory.C 308 2009-04-05 10:50:49Z heidrich $
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

#ifndef COLOR_TONECURVEFACTORY_C
#define COLOR_TONECURVEFACTORY_C

#include "ToneCurveFactory.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** default constructor */
ToneCurveFactory::ToneCurveFactory()
  : gammaOption( gamma, "\tgamma value\n", "--gamma", NULL ),
    thresholdOption( thres, "\tthreshold for linear vs. gamma curve\n",
		     "--threshold", NULL ),
    slopeOption( slope, "\tslope of the linar portion of the curve\n",
		 "--slope", NULL ),
    gainOption( gain, "\tgain for gamma curve\n", "--gain", NULL ),
    biasOption( bias, "\tbias for gamma curve\n", "--bias", NULL ),
    clampOption( clamp, "\tclamp values to 0...1?\n"
		 "\t(some standards allow for values outside this range)\n",
		 "--clamp", NULL, "--no-clamp", NULL ),
    standardOpt( gammaStandard,
		 "\tselect parameters according to image/video standard\n"
		 "\t(this will override specific gamma curve parameters)\n",
		 "--standard-gamma", "-sg" )
{
  Gamma::getStandardNames( gammaStandards );
  setGammaParams();
  standardOpt.addOptions( gammaStandards );
}


/** register the relevant command line options */
void
ToneCurveFactory::registerOptions( CommandlineParser &parser )
{
  //
  // register gamma curve options
  //
  parser.registerOption( &gammaOption );
  parser.registerOption( &thresholdOption );
  parser.registerOption( &slopeOption );
  parser.registerOption( &gainOption );
  parser.registerOption( &biasOption );
  parser.registerOption( &clampOption );

  parser.registerOption( &standardOpt );
}

/** whether or not the mapping is linear */
bool
ToneCurveFactory::isLinear()
{
  return gammaStandard== NULL && gamma== 1.0 && gain== 1.0 &&
    bias== 0.0 && thres== 0.0 && slope== 1.0 && clamp== true;
}

/** create a new tone curve using the current parameters */
ToneCurve *
ToneCurveFactory::makeToneCurve()
{
  // providing a standard curve?
  if( gammaStandard!= NULL )
    return new Gamma( gammaStandard );
  
  // default: create a user-specified gamma curve
  return new Gamma( gamma, gain, bias, thres, slope, clamp );
}


} /* namespace */

#endif /* COLOR_TONECURVEFACTORY_C */

