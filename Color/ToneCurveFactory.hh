// ==========================================================================
// $Id: ToneCurveFactory.hh 279 2009-02-15 01:50:46Z heidrich $
// factory for creating various tone curves
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

#ifndef COLOR_TONECURVEFACTORY_H
#define COLOR_TONECURVEFACTORY_H

/*! \file  ToneCurveFactory.hh
    \brief factory for creating various tone curves
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CommandlineParser.hh"
#include "Gamma.hh"

namespace MDA {

  /** \class ToneCurveFactory ToneCurveFactory.hh
      factory for creating various tone curves */
  
  class ToneCurveFactory {

  public:

    /** default constructor */
    ToneCurveFactory();
    
    /** register the relevant command line options */
    void registerOptions( CommandlineParser &parser );
    
    /** whether or not the mapping is linear */
    bool isLinear();
    
    /** create a new tone curve using the current parameters */
    ToneCurve *makeToneCurve();
    
    /** set parameters for various types of tone curves */
    inline void setGammaParams( double _gamma= 1.0, double _gain= 1.0,
				double _bias= 0.0, double _thres= 0.0,
				double _slope= 1.0, bool _clamp= true )
    {
      gamma= _gamma;	gain= _gain;	bias=_bias;
      thres= _thres;	slope=_slope;	clamp= _clamp;
      gammaStandard= NULL;
    }
    
  protected:
    
    /** list of supported gamma standards */
    list<const char *> gammaStandards;
    
    // Gamma curve parameters
    
    /** gamma */
    double gamma;

    /** gamma option */
    DoubleOption gammaOption;
    
    /** gamma linear gain */
    double gain;
    
    /** gamma linear gain option */    
    DoubleOption gainOption;
    
    /** gamma linear bias */
    double bias;
    
    /** gamma linear bias option */
    DoubleOption biasOption;
    
    /** gamma linear threshold */
    double thres;
    
    /** gamma linear threshold option */
    DoubleOption thresholdOption;
    
    /** gamma linear slope */
    double slope;
    
    /** gamma linear slope option */
    DoubleOption slopeOption;
    
    /** gamma clamp to positive values */
    bool clamp;
    
    /** gamma clamp to positive values option */
    BoolOption clampOption;
    
    /** gamma standard */
    const char *gammaStandard;
    
    /** gamma standard option */    
    StringSelectorOption standardOpt;
  };


} /* namespace */

#endif /* COLOR_TONECURVEFACTORY_H */

