// ==========================================================================
// $Id:$
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

#ifndef COLOR_GAMMA_C
#define COLOR_GAMMA_C

#include <string.h>
#include <math.h>
#include <iostream>

#include "Gamma.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

//
// local data for this file
//

class NamedGamma {
public:
  const char *name;
  double gamma;
  double gain;
  double bias;
  double threshold;
  double slope;
  bool clamp;
};

// names of standard gammas. first entry is the default
NamedGamma standardGammas[] = {
  { "2.2",		2.2, 1.000,  0.000, 0.00000,  1.00, true }, // default
  { "sRGB",		2.4, 1.055, -0.055, 0.00304, 12.92, true },
  { "Rec709",      1.0/0.45, 1.099, -0.099, 0.01800,  4.50, true },
  { "xvYCC",       1.0/0.45, 1.099, -0.099, 0.01800,  4.50, false },
  { "Apple",		1.8, 1.000,  0.000, 0.00000,  1.00, true },
  { "" } // end of list
};

//
// Gamma methods
//


/** set gamma paramters explictily */
void
Gamma::setGamma( double _gamma, double _gain, double _bias,
		 double _threshold, double _slope, bool _clamp )
{
  standard[0]= '\0';
  gamma= _gamma;
  gain= _gain;
  bias= _bias;
  threshold= _threshold;
  slope= _slope;
  clamp= _clamp;
}

/** a standard gamma curve */
void
Gamma::setGamma( const char *name )
{
  unsigned i= 0;
  // find standard with matching name, or last default entry
  while( standardGammas[i].name[0]!= '\0' )
    if( !strcmp( standardGammas[i].name, name ) )
      break;
    else
      i++;
  if( standardGammas[i].name[0]== '\0' )
  {
    i= 0;
    cerr << "Gamma standard \"" << name << "\" not found - using \""
	 << standardGammas[i].name << '\"' << endl;
  }
  
  // copy specifics of the standard
  strcpy( standard, standardGammas[i].name );
  gamma= standardGammas[i].gamma;
  gain= standardGammas[i].gain;
  bias= standardGammas[i].bias;
  threshold= standardGammas[i].threshold;
  slope= standardGammas[i].slope;
  clamp= standardGammas[i].clamp;
}
    
/** conversion from linear to non-linear space */
double
Gamma::fromLinear( double linearValue )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( linearValue< 0.0 )
    if( clamp )
      linearValue= 0.0;
    else
    {
      inversion= -1.0;
      linearValue= -linearValue;
    }
  // clamp against 1
  if( linearValue> 1.0 && clamp )
    linearValue= 1.0;
  
  if( linearValue<= threshold )
    return inversion * linearValue * slope;
  else
    return inversion * (gain * pow( linearValue, 1.0/gamma ) + bias);
}
    
/** conversion from non-linear to linear space */
double
Gamma::toLinear( double nonlinearValue )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( nonlinearValue< 0.0 )
    if( clamp )
      nonlinearValue= 0.0;
    else
    {
      inversion= -1.0;
      nonlinearValue= -nonlinearValue;
    }
  // clamp against 1
  if( nonlinearValue> 1.0 && clamp )
    nonlinearValue= 1.0;
  
  if( nonlinearValue<= slope*threshold )
    return inversion * nonlinearValue / slope;
  else
    return inversion * pow( (nonlinearValue - bias) / gain, gamma );
}


/** return a list of the names of known standard gammas */
void
Gamma::getStandardNames( list<const char *> &names )
{
  unsigned i= 0;
  while( standardGammas[i].name[0]!= '\0' )
    names.push_back( standardGammas[i++].name );
}


} /* namespace */

#endif /* COLOR_GAMMA_C */

