// ==========================================================================
// $Id:$
// simple gamma tone curves according to various standards
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

#ifndef COLOR_GAMMA_H
#define COLOR_GAMMA_H

/*! \file  Gamma.hh
    \brief simple gamma tone curves according to various standards
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>

#include "MDA/Base/CommandlineParser.hh"
#include "ToneCurve.hh"

namespace MDA {

  using namespace std;
  
  /** \class Gamma Gamma.hh
      simple gamma tone curves according to various standards
      
      The specific formula used is:
      
      \f[I_{out}=
      \left\{\begin{array}{ll}
      slope\cdot I_{in}&;\mathrm{if}\ I_{in}<threshold\\
      gain\cdot I_{in}^{1/gamma}+bias&;\mathrm{else}
      \end{array}\right.
      \f]
  */
  class Gamma: public ToneCurve {

  public:

    /** default constructor */
    inline Gamma( double gamma= 1.0, double gain= 1.0, double bias= 0.0,
		  double threshold= 0.0, double slope= 1.0, bool clamp= true )
    {
      setGamma( gamma, gain, bias, threshold, slope, clamp );
    }
    
    /** A standard gamma curve */
    inline Gamma( const char *name )
    {
      setGamma( name );
    }
    
    /** set gamma paramters explicitly */
    void setGamma( double gamma= 1.0, double gain= 1.0, double bias= 0.0,
		   double threshold= 0.0, double slope= 1.0, bool clamp= true );
    
    /** set bgamma parameters from named standard */
    void setGamma( const char *name );
    
    /** conversion from linear to non-linear space */
    virtual double fromLinear( double linearVal );
    
    /** conversion from non-linear to linear space */
    virtual double toLinear( double nonlinearVal );
    
    /** return a list of the names of known standard gammas
	(the strings are apppended to the end of the existing list) */
    static void getStandardNames( list<const char *> &names );

  protected:
    
    /** ostream operator can read components of the gamma class */
    friend ostream &operator<<( ostream &os, const Gamma &gamma );
    
    /** name of gamma standard, if set */
    char standard[128];
    
    /** gamma exponent for the non-linear part of the curve */
    double gamma;
    
    /** gain for non-linear portion of the curve */
    double gain;
    
    /** bias for non-linear part of curve */
    double bias;
    
    /** threshold between linear and non-linear part */
    double threshold;
    
    /** slope of linear part of curve */
    double slope;
    
    /** whether or not to clamp to the range 0..1 */
    bool clamp;
  };

  
} /* namespace */

#endif /* COLOR_GAMMA_H */

