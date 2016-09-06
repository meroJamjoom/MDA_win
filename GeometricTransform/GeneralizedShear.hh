// ==========================================================================
// $Id:$
// A shear, extended with a 1D scale/translation
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

#ifndef GEOMETRICTRANSFORM_GENERALIZEDSHEAR_H
#define GEOMETRICTRANSFORM_GENERALIZEDSHEAR_H

/*! \file  GeneralizedShear.hh
    \brief A shear, extended with a 1D scale/translation
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "GeometricTransformation.hh"


namespace MDA {

  /** \class GeneralizedShear GeneralizedShear.hh
      a shear, extended with a 1D scale/translation */
  template<class T>
  class GeneralizedShear: public GeometricTransformation<T> {

  public:

    /** default constructor */
    inline GeneralizedShear( Resampler<T> *_resampler= NULL )
      : GeometricTransformation<T>( _resampler ), axis( -1 )
    {}

    /** set shear parameters */
    void setShear( unsigned _axis, vector<double> &coeff )
    {
      coefficients= coeff;
      axis= _axis;
    }
    
    /** apply the transform to a single channel in one array */
    virtual bool apply( Array<T> &srcArray, int srcChannel,
                        Array<T> &dstArray, int dstChannel );
    
  protected:
    
    /** vector of coefficients along each axis (+1 for translation) */
    vector<double> coefficients;
    
    /** the axis to apply the shear to */
    unsigned axis;
    
  };


} /* namespace */

#endif /* GEOMETRICTRANSFORM_GENERALIZEDSHEAR_H */

