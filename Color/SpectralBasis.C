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

#ifndef COLOR_SPECTRALBASIS_C
#define COLOR_SPECTRALBASIS_C

#include "MDA/Base/Errors.hh"
#include "MDA/LinearAlgebra/LinAlg.hh"
#include "SpectralBasis.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

    
/** XYZ coordinates of a linerar combination of the basis */
void
SpectralBasis::toXYZ( const Vector &coeff, Vector &XYZ )
{
  unsigned dimensions= coeff.getSize();
  errorCond( dimensions== basis.size(),
	     "  #coefficents does not match #basis functions" );
  
  if( dimensions== 0 )
    XYZ.zero();
  else
  {
    Spectrum accumulator( *basis[0] );
    accumulator*= coeff[0];
    for( unsigned i= 1 ; i< dimensions ; i++ )
      accumulator.addScalarTimesSpectrum( coeff[i], *basis[i] );
    accumulator.toXYZ( XYZ );
  }
}
  
/** XYZ coordinates for dual basis (this method can be used to get
    tristimulus color from a multiband measurement using the
    spectral basis as color filter) */
void
SpectralBasis::dualToXYZ( const Vector &dualCoeff, Vector &XYZ )
{
  unsigned i, j;
  
  unsigned dimensions= dualCoeff.getSize();
  errorCond( dimensions== basis.size(),
	     "  #coefficents does not match #basis functions" );
  
  if( dimensions== 0 )
  {
    XYZ.zero();
    return;
  }
  // create dual basis matrix
  if( !dualInitialized )
  {
    // matrix of dot products of all basis spectra
    Matrix dps( dimensions, dimensions );
    for( i= 0 ; i< dimensions ; i++ )
      for( j= i ; j< dimensions ; j++ )
      {
	Spectrum dp( *basis[i] );
	dp*= *basis[j];
	dps[i][j]= dps[j][i]= dp.integral();
      }
    
    // dualBasis is the inverse of dps
    dualBasis= inverse( dps );
  }
  
  // convert dual coefficients to primary coefficients, then call
  // function for XYZ conversion from primary coefficients
  Vector coeff( dualCoeff.getSize() );
  dualBasis.rightMultiply( dualCoeff, coeff );
  toXYZ( coeff, XYZ );
}


} /* namespace */

#endif /* COLOR_SPECTRALBASIS_C */

