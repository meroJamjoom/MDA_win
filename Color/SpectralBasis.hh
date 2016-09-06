// ==========================================================================
// $Id:$
// colors as a linear combination of basis spectra
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

#ifndef COLOR_SPECTRALBASIS_H
#define COLOR_SPECTRALBASIS_H

/*! \file  SpectralBasis.hh
    \brief colors as a linear combination of basis spectra
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>

#include "MDA/LinearAlgebra/LinAlg.hh"
#include "MDA/Array/Array.hh"
#include "Spectrum.hh"


namespace MDA {

  using namespace std;
  
  /** \class SpectralBasis SpectralBasis.hh
      colors as a linear combination of basis spectra */
  
  class SpectralBasis {

  public:

    /** default constructor */
    inline SpectralBasis()
      : dualInitialized( false ), dualBasis()
    {}
    
    /** constructor from Array channels (float version)
     *  \param a the array
     *  \param channels list of channels
     *  \param aRange the range of wavelengths represented by the data
     *  values in the channels. The samples in the array are expected
     *  to be spread uniformly across this interval
     *  \param range the final range of the spectrum we create
     */
    inline SpectralBasis( Array<double> &array, const ChannelList &channels,
			  const IntRange &aRange, const IntRange &range )
      : dualInitialized( false ), dualBasis()
    {
      unsigned numBasis= channels.vec.size();
      basis.reserve( numBasis );
      for( unsigned i= 0 ; i< numBasis ; i++ )
	basis.push_back( new Spectrum( array, channels.vec[i], aRange, range ));
    }
    
    /** constructor from Array channels (double version)
     *  \param a the array
     *  \param channels list of channels
     *  \param aRange the range of wavelengths represented by the data
     *  values in the channels. The samples in the array are expected
     *  to be spread uniformly across this interval
     *  \param range the final range of the spectrum we create
     */
    inline SpectralBasis( Array<float> &array, const ChannelList &channels,
			  const IntRange &aRange, const IntRange &range )
      : dualInitialized( false ), dualBasis()
    {
      unsigned numBasis= channels.vec.size();
      basis.reserve( numBasis );
      for( unsigned i= 0 ; i< numBasis ; i++ )
	basis.push_back( new Spectrum( array, channels.vec[i], aRange, range ));
    }
    
    /** destructor */
    inline ~SpectralBasis()
    {
      for( unsigned i= 0 ; i< basis.size() ; i++ )
	delete basis[i];
    }
    
    /** XYZ coordinates of a linear combination of the basis */
    virtual void toXYZ( const Vector &coeff, Vector &XYZ );
    
    /** XYZ coordinates for dual basis (this method can be used to get
	tristimulus color from a multiband measurement using the
	spectral basis as color filter) */
    virtual void dualToXYZ( const Vector &dualCoeff, Vector &XYZ );
    
  protected:
    
    /** basis spectra */
    vector<Spectrum *> basis;
    
    /** whether the dual has beeen computed already */
    bool dualInitialized;
    
    /** coefficients of the dual basis, represented in matrix form */
    Matrix dualBasis;
  };


} /* namespace */

#endif /* COLOR_SPECTRALBASIS_H */

