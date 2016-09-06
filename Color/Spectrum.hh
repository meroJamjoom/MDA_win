// ==========================================================================
// $Id:$
// a spectral representation of light
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

#ifndef COLOR_SPECTRUM_H
#define COLOR_SPECTRUM_H

/*! \file  Spectrum.hh
    \brief a spectral representation of light
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/Range.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "MDA/Array/Array.hh"

namespace MDA {

  /** \class Spectrum Spectrum.hh
      a spectral representation of light */
  
  class Spectrum {

  public:

    /** constructor (default with empty range) */
    inline Spectrum( IntRange _range= IntRange( 0, -1 ) )
      : range( _range )
    {
      int numElem= range.val.second - range.val.first + 1;
      if( numElem> 0 )
      {	
	data= new double[numElem];
	bzero( data, numElem*sizeof( double ) );
      }
      else
	data= NULL;
    }
    
    /** constructor from named spectrum */
    Spectrum( const char *name )
    {
      data= NULL;
      range.val.first= 0;
      range.val.second= -1;
      set( name );
    }
    
    /** constructor from an array channel (float version)
     *  \param a the array
     *  \param channel channel index
     *  \param aRange the range of wavelengths represented by the data
     *  values in the channel. The samples in the array are expected
     *  to be spread uniformly across this interval
     *  \param range the final range of the spectrum we create
     */
    Spectrum( Array<float> &a, unsigned channel,
	      const IntRange &aRange, const IntRange &range );
    
    /** constructor from an array channel (double version)
     *  \param a the array
     *  \param channel channel index
     *  \param aRange the range of wavelengths represented by the data
     *  values in the channel. The samples in the array are expected
     *  to be spread uniformly across this interval
     *  \param range the final range of the spectrum we create
     */
    Spectrum( Array<double> &a, unsigned channel,
	      const IntRange &aRange, const IntRange &range );
    
    /** blackbody spectrum (in W/m^2)
     *  \param T: the temperature in Kelvin
     *  \param r: spectral range of interest (in nm)
     */
    Spectrum( double T, IntRange r= IntRange( 360, 830 ) );
    
    /** copy constructor */
    inline Spectrum( const Spectrum &orig )
      : range( orig.range )
    {
      int numElem= range.val.second - range.val.first + 1;
      if( numElem> 0 )
      {	
	data= new double[numElem];
	memcpy( data, orig.data, numElem*sizeof( double ) );      
      }
      else
	data= NULL;
    }
    
    /** destructor */
    inline ~Spectrum()
    {
      if( data!= NULL )
	delete [] data;
    }
    
    /** set the whole spectrum to a known value */
    void set( const char *name );
    
    /** get spectral value for one wavelength */
    inline double get( unsigned lambda ) const
    {
      if( lambda< range.val.first || lambda> range.val.second )
	return 0.0;
      else
	return data[lambda-range.val.first];
    }
    
    /** set spectral value for one wavelength
     *  \returns true on success, false if wavelength out of range
     */
    inline bool set( unsigned lambda, double value )
    {
      if( lambda< range.val.first || lambda> range.val.second )
	return false;
      data[lambda-range.val.first]= value;
      return true;
    }
    
    /** return CIE XYZ values
     *  \param XYZ: result vector
     *  \param use1964: whether to use the 10 degree (1964) matching
     *  functions or the 2 degree (1931) ones
     */
    void toXYZ( Vector &XYZ, bool use1964= true );
    
    /** spectral integral */
    inline double integral()
    {
      double result= 0.0;
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	result+= data[i];
      return result;
    }
    
    /** multiply with scalar */
    inline Spectrum &operator*=( double scalar)
    {
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	data[i]*= scalar;
      return *this;
    }
    
    /** divide by scalar */
    inline Spectrum &operator/=( double scalar)
    {
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	data[i]/= scalar;
      return *this;
    }
    
    /** add scalar */
    inline Spectrum &operator+=( double scalar)
    {
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	data[i]+= scalar;
      return *this;
    }
    
    /** subtract scalar */
    inline Spectrum &operator-=( double scalar)
    {
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	data[i]-= scalar;
      return *this;
    }
    
    /** clamp negative amplitudes to 0 */
    inline Spectrum &clampNegative()
    {
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	if( data[i]< 0.0 )
	  data[i]= 0.0;
      return *this;
    }
    
    /** normalize such that the power (i.e integral over wavelengths) is 1 */
    inline Spectrum &normalizePower()
    {
      return *this/= integral();
    }
    
    /** normalize such that the maximum value is 1 */
    inline Spectrum &normalizeAmplitude()
    {
      // assuming there is at least one positive value
      double max= -1.0;
      int numElem= range.val.second - range.val.first + 1;
      for( int i= 0 ; i< numElem ; i++ )
	max= data[i]> max ? data[i] : max;
      for( int i= 0 ; i< numElem ; i++ )
	data[i]/= max;
      return *this;
    }
        
    /** multiply another spectrum onto this one */
    inline Spectrum &operator*=( const Spectrum &other )
    {
      // if the ranges don't match, use the more sophisticated
      // multiplication method, otherwise just multiply in-space
      if( range.val!= other.range.val )
      {
	// save current data temporarily
	double *oldData= data;
	data= NULL;
	
	multiply( oldData, range, other.data, other.range, data, range );
	// clean up
	delete [] oldData;
      }
      else
      {
	int numElem= range.val.second - range.val.first + 1;
	for( int i= 0 ; i< numElem ; i++ )
	  data[i]*= other.data[i];
      }
      
      return *this;
    }
    
    /** add another spectrum onto this one */
    inline Spectrum &operator+=( const Spectrum &other )
    {
      // if the ranges don't match, use the more sophisticated
      // multiplication method, otherwise just multiply in-space
      if( range.val!= other.range.val )
      {
	// save current data temporarily
	double *oldData= data;
	data= NULL;
	add( oldData, range, other.data, other.range, data, range );
	
	// clean up old array
	delete [] oldData;
      }
      else
      {
	int numElem= range.val.second - range.val.first + 1;
	for( int i= 0 ; i< numElem ; i++ )
	  data[i]+= other.data[i];
      }
      
      return *this;
    }
    
    /** *this = *this + scalar * other */
    inline Spectrum &addScalarTimesSpectrum( double scalar,
					     const Spectrum &other )
    {
      // if the ranges don't match, use the more sophisticated
      // multiplication method, otherwise just multiply in-space
      if( range.val!= other.range.val )
      {
	// save current data temporarily
	double *oldData= data;
	data= NULL;
	add( oldData, range, other.data, other.range, data, range, scalar );
	
	// clean up old array
	delete [] oldData;
      }
      else
      {
	int numElem= range.val.second - range.val.first + 1;
	for( int i= 0 ; i< numElem ; i++ )
	  data[i]+= scalar * other.data[i];
      }
      
      return *this;
    }
    
    /** product of two spectra */
    inline Spectrum operator*( const Spectrum &other )
    {
      Spectrum newSpectrum( IntRange( 0, -1 ) ); // initially empty spectrum

      multiply( data, range, other.data, other.range,
		newSpectrum.data, newSpectrum.range );
      
      return newSpectrum;
    }

    /** sum of two spectra */
    inline Spectrum operator+( const Spectrum &other )
    {
      Spectrum newSpectrum( IntRange( 0, -1 ) ); // initially empty spectrum
      
      add( data, range, other.data, other.range,
	   newSpectrum.data, newSpectrum.range );
      
      return newSpectrum;
    }
    
  protected:
    
    /** multiply two spectral arrays into an existing third array */
    void multiply( const double *s1, const IntRange &r1,
		   const double *s2, const IntRange &r2,
		   double *&d, IntRange &rd );
    
    /** add two spectral arrays into an existing third array (with an
	optional scalar multiplier for the second array) */
    void add( const double *s1, const IntRange &r1,
	      const double *s2, const IntRange &r2,
	      double *&d, IntRange &rd,
	      double s2factor= 1.0 );
    
    /** value range */
    IntRange range;
    
    /** data vector */
    double *data;
    
  };


} /* namespace */

#endif /* COLOR_SPECTRUM_H */

