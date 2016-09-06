// ==========================================================================
// $Id: ResamplerThreading.C 712 2010-04-23 09:16:34Z martinle $
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

#ifndef RESAMPLING_RESAMPLERTHREADING_C
#define RESAMPLING_RESAMPLERTHREADING_C

#include "ResamplerThreading.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc oes not yield
  // problems with other packages!
  using namespace std;

/** resample a single line */
template <class T>
void
LinearResamplerJob<T>::execute( int threadID )
{
  resampler->resampleLinear( inLine, inCount, inStride,
			     outLine, outCount, outStride, a, b );
}


/** resample a single line */
template <class T>
void
RationalResamplerJob<T>::execute( int threadID )
{
  LinearResamplerJob<T>::resampler->
    resampleRational( LinearResamplerJob<T>::inLine,
		      LinearResamplerJob<T>::inCount,
		      LinearResamplerJob<T>::inStride,
		      LinearResamplerJob<T>::outLine,
		      LinearResamplerJob<T>::outCount,
		      LinearResamplerJob<T>::outStride,
		      LinearResamplerJob<T>::a,
		      LinearResamplerJob<T>::b,
		      c, d );
}


/** resample a single line */
template <class T>
void
IrregularResamplerJob<T>::execute( int threadID )
{
  LinearResamplerJob<T>::resampler->
    resampleIrregular( LinearResamplerJob<T>::inLine,
		       LinearResamplerJob<T>::inCount,
		       LinearResamplerJob<T>::inStride,
		       LinearResamplerJob<T>::outLine,
		       LinearResamplerJob<T>::outCount,
		       LinearResamplerJob<T>::outStride,
		       samples, sampleStride );
}


// explicit template instation code
template class LinearResamplerJob<float>;
template class LinearResamplerJob<double>;
template class LinearResamplerJob<unsigned char>;
template class RationalResamplerJob<float>;
template class RationalResamplerJob<double>;
template class RationalResamplerJob<unsigned char>;
template class IrregularResamplerJob<float>;
template class IrregularResamplerJob<double>;
template class IrregularResamplerJob<unsigned char>;



} /* namespace */

#endif /* RESAMPLING_RESAMPLERTHREADING_C */

