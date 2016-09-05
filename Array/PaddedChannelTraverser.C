// ==========================================================================
// $Id: PaddedChannelTraverser.C 384 2009-09-25 19:48:11Z heidrich $
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

#ifndef ARRAY_PADDEDCHANNELTRAVERSER_C
#define ARRAY_PADDEDCHANNELTRAVERSER_C

#include "PaddedChannelTraverser.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** constructor */

PaddedChannelTraverser::PaddedChannelTraverser( CoordinateVector &_dim,
						unsigned long _padding )
  : dim( _dim ), iter( _dim ),
    padding( _padding ), paddedDim( pad( _dim, _padding ) )
{
  scanlinePixels= _dim.vec[0];
  numScanlines= scanlines( _dim );
  numPixels= numScanlines * scanlinePixels;
}



} /* namespace */

#endif /* ARRAY_PADDEDCHANNELTRAVERSER_C */

