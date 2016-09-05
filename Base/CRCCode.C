// ==========================================================================
// $Id: CRCCode.C 661 2010-03-22 01:19:14Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef BASE_CRCCODE_C
#define BASE_CRCCODE_C

#include "CRCCode.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** encode a data value */
unsigned long
CRCCode::encode( unsigned long data ) const
{
  unsigned long i;
  
  // bounds check...
  if( data>= (1<<dBits) )
  {
    warning( "data exceeds allocated bit depth" );
    data&= (1<<dBits)-1;
  }
  
  // polygon division; crc code is remainder
  unsigned long crc= data << cBits;
  unsigned long div= poly << dBits;
  while( div>= poly )
  {
    if( (crc^div)< crc )
      crc^= div;
    div>>= 1;
  }
  
  // pack crc code in cBits LSBs
  return (data << cBits) | crc;
}



} /* namespace */

#endif /* BASE_CRCCODE_C */

