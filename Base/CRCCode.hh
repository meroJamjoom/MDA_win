// ==========================================================================
// $Id: CRCCode.hh 661 2010-03-22 01:19:14Z heidrich $
// Cyclic Redundancy Check codes
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

#ifndef BASE_CRCCODE_H
#define BASE_CRCCODE_H

/*! \file  CRCCode.hh
    \brief Cyclic Redundancy Check codes
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Errors.hh"

namespace MDA {

  /** \class CRCCode CRCCode.hh
      Cyclic Redundancy Check codes */
  
  class CRCCode {

  public:

    /** default constructor */
    CRCCode( unsigned dataBits, unsigned crcBits, unsigned long generatorPoly )
      : dBits( dataBits ), cBits( crcBits ),
	poly( generatorPoly | (1 << crcBits) )
    {
      errorCond( dBits+cBits<= sizeof(unsigned long)*8,
		 "  sum of data and CRC bits must fit into an unsigned long" );
      errorCond( generatorPoly< (1 << (crcBits+1)),
		 "  degree of generator polygon is too high" );
    }
    
    /** encode a data value */
    unsigned long encode( unsigned long data ) const;
    
    /** decode a value
	\param code: the value to be decoded
	\param data: the decoded data (result)
	\returns whether or not the deconding was succesful */
    inline bool decode( unsigned long code, unsigned long &data ) const
    {
      data= code >> cBits;
      return encode( data ) == code;
    }
    
    /** validate a code */
    inline bool validate( unsigned long code ) const
    {
      unsigned long tmp;
      return decode( code, tmp );
    }
    
  protected:
      
    /** number of bits in a data word */
    unsigned dBits;
    
    /** number of bits in the code section */
    unsigned cBits;
    
    /** generator polygon (including leading 1) */
    unsigned long poly;
      
  };


} /* namespace */

#endif /* BASE_CRCCODE_H */

