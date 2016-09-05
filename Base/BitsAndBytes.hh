// ==========================================================================
// $Id: BitsAndBytes.hh 661 2010-03-22 01:19:14Z heidrich $
// Functions dealing with bit representations and endianness
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2006-, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef BASE_BITSANDBYTES_H
#define BASE_BITSANDBYTES_H

/*! \file  BitsAndBytes.hh
    \brief Functions dealing with bit representations and endianness
 */

#if defined(_WIN32) || defined(_WIN64)
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#define strncasecmp _strnicmp
#define strcasecmp _stricmp

#endif

#include <iostream>


#include "Errors.hh"
#include "CommandlineParser.hh"

namespace MDA {

  using namespace std;
  
  /** the endianness of a CPU or data field */
  enum ByteOrder
  {
    UnknownEndian,
    LittleEndian,
    BigEndian
  };
  
  /** endianness of the CPU */
  inline ByteOrder cpuByteOrder()
  {
    int test= 0xFFAA8844;
    unsigned char *bytes= (unsigned char *)&test;
    
    if( bytes[0]== 0xFF )
    {
      errorCond( bytes[1]== 0xAA && bytes[2]== 0x88 && bytes[3]== 0x44,
		 "  Byte order seems to be neither big nor little endian" );
      return BigEndian;
    }
    else
    {
      errorCond( bytes[3]== 0xFF && bytes[2]== 0xAA &&
		 bytes[1]== 0x88 && bytes[0]== 0x44,
		 "  Byte order seems to be neither big nor little endian" );
      return LittleEndian;
    }
  }
  
  /** write endianness to ostream */
  ostream &operator<<( ostream &os, ByteOrder byteOrder );
  
  /** read endianness from istream */
  istream &operator>>( istream &is, ByteOrder &byteOrder );
  
  
  /** swap endianness (in place) */
  void *swapByteOrder( void *bytes, unsigned numBytes );

  /** swap the endianness of a whole array (in place) */
  void *swapArrayByteOrder( void *bytes, unsigned long entries,
			    unsigned bytesPerEntry );

  /** is a given number a power of two? */
  bool isPowerOf2( unsigned long val );

  /** the smallest power of two greater or equal to the given value */
  unsigned long nextPowerOf2( unsigned long val );

  /** Hamming distance (#bits that differ) between 2 unsigned longs
      (with one argument: count number of 1-bits)
   */
  unsigned hammingDist( unsigned long v1, unsigned long v2= 0ul );
  

  /** \class ByteOrderOption BitsAndBytes.hh
      parser for byte order id */
  class ByteOrderOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to DataType object */
    ByteOrderOption( ByteOrder &o,
		     const char *msg="\tbyte order (little, big, or unknown)\n",
		     const char *longOpt= "--byte-order",
		     const char *shortOpt= "-b" )
      : CommandlineOption( msg, longOpt, shortOpt ), order( o )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
  
    /** reference to the DataType object */
    ByteOrder &order;
  };
  
  
  
} /* namespace */



#endif /* BASE_BITSANDBYTES_H */

