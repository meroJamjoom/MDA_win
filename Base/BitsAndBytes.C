// ==========================================================================
// $Id: BitsAndBytes.C 648 2010-03-18 03:55:39Z heidrich $
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

#ifndef BASE_BITSANDBYTES_C
#define BASE_BITSANDBYTES_C

#include <iostream>
#include <string.h>

#include "BitsAndBytes.hh"

namespace MDA {

using namespace std;


static const int maxByteOrderStringLength= 8;
static const char byteOrderStrings[3][maxByteOrderStringLength]=
{
  "unknown", "little", "big"
};


/** write byte order to ostream */
ostream
&operator<<( ostream &os, ByteOrder byteOrder )
{
  return os << byteOrderStrings[byteOrder];
}
  
/** read byte order from istream */
istream
&operator>>( istream &is, ByteOrder &byteOrder )
{
  char buffer[maxByteOrderStringLength];

  is >> ws;
  is.width( maxByteOrderStringLength );
  is >> buffer;
  byteOrder= UnknownEndian;
  
  for( int i= 0 ; i< 3 ; i++ )
    if( !strncasecmp( buffer, byteOrderStrings[i], maxByteOrderStringLength ) )
    {
      byteOrder= (ByteOrder)i;
      return is;
    }

  // if we get here, things have gone wrong...
  warning( "  could not read byte order!" );
  return is;
}


/** swap byte order (in place) */
void *
swapByteOrder( void *bytes, unsigned numBytes )
{
  unsigned char *mem= (unsigned char *)bytes;
  unsigned char h;
  
  for( long i= 0 ; i< numBytes/2 ; i++ )
  {
    h= mem[i]; mem[i]= mem[numBytes-1-i]; mem[numBytes-1-i]= h;
  }

  return bytes;
}


/** swap the byte order of a whole array (in place) */
void *
swapArrayByteOrder( void *bytes, unsigned long entries, unsigned bytesPerEntry)
{
  unsigned char *mem= (unsigned char *)bytes;
  unsigned char h;
  long i, j, k;
  
  for( j= k= 0 ; j< entries ; j++, k+= bytesPerEntry )
    for( i= 0 ; i< bytesPerEntry/2 ; i++ )
    {
      h= mem[k+i];
      mem[k+i]= mem[k+bytesPerEntry-1-i];
      mem[k+bytesPerEntry-1-i]= h;
    }

  return bytes;
}


/** is a given number a power of two? */
bool
isPowerOf2( unsigned long val )
{
  for( ; val> 1 ; val>>= 1 )
    if( val&1 )
      return false;
  return val==1;
}

/** the smallest power of two greater or equal to the given value */
unsigned long
nextPowerOf2( unsigned long val )
{
  unsigned long ret;
  for( ret= 1 ; ret< val ; ret<<= 1 )
    ;
  return ret;
}

/** Hamming distance (#bits that differ) between 2 unsigned longs */
unsigned
hammingDist( unsigned long v1, unsigned long v2 )
{
  unsigned long val= v1 ^ v2;
  unsigned hamming= 0u;
  do
    if( val & 1 )
      hamming++;
  while( (val>>= 1)> 0ul );
  
  return hamming;
}



  //
  // ByteOrderOption functions
  // 

/** the actual parsing function */
bool
ByteOrderOption::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> order;
  return !optStr.fail();
}

/** output usage string */
void
ByteOrderOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <order>\n"
     << helpTxt
     << "\tCurrently: " << order << "\n\n";
}


} /* namespace */

#endif /* BASE_BITSANDBYTES_C */

