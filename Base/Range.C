// ==========================================================================
// $Id: Range.C 411 2009-10-14 23:27:26Z heidrich $
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

#ifndef BASE_RANGE_C
#define BASE_RANGE_C

#include <iostream>
#include <limits.h>
#include <float.h>

#include "Errors.hh"
#include "Range.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

// full and empty ranges for all types
// (only if not included just for template code)
#ifndef BASE_RANGE_TEMPLATE
template<> const Range<int>
Range<int>::Empty( INT_MAX, INT_MIN );
template<> const Range<int>
Range<int>::Full( INT_MIN, INT_MAX );
template<> const Range<unsigned>
Range<unsigned>::Empty( INT_MAX, 0 );
template<> const Range<unsigned>
Range<unsigned>::Full( 0, INT_MAX );
template<> const Range<long>
Range<long>::Empty( LONG_MAX, LONG_MIN );
template<> const Range<long>
Range<long>::Full( LONG_MIN, LONG_MAX );
template<> const Range<unsigned long>
Range<unsigned long>::Empty( LONG_MAX, 0ul );
template<> const Range<unsigned long>
Range<unsigned long>::Full( 0ul, LONG_MAX );
template<> const Range<float>
Range<float>::Empty( FLT_MAX, FLT_MIN );
template<> const Range<float>
Range<float>::Full( FLT_MIN, FLT_MAX );
#endif


/** write range to ostream */
template <class T> ostream &
operator<<( ostream &os, const Range<T> &range )
{
  return os << range.val.first << ':' << range.val.second;
}
  
/** read range from istream */
template <class T> istream &
operator>>( istream &is, Range<T> &range )
{
  is >> range.val.first;
  if( is.peek()== ':' )
  {
    is.get();
    is >> range.val.second;
  }
  else
    is.setstate( ios::failbit );

  return is;
}

/** clamp a value to a range */
template <class T> T
clampToRange( const T &value, const Range<T> &range )
{
  return value< range.val.first ? range.val.first : 
    value > range.val.second ? range.val.second : value;
}



/** write range list to ostream */
template <class T> ostream &
operator<<( ostream &os, const RangeList<T> &rangeL )
{
  int numEntries= rangeL.vec.size();
  if( numEntries> 0 )
  {
    os << rangeL.vec[0];
    for( int i = 1 ; i< numEntries ; i++ )
      os << ',' << rangeL.vec[i];
  }
  
  return os;
}
  
/** read range list from istream */
template <class T> istream &
operator>>( istream &is, RangeList<T> &rangeL )
{
  Range<T> range;
  
  rangeL.vec.clear();
  do
  {
    // read one raneg, and append it to the vector
    is >> range;
    if( is.fail() )
      break;
    rangeL.vec.push_back( range );
    
    // stop parsing if we no longer get commas
    if( is.eof() )
      break;
    if( is.peek()== ',' )
      is.get();
    else
      break;
  }
  while( !is.fail() );

  warnCond( !is.fail(), "  could not read range!" );
  return is;
}



} /* namespace */

#endif /* BASE_RANGE_C */

