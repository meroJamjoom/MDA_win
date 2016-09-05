// ==========================================================================
// $Id: Range.hh 409 2009-10-14 21:37:44Z heidrich $
// range of values and list of ranges
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

#ifndef BASE_RANGE_H
#define BASE_RANGE_H

/*! \file  Range.hh
    \brief range of values and list of ranges
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <utility>
#include <vector>

namespace MDA {

  using namespace std;
  
  /** \class Range Range.hh
      range of values */
  template <class T> class Range {
  public:
    /** constructor from pair of values */
    inline Range( T a= 0, T b= 0 )
    {
      val.first= a;
      val.second= b;
    }
    /** the actual values */
    pair<T,T> val;
    /** empty range in type T */
    static const Range Empty;
    /** full range in type T */
    static const Range Full;
  };
  
  /** range of integers */
  typedef Range<int> IntRange;
  
  /** range of non-negative integers */
  typedef Range<unsigned int> UIntRange;
  
  /** range of long integers */
  typedef Range<long> LongRange;
  
  /** range of non-negative long integers */
  typedef Range<unsigned long> ULongRange;
    
  /** range of floats */
  typedef Range<float> FloatRange;
  
  /** clamp a value to a range */
  template <class T> T clampToRange( const T &value, const Range<T> &range );
  
  /** write range to ostream */
  template <class T> ostream &operator<<( ostream &os, const Range<T> &range );
  
  /** read range from istream */
  template <class T> istream &operator>>( istream &is, Range<T> &range );
  
  
  /** \class RangeList Range.hh
      list of ranges */
  template <class T> class RangeList {
  public:
    /** constructor for empty range vector */
    inline RangeList() {}
    /** constructor for range vector of a certain dimension
	(initialized to Empty) */
    inline RangeList( unsigned dimension )
      : vec( dimension )
    {
      for( unsigned i= 0 ; i< dimension ; i++ )
	vec[i]= Range<T>::Empty;
    }
    /** the actual list */
    vector< Range<T> > vec;
  };
  
  /** list of IntRange */
  typedef RangeList<int> IntRangeList;
  
  /** list of UIntRange */
  typedef RangeList<unsigned int> UIntRangeList;
  
  /** list of LongRange */
  typedef RangeList<long> LongRangeList;
  
  /** list of ULongRange */
  typedef RangeList<unsigned long> ULongRangeList;
  
  /** List of FloatRange */
  typedef RangeList<float> FloatRangeList;
  
  /** write range list to ostream */
  template <class T> ostream &operator<<( ostream &os,
					  const RangeList<T> &rangeL );
  
  /** read range list from istream */
  template <class T> istream &operator>>( istream &is, RangeList<T> &rangeL );
  

} /* namespace */


#ifndef BASE_RANGE_C
#ifndef BASE_RANGE_TEMPLATE
#define BASE_RANGE_TEMPLATE
#include "Range.C"
#endif
#endif


#endif /* BASE_RANGE_H */

