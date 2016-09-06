// ==========================================================================
// $Id: MedianTable.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef FILTERS_MEDIANTABLE_C
#define FILTERS_MEDIANTABLE_C

#include <iostream>

#include "MDA/Base/Errors.hh"
#include "MedianTable.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions

/** set table range (also clears the data structure ) */
template<class T>
void
MedianTable<T>::setRange( T min, T max, int size )
{
  tableSize= size;
  minVal= min;
  maxVal= max;
  mult= (tableSize-1) / (double)(maxVal-minVal);
  needsAllocation= true;
}


//
// ContinuousMedianTable members
//

/** add a value to the data structure */
template<class T>
void
ContinuousMedianTable<T>::add( T val, short index )
{
  typename list<HashEntry>::iterator iter;
  
  // skip all list entries less than the specified value
  for( iter= table[index].begin() ;
       iter!= table[index].end() && val> iter->value ;
       iter++ )
    ;
  
  if( iter!= table[index].end() && val== iter->value )
    // entry with same value already exists
    iter->count++;
  else
    // create new entry
    table[index].insert( iter, HashEntry( val ) );

  // update counters
  numElements++;
  counters[index]++;
  short lastHash= this->getHash( lastMedian );
  if( index< lastHash )
    numSmaller++;
  if( index> lastHash )
    numLarger++;
}
    
/** remove a value from the data structure */
template<class T>
bool
ContinuousMedianTable<T>::remove( T val, short index )
{
  typename list<HashEntry>::iterator iter;
  
  // skip all list entries less than the specified value
  for( iter= table[index].begin() ;
       iter!= table[index].end() && val> iter->value ;
       iter++ )
    ;
  
  // the current entry should now be the one we seek. If not, we try
  // to remove an entry that has never been added
  if( iter== table[index].end() || iter->value!= val )
  {
    warning( "Not found!!\n" );
    return false;
  }
  // decrement counter, and delete entry if counter reaches 0
  if( --(iter->count)<= 0 )
    table[index].erase( iter );
  
  // update counters
  numElements--;
  counters[index]--;
  short lastHash= this->getHash( lastMedian );
  if( index< lastHash )
    numSmaller--;
  if( index> lastHash )
    numLarger--;
  
  return true;
}


/** clear the data structure */
template<class T>
void
ContinuousMedianTable<T>::clear()
{
  // (re)allocate memory if necessary
  if( MedianTable<T>::needsAllocation )
  {
    if( counters!= NULL )
      delete [] counters;
    if( table!= NULL )
      delete [] table;
    table= new list<HashEntry>[MedianTable<T>::tableSize];
    counters= new int[MedianTable<T>::tableSize];
    MedianTable<T>::needsAllocation= false;
  }
  
  // initialize the table
  for( int i= 0 ; i< MedianTable<T>::tableSize ; i++ )
  {
    table[i].clear();
    counters[i]= 0;
  }
  numElements= numSmaller= numLarger= 0;
  lastMedian= (MedianTable<T>::minVal+MedianTable<T>::maxVal)/2.0;
}

/** find the median of all values currently in the data structure */
template<class T>
T
ContinuousMedianTable<T>::getMedian()
{
  // find correct table index
  short index= this->getHash( lastMedian );
  while( numSmaller> numElements/2 && index> 0 )
  {
    numLarger+= counters[index];
    numSmaller-= counters[--index];
  }
  while( (numLarger> numElements/2 || counters[index]== 0) &&
	 index< MedianTable<T>::tableSize-1 )
  {
    numSmaller+= counters[index];
    numLarger-= counters[++index];
  }

  // once we have the right table index we search the corresponding
  // list
  int skip= (numElements+1)/2 - numSmaller;
  typename list<HashEntry>::iterator iter;
  for( iter= table[index].begin() ; iter!= table[index].end() ; iter++ )
    if( iter->count>= skip )
      return lastMedian= iter->value;
    else
      skip-= iter->count;
  
  // should never be reached
  warning( "  median not found" );
  return lastMedian;  // keeps compiler happy
}








//
// QuantizedMedianTable members
//

/** find the median of all values currently in the data structure */
template<class T>
T
QuantizedMedianTable<T>::getMedian()
{
  // find correct table index
  while( numSmaller> numElements/2 && lastMedian> 0 )
    numSmaller-= table[--lastMedian];
  while( (numElements/2> numSmaller+table[lastMedian] ||
	  table[lastMedian]== 0 ) &&
	 lastMedian< MedianTable<T>::tableSize-1 )
    numSmaller+= table[lastMedian++];
  
  return values[lastMedian];
}


template class MedianTable<float>;
template class MedianTable<double>;
template class ContinuousMedianTable<float>;
template class ContinuousMedianTable<double>;
template class QuantizedMedianTable<float>;
template class QuantizedMedianTable<double>;


} /* namespace */

#endif /* FILTERS_MEDIANTABLE_C */

