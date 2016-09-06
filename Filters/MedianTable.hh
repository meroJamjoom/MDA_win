// ==========================================================================
// $Id: MedianTable.hh 999 2014-05-28 15:07:31Z heidrich $
// a table useful for the incremental computation of medians
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

#ifndef FILTERS_MEDIANTABLE_H
#define FILTERS_MEDIANTABLE_H

/*! \file  MedianTable.hh
    \brief a table useful for the incremental computation of medians
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>

namespace MDA {
  
  using namespace std;
  
  /** \class MedianTable MedianTable.hh
      a table useful for the incremental computation of medians */
  template<class T>
  class MedianTable {

  public:
    
    /** default constructor */
    MedianTable( T min= 0.0, T max= 1.0, int size= 256 )
    {
      setRange( min, max, size );
    }
    
    /** destructor */
    virtual ~MedianTable()
    {}
    
    /** set table range (also clears the data structure ) */
    void setRange( T min, T max, int size= 256 );
    
    /** clear the data structure */
    virtual void clear()
    {}
    
    /** add a value to the data structure */
    inline void add( T val )
    {
      add( val, getHash( val ) );
    }
    
    /** add a value to the data structure (using precomputed hash) */
    virtual void add( T val, short hashIndex )= 0;
    
    /** remove a value from the data structure */
    inline bool remove( T val )
    {
      return remove( val, getHash( val ) );
    }
    
    /** remove a value from the data structure (using precomputed hash) */
    virtual bool remove( T val, short hashIndex )= 0;
    
    /** replace a value in the data structure */
    inline bool replace( T oldVal, T newVal )
    {
      return replace( oldVal, getHash( oldVal ), newVal, getHash( newVal ) );
    }
    
    /** replace a value in the data structure (using precomputed hashs) */
    virtual bool replace( T oldVal, short oldHash, T newVal, short newHash )=0;
    
    /** find the median of all values currently in the data structure */
    virtual T getMedian()= 0;
    
    /** compute table index for a given value */
    inline short getHash( T value )
    {
      short index= (int)((value-minVal) * mult);
      return index>= tableSize-1 ? tableSize-1 : (index< 0 ? 0 : index);
    }
    
  protected:

    /** the smallest value in the range */
    T minVal;
      
    /** the largest value in the range */
    T maxVal;
    
    /** entry-to-entry increment in table */
    T mult;
    
    /** the number of entries in the table */
    int tableSize;
    
    /** whether the memory needs to be allocated before use, or
	whether it exists already */
    bool needsAllocation;
  };
  
  

  /** \class ContinuousMedianTable MedianTable.hh
      a hash datastructure for fast median filtering (contiunuous data) */
  template<class T>
  class ContinuousMedianTable: public MedianTable<T> {
    
  public:
    
    /** default constructor */
    ContinuousMedianTable( T min= 0.0, T max= 1.0, int size= 256 )
      : table( NULL ), counters( NULL ), MedianTable<T>( min, max, size )
    {}
    
    /** destructor */
    virtual ~ContinuousMedianTable()
    {
      delete [] counters;
      delete [] table;
    }
    
    /** add a value to the data structure */
    virtual void add( T val, short hashValue );
    
    /** remove a value from the data structure */
    virtual bool remove( T val, short hashValue );
    
    /** remove a value from the data structure */
    virtual bool replace( T oldVal, short oldHash, T newVal, short newHash )
    {
      bool status= remove( oldVal, oldHash );
      add( newVal, newHash );
      return status;
    }
    
    /** clear the data structure */
    virtual void clear();
    
    /** find the median of all values currently in the data structure */
    virtual T getMedian();
    
    
  protected:
    
    /** an entry in the data structure */
    struct HashEntry {
    public:
      HashEntry( T val )
	: value(val), count(1)
      {}
      T value;
      int count;
    };
    
    /** the table */
    list<HashEntry> *table;
    
    /** counters for the number of elements per table entry */
    int *counters;
    
    /** the last reported median value */
    T lastMedian;
    
    /** number of elements in the table */
    int numElements;
    
    /** number of elements to the left of the last reported median */
    int numSmaller;
    
    /** number of elements to the right of the last reported median */
    int numLarger;
    
  };
  
  
  /** \class QuantizedMedianTable MedianTable.hh
      a hash datastructure for fast median filtering (quantized data) */
  template<class T>
  class QuantizedMedianTable: public MedianTable<T> {
    
  public:
    
    /** default constructor */
    QuantizedMedianTable( T min= 0.0, T max= 1.0, int size= 256 )
      : table( NULL ), values( NULL ), MedianTable<T>( min, max, size )
    {}
    
    /** destructor */
    virtual ~QuantizedMedianTable()
    {
      delete [] values;
      delete [] table;
    }
    
    /** add a value to the data structure */
    virtual void add( T val, short hashValue )
    {
      table[hashValue]++;
      numElements++;
      if( hashValue< lastMedian )
	numSmaller++;
    }
    
    /** remove a value from the data structure */
    virtual bool remove( T val, short hashValue )
    {
      table[hashValue]--;
      numElements--;
      if( hashValue< lastMedian )
	numSmaller--;
      return true;
    }
    
    /** replace a value in the data structure */
    virtual bool replace( T oldVal, short oldHash, T newVal, short newHash )
    {
      table[oldHash]--;
      table[newHash]++;
      if( oldHash< lastMedian )
	numSmaller--;
      if( newHash< lastMedian )
	numSmaller++;
      return true;
    }
    
    /** clear the data structure */
    virtual void clear()
    {
      // (re)allocate memory if necessary
      if( MedianTable<T>::needsAllocation )
      {
	if( values!= NULL )
	  delete [] values;
	if( table!= NULL )
	  delete [] table;
	table= new int[MedianTable<T>::tableSize];
	values= new T[MedianTable<T>::tableSize];
	MedianTable<T>::needsAllocation= false;
      }
      
      // then initialize memory
      for( int i= 0 ; i< MedianTable<T>::tableSize ; i++ )
      {
	table[i]= 0;
	values[i]=  MedianTable<T>::minVal + i/MedianTable<T>::mult;
      }
      numElements= numSmaller= 0;
      lastMedian= MedianTable<T>::tableSize/2;
    }
    
    /** find the median of all values currently in the data structure */
    virtual T getMedian();
    
    
  protected:
    
    /** the table */
    int *table;
    
    /** counters for the number of elements per table entry */
    T *values;
    
    /** the last reported median value */
    int lastMedian;
    
    /** number of elements in the table */
    int numElements;
    
    /** number of elements to the left of the last reported median */
    int numSmaller;
    
  };

} /* namespace */

#endif /* FILTERS_MEDIANTABLE_H */

