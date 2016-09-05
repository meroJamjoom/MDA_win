// ==========================================================================
// $Id: Array.hh 653 2010-03-20 23:43:18Z heidrich $
// An n-dimensional array, separated into channels
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

#ifndef ARRAY_ARRAY_H
#define ARRAY_ARRAY_H

/*! \file  Array.hh
    \brief An n-dimensional array, separated into channels
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>

#include <MDA/Base/Types.hh>
#include <MDA/Base/ChannelList.hh>
#include "MemoryObject.hh"
#include "CoordinateArithmetic.hh"

namespace MDA {

  using namespace std;
  
  /** \class Array Array.hh
      An n-dimensional array, separated into channels */
  template<class T>
  class Array {

  public:
    
    
    /** \class Channel Array.hh
	One channel of an n-dimensional array */
    class Channel: public MemoryObject {
      
    public:
      
      /** constructor */
      Channel( const Array<T> *array, T background );
      
      /** read-only access to channel data */
      inline T operator[]( unsigned long index ) const
      {
	return data[index];
      }
      
      /** read/write access to channel data */
      inline T& operator[]( unsigned long index )
      {
	return data[index];
      }
      
      /** set background value */
      inline void setBackground( T value )
      {
	// data[-1] is still inside the allocated raw memory; it is
	// used for storing the background color
	data[-1]= value;
      }
      
      /** return background value */
      inline T getBackground()
      {
	// data[-1] is still inside the allocated raw memory; it is
	// used for storing the background color
	return data[-1];
      }
      
    protected:
      
      /** pointer to the array that this channel belongs to */
      const Array* array;
      
      /** the channel data - same as the double * of MemoryObject, but
	  of type T */
      T *data;
      
    private:
      
      /** assignment operator: doesn't exist */
      Channel &operator=( const Channel &ch );
      
      /** no copy constructor, either */
      Channel( const Channel &ch );
      
      /** destructor is private - clients use deref() */
      virtual ~Channel()
      {}
      
    };
    
    
    /** default constructor */
    inline Array()
      : nativeType( type )
    {}
    
    /** constructor for specific dimensions */
    inline Array( const CoordinateVector &dimensions )
      : dim( dimensions ), nativeType( type )
    {}
    
    /** copy constructor */
    inline Array( const Array &orig )
      : dim( orig.dim ), nativeType( orig.nativeType )
    {
      // create references to each of the channels
      for( unsigned i= 0 ; i< orig.channels.size() ; i++ )
	channels.push_back( (Channel *)(orig.channels[i]->ref()) );
    }
    
    /** destructor */
    ~Array();
    
    /** return "native" data type of the array */
    inline DataType getNativeType() const
    {
      return nativeType;
    }
    
    /** set "native" data type */
    inline void setNativeType( DataType newType )
    {
      nativeType= newType;
    }
    
    /** return the dimansion of the array */
    inline const CoordinateVector &getDimension() const
    {
      return dim;
    }
    
    /** return number of channels in the array */
    inline int getNumChannels() const
    {
      return channels.size();
    }
    
    /** get a scanline from the Array data structure (all channels) */
    void getScanline( unsigned long index, T* scanline ) const;
    
    /** update a scanline in the Array data structure (all channels) */
    void updateScanline( unsigned long index, const T* scanline );
    
    /** read array from MDA stream */
    bool read( istream &is= cin );
    
    /** read array from MDA file */
    bool read( const char *fileName );
    
    /** write array to MDA stream */
    bool write( ostream &os= cout, DataType outType= type );
    
    /** write array to MDA file */
    bool write( char *fileName, DataType outType= type );
    
    /** return a channel */
    inline Channel *operator[]( unsigned long i )
    {
      return i< channels.size() ? channels[i] : NULL;
    }
    
    /** add a new channel */
    inline unsigned addChannel( T background= 0 )
    {
      Channel *h= new Channel( this, background );
      channels.push_back( h );
      return channels.size()-1; // return index of new channel
    }
    
    /** produce a ChannelList with the IDs of all channels */
    inline ChannelList allChannels()
    {
      unsigned numChannels= channels.size();
      ChannelList ch( numChannels ); // empty, but with memory allocated
      for( unsigned i= 0 ; i< numChannels ; i++ )
	ch.vec.push_back( i );
      return ch;
    }
    
    /** produce a AxisList with the IDs of all axes */
    inline AxisList allAxes()
    {
      unsigned numAxes= dim.vec.size();
      AxisList axes( numAxes ); // empty, but with memory allocated
      for( unsigned i= 0 ; i< numAxes ; i++ )
	axes.vec.push_back( i );
      return axes;
    }
    
    /** remove a channel */
    inline void deleteChannel( unsigned index )
    {
      if( index< channels.size() )
      {
	channels[index]->deref();
	channels.erase( channels.begin()+index );
      }
    }
    
    /** swap two channels */
    inline void swapChannels( unsigned index1, unsigned index2 )
    {
      if( index1< channels.size() && index2< channels.size() )
      {
	Channel *tmp= channels[index1];
	channels[index1]= channels[index2];
	channels[index2]= tmp;
      }
    }
    
  protected:
    
    /** a list of all channels associated with this array */
    vector<Channel *> channels;
    
    /** dimensions of the array */
    CoordinateVector dim;
    
    /** "native" data type; e.g. the type read from a file */
    DataType nativeType;
    
    /** the DataType corresponding to the template argument T */
    static const DataType type;
  };
  
  typedef Array<float> FloatArray;
  typedef Array<double> DoubleArray;
  typedef Array<float>::Channel FloatChannel;
  typedef Array<double>::Channel DoubleChannel;

  
  
} /* namespace */


#endif /* ARRAY_ARRAY_H */

