// ==========================================================================
// $Id: Array.C 665 2010-03-23 05:20:27Z martinle $
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

#ifndef ARRAY_ARRAY_C
#define ARRAY_ARRAY_C

#include <fstream>

#include "MDAFileIO.hh"
#include "Array.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** mapping of template argument T to DataType
    (this list needs to be augmented for other possible template types)
*/
template<> const DataType Array<unsigned char>::type=	UByte;
template<> const DataType Array<char>::type=		Byte;
template<> const DataType Array<unsigned short>::type=	UShort;
template<> const DataType Array<short>::type=		Short;
template<> const DataType Array<unsigned int>::type=	UInt;
template<> const DataType Array<int>::type=		Int;
template<> const DataType Array<float>::type=		Float;
template<> const DataType Array<double>::type=		Double;

//
// Channel members
//

/** Channel constructor */
template<class T>
Array<T>::Channel::Channel( const Array<T> *a, T b )
  : MemoryObject(), array( a )
{
  unsigned i;
  unsigned long size= 1;
  CoordinateVector dim= array->getDimension();
  
  for( i= 0 ; i< dim.vec.size() ; i++ )
    size*= dim.vec[i];
  
  // convert size from units of sizeof( T ) to units of sizeof( double )
  size= (sizeof( T ) * size + sizeof( double ) - 1) / sizeof( double );
  allocate( size );
  data= (T *)MemoryObject::data;
  
  // data[-1] is used to store the background value
  // this is always inside the memory region allocated by MemoryObject
  data[-1]= b;
}


//
// Array members
//


/** Array destructor */
template<class T>
Array<T>::~Array()
{
  for( int i= channels.size()-1 ; i>= 0 ; i-- )
    channels[i]->deref();
}

/** get a scanline from the Array data structure (all channels) */
template<class T>
void
Array<T>::getScanline( unsigned long index, T* scanline ) const
{
  int numChannels= channels.size();
  unsigned long numEntries= dim.vec[0];
  unsigned long channelIndex= index*numEntries;
  
  for( unsigned long i= 0 ; i< numEntries ; i++, channelIndex++ )
    for( unsigned j= 0 ; j< numChannels ; j++ )
      *scanline++= (*channels[j])[channelIndex];
}

/** update a scanline in the Array data structure (all channels) */
template<class T>
void
Array<T>::updateScanline( unsigned long index, const T* scanline )
{
  int numChannels= channels.size();
  unsigned long numEntries= dim.vec[0];
  unsigned long channelIndex= index*numEntries;
  
  for( unsigned long i= 0 ; i< numEntries ; i++, channelIndex++ )
    for( unsigned j= 0 ; j< numChannels ; j++ )
      (*channels[j])[channelIndex]= *scanline++;
}

/** read array from MDA stream */
template<class T>
bool
Array<T>::read( istream &is )
{
  unsigned long i;
  
  // delete any channels that exist right now
  for( i= channels.size() ; i> 0 ; )
    channels[--i]->deref();
  channels.clear();
  
  // read input header
  MDAReader reader;
  reader.connect( is );
  bool result =  reader.readHeader();
  if( !warnCond( result, "  cannot read file header\n" ) )
    return false;
  dim=		reader.getDim();
  nativeType=	reader.getType();
  unsigned int  numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  
  // create as many channels as needed
  for( i= 0 ; i< numChannels ; i++ )
    addChannel();
  
  // temporary memory holding one scanline already copnverted to type T
  unsigned long numEntries= dim.vec[0]*numChannels;
  T *scanline= new T[numEntries];
  
  // read individual scanlines and add them to the channel data
  for( i= 0 ; i< numScanlines ; i++ )
  {
    // read next scanline and convert to T
    typeConvert( reader.readScanline(), nativeType, scanline, type, numEntries);
    updateScanline( i, scanline );
  }
  
  // clean up
  delete [] scanline;
  reader.disconnect();
  return true;
}
    
/** read array from MDA file */
template<class T>
bool
Array<T>::read( const char *fileName )
{
  ifstream is;
  is.open( fileName, ifstream::in & ifstream::binary );
  bool result= is.good() && read( is );
  is.close();
  
  return result;
}
    
/** write array to MDA stream */
template<class T>
bool
Array<T>::write( ostream &os, DataType outType )
{
  // write in native type if the user hasn't specified a data type
  if( outType== UndefinedType )
    outType= nativeType;
  
  MDAWriter writer;
  writer.connect( os );
  bool result = writer.writeHeader( dim, channels.size(), outType );
  if( !warnCond( result,
		 "  cannot write file header\n" ) )
    return false;
  unsigned long numScanlines= writer.getNumScanlinesLeft();
  
  // temporary memory holding one scanline of type T
  unsigned long numEntries= dim.vec[0]*channels.size();
  T *scanline= new T[numEntries];
  // and another one to hold the scanline in the output format
  char *scanlineOut= new char[writer.getScanlineSize()];
  
  // read individual scanlines and add them to the channel data
  for( unsigned long i= 0 ; i< numScanlines ; i++ )
  {
    // get next scanline from array
    getScanline( i, scanline );
    // convert it to outType
    typeConvert( scanline, type, scanlineOut, outType, numEntries );
    // and write to stream
    writer.writeScanline( scanlineOut );
  }
  
  delete [] scanline;
  delete [] scanlineOut;
  return writer.disconnect();
}

/** write array to MDA file */
template<class T>
bool
Array<T>::write( char *fileName, DataType outType )
{
  ofstream os;
  os.open( fileName, ofstream::out & ifstream::binary );
  bool result= os.good() && write( os );
  os.close();
  
  return result;
}
    

/* we use explicit instantiation for this class */
template class Array<float>;
template class Array<double>;
template class Array<unsigned char>;


} /* namespace */

#endif /* ARRAY_ARRAY_C */

