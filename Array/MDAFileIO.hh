// ==========================================================================
// $Id: MDAFileIO.hh 309 2009-04-05 10:51:37Z heidrich $
// Input and output of MDA streams
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

#ifndef MDAFILEIO_H
#define MDAFILEIO_H

/*! \file  MDAFileIO.hh
    \brief Input and output of MDA streams
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <fstream>


#include "MDA/Base/BitsAndBytes.hh"
#include "MDA/Base/Types.hh"
#include "MDA/Base/CoordinateVector.hh"
#include "MDA/Base/ChannelList.hh"


namespace MDA {

  using namespace std;
  
  /** \class MDA_IO MDAFileIO.hh
      Abstract base class for reading/writing MDA streams */
  class MDA_IO {
    
  public:
    
    /** connect to provided file */
    virtual bool connect( const char *fileName )= 0;
    
    /** query dimensionality */
    inline const CoordinateVector &getDim() const
    {
      return dim;
    }
    
    /** query data type  */
    inline DataType getType() const
    {
      return dataType;
    }
    
    /** query number of channels */
    inline unsigned int getNumChannels() const
    {
      return numChannels;
    }
    
    /** return full list of all channels */
    inline ChannelList allChannels() const
    {
      ChannelList result;
      for( unsigned i= 0 ; i< numChannels ; i++ )
	result.vec.push_back( i );
      return result;
    }
    
    /** return the size, in bytes, of a scanline
	(only meaningful if connected) */
    inline unsigned long getScanlineSize() const
    {
      return scanlineSize;
    }
    
    /** return the number of scanlines that still need to be processed
	(only meaningful if connected) */
    inline unsigned long getNumScanlinesLeft() const
    {
      return numScanlinesLeft;
    }
    
    /** query data type */
    inline DataType getDataType() const
    {
      return dataType;
    }
    
    
  protected:
    
    /** calculate scanline size and number of scanlines */
    void setupScanlineInfo();
    
    /** size (in bytes) of a scanline */
    unsigned long scanlineSize;
    
    /** number of scanlines left */
    unsigned long numScanlinesLeft;
    
    /** dimensionality */
    CoordinateVector dim;
    
    /** number of channels */
    unsigned int numChannels;
    
    /** data type of the stream */
    DataType dataType;

    /** byte order of data stream */
    ByteOrder byteOrder;
  };



  /** \class MDAReader MDAFileIO.hh
      Reading MDA streams */
  
  class MDAReader: public MDA_IO {

  public:

    /** default constructor */
    MDAReader();

    /** destructor */
    ~MDAReader();
    
    /** connect to provided istream */
    bool connect( istream &is );
    
    /** constructor: open provided file */
    virtual bool connect( const char *fileName );
    
    /** is the reader connected to a vlaid stream? */
    inline bool isConnected() const
    {
      return (inStream!= NULL && inStream->good());
    }
    
    /** query if object is ready for reading scanlines */
    inline bool isReady() const
    {
      return (numScanlinesLeft!= 0) && isConnected();
    }
    
    /** return the number of scanlines that remain to be read */
    inline unsigned long getNumScanlinesLeft() const
    {
      return numScanlinesLeft;
    }
    
    /** disconnect from a stream */
    void disconnect();
    
    /** read MDA header from the istream we are connected to */
    bool readHeader();
    
    /** read one scanline (NULL: not ready) */
    void *readScanline();
    
    /** skip one scanline */
    bool skipScanline();

    /** true if the MDA file has opposite byte order to the CPU */
    bool endianSwap();

    /** allow external code to directly access the data pointer */
    istream* getStreamPointer();
    
  protected:
    
    /** istream to read the data from */
    istream *inStream;
    
  private:
    
    /** a buffer for one scanline */
     char *scanline;
    
    /** did we open the istream ourselves? */
    bool  ownIstream;
    
    /** private istream, created from file */
    ifstream myIstream;
  };



  /** \class MDAWriter MDAFileIO.hh
      Writing MDA streams */
  
  class MDAWriter: public MDA_IO {

  public:

    /** default constructor */
    MDAWriter();
    
    /** denstructor */
    ~MDAWriter();
    
    /** connect to provided ostream */
    bool connect( ostream &os );
    
    /** constructor: open provided file */
    virtual bool connect( const char *fileName );
    
    /** is the reader connected to a vlaid stream? */
    inline bool isConnected()
    {
      return (outStream!= NULL && outStream->good());
    }
    
    /** query if object is ready for reading scanlines */
    inline bool isReady()
    {
      return (numScanlinesLeft!= 0)  && isConnected();
    }
    
    /** return the number of scanlines that remain to be written */
    inline unsigned long getNumScanlinesLeft()
    {
      return numScanlinesLeft;
    }
    
    /** disconnect from a stream (returns false if write incomplete) */
    bool disconnect();
    
    /** read MDA header from the istream we are connected to */
    bool writeHeader( const CoordinateVector &dimension,
		      unsigned int channels, DataType dType,
		      ByteOrder endianness= UnknownEndian );
    
    /** read one scanline (false: not ready) */
    bool writeScanline( void *scanline );
    
    
  protected:

    /** ostream to write the data to */
    ostream *outStream;
    
  private:
    
    /** did we open the istream ourselves? */
    bool ownOstream;
    
    /** private istream, created from file */
    ofstream myOstream;
  };
} /* namespace */



#endif /* MDAFILEIO_H */

