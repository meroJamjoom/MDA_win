// ==========================================================================
// $Id:$
// Support OpenCV library with MDA files
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: atcheson (Bradley Atcheson)
// Email:   atcheson@cs.ubc.ca
// ==========================================================================

#ifndef BRIDGE_OPENCVBRIDGE_H
#define BRIDGE_OPENCVBRIDGE_H

/*! \file  OpenCVBridge.hh
    \brief Support OpenCV library with MDA files
 */

#if defined ( _WIN32) || defined (_WIN64)

// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>
#include <cv.h>
#include <highgui.h>
#include "MDA/Array/MDAFileIO.hh"


// Note that I wrote this bridge first, and was just experimenting with
// interfaces. The CImg has a nicer structure I think.


namespace MDA {

  /** \class OpenCVReadBridge OpenCVBridge.hh
   *  Connect OpenCV library to MDA files. OpenCV deals primarily with 2D
   *  images, so we assume that the first two dimensions of an MDA correspond
   *  to a 2D "frame" and these are extracted, in sequence, using this bridge.
   *  For example, a 640x480x10x99 MDA would have 10x99 = 990 frames (this
   *  could be a sequence of 99 volumetric data arrays). Note that OpenCV
   *  is limited to between 1 and 4 channels, and many functions require
   *  uint8 or float32 data formats, so you may need to convert the MDA first
   *  or else make sure to convert the OpenCV arrays.
   */
  class OpenCVReadBridge {

  public:

    /** Default constructor */
    OpenCVReadBridge();

    /** Connect to an MDA stream */
    bool connect( istream& is=cin );
    /** Connect to an MDA file */
    bool connect( const char* filename );
    /** Disconnect from current MDA object */
    void disconnect();

    /** Allocate memory for a new OpenCV buffer to hold one frame */
    IplImage* createEmptyFrame();
    /** Extract next frame from the MDA object and copy into OpenCV buffer
     *  \param "img" destination OpenCV buffer. An new empty buffer will be
     *               allocated if this is NULL. 
     *  \return pointer to OpenCV array, or NULL if no more frames in MDA */
    IplImage* extractNextFrame( IplImage* img=NULL );
    /** Skip past next frame in MDA object */
    void skipNextFrame();

    /** Total number of frames in MDA */
    unsigned getNumFrames();

    /** True if any more frames exist in sequence */
    bool moreFrames();

    /** How many channels in MDA
     *  Note that all of these member access functions will give useless
     *  results if you try to call them before a successful connection
     *  to an MDA is mda */
    unsigned getNumChannels() {
      return channels;
    }

    /** Dimensions of MDA */
    CoordinateVector getDim() {
      return reader.getDim();
    }

    /** Data type of the connected MDA stream */
    DataType getDataType() {
      return reader.getDataType();
    }

  protected:

    /** Check that MDA object is valid */
    bool validateHeader();

    /** Read frame dimensions and data format from MDA */
    void configure();

  protected:

    MDAReader reader;

    /** Width of frame */
    unsigned width;
    /** Height of frame */
    unsigned height;
    /** Number of channels per frame */
    unsigned channels;
    /** Data format of OpenCV array */
    unsigned depth;
    /** MDA data format to convert to so as to be compatible with 'depth' */
    DataType targetType;

  };


  /** \class OpenCVWriteBridge OpenCVBridge.hh
   *  Convert OpenCV arrays to MDA streams. Since OpenCV happens to store
   *  multichannel data in an interleaved format, and MDA files can store
   *  all of the native OpenCV types, this is quite a simple and efficient
   *  operation - we can just copy the scanlines across directly to the
   *  output stream.
   */
  class OpenCVWriteBridge
  {

  public:

    /** Constructor. This ensures that the output type matches that of the
     *  OpenCV buffers we will be converting from.
     *  \param "img" One of the OpenCV frames that we will be writing
     *  \param "nFrames" The number of 2D frames that will be output
     */
    OpenCVWriteBridge( IplImage* img, unsigned nFrames );

    /** Connect to the output MDA stream */
    bool connect( ostream& os=cout );
    /** Connect to an output MDA file */
    bool connect( const char* filename );
    /** Disconnect from current MDA object */
    void disconnect();

    /** Append the OpenCV image to the output MDA stack. This automatically
     *  ignores any padding that OpenCV may have added to the scanlines.
     *  \param "img" The 2D OpenCV image to write.
     *  \return False on error.
     */
    bool writeFrame( const IplImage* img );

  protected:

    MDAWriter writer;

    /** Dimensions of output MDA stack (always 3D) */
    CoordinateVector dim;

    /** Number of channels in output MDA */
    unsigned channels;

    /** Data type of the input OpenCV frames */
    unsigned depth;

    /** Data type of output MDA */
    DataType type;

  };

} /* namespace */



#endif /* BRIDGE_OPENCVBRIDGE_H */

