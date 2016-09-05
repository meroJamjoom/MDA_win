// ==========================================================================
// $Id: PaddedChannelTraverser.hh 384 2009-09-25 19:48:11Z heidrich $
// traversal class for pixels in padded arrays/channels
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef ARRAY_PADDEDCHANNELTRAVERSER_H
#define ARRAY_PADDEDCHANNELTRAVERSER_H

/*! \file  PaddedChannelTraverser.hh
    \brief traversal class for pixels in padded arrays/channels
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "CoordinateArithmetic.hh"

namespace MDA {

  /** \class PaddedChannelTraverser PaddedChannelTraverser.hh
      traversal class for pixels in padded arrays/channels */
  
  class PaddedChannelTraverser {

  public:

    /** default constructor */
    PaddedChannelTraverser( CoordinateVector &_dim, unsigned long padding );
    
    /** start traversal */
    inline void begin()
    {
      iter.begin();
      uOff= scanlinePos= 0;
      pOff= pixelOffset( iter.getPos(), paddedDim, padding );
    }
    
    /** check if we are done */
    inline bool isAtEnd() const
    {
      return uOff>= numPixels;
    }
    
    /** increment iterator */
    inline void operator++()
    {
      uOff++;
      if( ++scanlinePos>= scanlinePixels )
      {
	// increment position by one scanline, recompute padded offset
	iter.incrComp( 1 );
	pOff= pixelOffset( iter.getPos(), paddedDim, padding );
	scanlinePos= 0;
      }
      else
	pOff++;
    }
    
    /** go to start of next scanline */
    inline void nextScanline()
    {
      iter.incrComp( 1 );
      uOff= pixelOffset( iter.getPos(), dim );
      pOff= pixelOffset( iter.getPos(), paddedDim, padding );
    }
      
    /** return unpadded array offset */
    inline unsigned long uOffset() const
    {
      return uOff;
    }
    
    /** return padded array offset */
    inline unsigned long pOffset() const
    {
      return pOff;
    }
    
  protected:
    
    //
    // iterator state
    //
    
    /** unpadded array offset */
    unsigned long uOff;
    
    /** padded offset */
    unsigned long pOff;
    
    /** current position within scanline */
    unsigned long scanlinePos;
    
    //
    // basic dimensions
    //
        
    /** unpadded dimension */
    CoordinateVector dim;
    
    /** padding radius */
    unsigned long padding;
    
    /** padded dimension */
    CoordinateVector paddedDim;
    
    /** number of scanlines */
    unsigned long numScanlines;
    
    /** number of pixels per scanline */
    unsigned long scanlinePixels;
    
    /** number of total pixels */
    unsigned long numPixels;
    
    /** iterator for unapdded array */
    CoordinateVectorIter iter;
    
  };


} /* namespace */

#endif /* ARRAY_PADDEDCHANNELTRAVERSER_H */

