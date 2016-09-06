// ==========================================================================
// $Id: FloodFill.hh 301 2009-03-25 04:16:15Z heidrich $
// Flood fill of one channel using another as boundary mask
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

#ifndef FILTERS_FLOODFILL_H
#define FILTERS_FLOODFILL_H

/*! \file  FloodFill.hh
    \brief Flood fill of one channel using another as boundary mask
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"

namespace MDA {

  /** \class FloodFill FloodFill.hh
      Flood fill of one channel using another as boundary mask */
  template<class T>
  class FloodFill: public Filter<T> {

  public:
    
    FloodFill( T _boundaryVal= 1.0, T _fillVal= 1.0 )
      : boundaryVal( _boundaryVal ), fillVal( _fillVal )
    {}
    
    /** apply flood fill to a number of dimensions and channels
	\bug the function does support operations on a subset of axes,
	but right now the axis 0 always needs to be part of the axis
	list */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes );
    
  protected:
    
    /**  \class ScanSegment FloodFill.hh
	 a private class of scanline segments (describing jobs in the
	 flood fill algorithm) */
    class ScanSegment {
    public:
      /** constructor */
      ScanSegment( CoordinateVector &_pos, unsigned long _lineIndex,
		   unsigned long _minX, unsigned long _maxX )
	: pos( _pos ), lineIndex( _lineIndex ), minX( _minX ), maxX( _maxX )
      {}
      
      /** position vector for the scanline */
      CoordinateVector pos;
      /** scanline position as an array offset */
      unsigned long lineIndex;
      /** minimum x of the scanline span */
      unsigned long minX;
      /** maximum x of the scanline span */
      unsigned long maxX;
    };
    
    /** flood fill a single channel from a given (set of) seed segment(s) */
    void growSeed( T *channel, T *maskPtr,
		   list<ScanSegment> &seed, CoordinateVector &dim,
		   AxisList &axes, CoordinateVector &incr );
    
    /** the mask value indicating a boundary */
    double boundaryVal;
    
    /** the fill value for the floodfill */
    double fillVal;
  };


} /* namespace */



#endif /* FILTERS_FLOODFILL_H */

