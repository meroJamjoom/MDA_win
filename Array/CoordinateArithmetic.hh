// ==========================================================================
// $Id: CoordinateArithmetic.hh 953 2012-05-08 21:20:16Z krim $
// helper functions for common calculations using array dimensions
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

#ifndef ARRAY_COORDINATEARITHMETIC_H
#define ARRAY_COORDINATEARITHMETIC_H

/*! \file  CoordinateArithmetic.hh
    \brief helper functions for common calculations using array dimensions
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <MDA/Base/CoordinateVector.hh>

namespace MDA {

  //
  // helpers for coordinate vector arithmetic 
  //
  
  /** number of pixels in a channel of a certain dimension */
  inline unsigned long size( const CoordinateVector &dim )
  {
    unsigned long size= 1;
    for( unsigned i= dim.vec.size() ; i> 0 ; )
      size*= dim.vec[--i];
    return size;
  }
  
  /** number of scanlines in a channel of a certain dimension */
  inline unsigned long scanlines( const CoordinateVector &dim )
  {
    unsigned long size= 1;
    for( unsigned i= dim.vec.size() ; i> 1 ; )
      size*= dim.vec[--i];
    return size;
  }
  
  /** compute the dimensions of a padded channel */
  inline CoordinateVector pad( const CoordinateVector &dim, unsigned long padding )
  {
    CoordinateVector padded= dim;
    for( unsigned i= dim.vec.size() ; i> 0 ; )
      padded.vec[--i]+= 2*padding;
    return padded;
  }
  
  /** compute offset for a position with respect to an (unpadded)
      channel of given dimensions */
  inline unsigned long pixelOffset( const CoordinateVector &pos,
				    const CoordinateVector &dim )
  {
    unsigned long off= 0;
    for( unsigned i= pos.vec.size() ; i> 0 ; )
    {
      --i;
      off= off*dim.vec[i] + pos.vec[i];
    }
    return off;
  }
  
  /** compute offset for a position with respect to a padded channel
      of given dimensions */
  inline unsigned long pixelOffset( const CoordinateVector &pos,
				    const CoordinateVector &paddedDim,
				    unsigned long padding )
  {
    unsigned long off= 0;
    for( unsigned i= pos.vec.size() ; i> 0 ; )
    {
      --i;
      off= off*paddedDim.vec[i] + pos.vec[i]+padding;
    }
    return off;
  }

} /* namespace */

#endif /* ARRAY_COORDINATEARITHMETIC_H */

