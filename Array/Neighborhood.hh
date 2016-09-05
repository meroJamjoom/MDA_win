// ==========================================================================
// $Id: Neighborhood.hh 382 2009-09-25 05:16:27Z heidrich $
// A representation of pixel neighborhoods through offsets and distance
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

#ifndef ARRAY_NEIGHBORHOOD_H
#define ARRAY_NEIGHBORHOOD_H

/*! \file  Neighborhood.hh
    \brief A representation of pixel neighborhoods throgh offsets and distances
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Array.hh"

namespace MDA {

  /** \class Neighborhood Neighborhood.hh

      A representation of pixel neighborhoods through offsets and
      distances that can be iterated over.

      The neighborhood is an N-dimensional box, distances can be
      computed using any p-norm.
  */
  class Neighborhood {
    
  public:
    
    /** \class Neighborhood::Pixel Neighborhood.hh
	Data for a single pixel in a neighborhood */
    class Pixel {
    public:
      /** (signed) array offset to the center pixel in the unpadded array */
      long uOff;
      /** (signed) array offset to the center pixel in the padded array */
      long pOff;
      /** distance from the center pixel according to some metric */
      double dist;
    };
    
    
    /** constructor for a hypercube */
    inline Neighborhood( CoordinateVector &dim, unsigned long padding,
			 unsigned long radius, double norm= 2.0 )
      : current( 0 )
    {
      init( dim, padding, radius, norm );
    }
    
    /** constructor for a hypercube along certain axes */
    inline Neighborhood( CoordinateVector &dim, unsigned long padding,
			 AxisList &axes, unsigned long radius, double norm=2.0 )
      : current( 0 )
    {
      init( dim, padding, axes, radius, norm );
    }
    
    /** constructor for an arbitrary box neighborhood */
    inline Neighborhood( CoordinateVector &dim, unsigned long padding,
			 LongRangeList &box, double norm= 2.0 )
      : current( 0 )
    {
      init( dim, padding, box, norm );
    }
    
    /** destructor */
    inline ~Neighborhood()
    {
      delete [] pixels;
    }
    
    //
    // basic access
    //
    
    /** k-th pixel (bundary checked - use iterator interface for speed) */
    inline const Pixel &operator[]( unsigned long k ) const
    {
      errorCond( k< numPixels, "Invalid pixel index!" );
      return pixels[k];
    }
    
    /** return number of pixels in neighborhood */
    inline unsigned long getNumPixels() const
    {
      return numPixels;
    }
    
    //
    // iterator interface
    //
    
    /** reset iterator */
    inline void begin( unsigned long startPos= 0 )
    {
      errorCond( startPos< numPixels, "Start position out of range!" );
      current= startPos;
    }
    
    /** check termination */
    inline bool isAtEnd() const
    {
      return current>= numPixels;
    }
    
    /** increment iterator */
    inline void operator++()
    {
      current++;
    }
    
    /** return current pixel */
    inline Pixel &operator*()
    {
      return pixels[current];
    }
    
  protected:
    
    /** subclasses can create objects with uninitialized
	neighborhoods, and them themselves */
    inline Neighborhood()
      : current( 0 ), numPixels( 0 )
    {}
    
    /** mask out certain pixels in the box during neighborhood
	construction. True means the pixel is part of the
	neighborhood, false that it gets removed.  This method should
	be overloaded in subclasses for specific behaviors. */
    virtual bool mask( CoordinateVector &pos ) const;
    
    /** actual initialization method used by constructors */
    void init( CoordinateVector &dim, unsigned long padding,
	       unsigned long radius, double p );
      
    /** actual initialization method used by constructors */
    void init( CoordinateVector &dim, unsigned long padding,
	       AxisList &axes, unsigned long radius, double p );
      
    /** actual initialization method used by constructors */
    void init( CoordinateVector &dim, unsigned long padding,
	       LongRangeList &box, double p );
      
    /** number of pixels in the neighborhood */
    unsigned long numPixels;
    
    /** array of pixels in the neighborhood */
    Pixel *pixels;
    
    /** current index for iterator interface */
    unsigned long current;
  };
  

  
  /** \class BallNeighborhood Neighborhood.hh
      
      A neighborhood representing a ball of a certain radius
      (according to any p-norm).
  */
  class BallNeighborhood: public Neighborhood {
    
  public:
    
    /** constructor for a hypercube */
    inline BallNeighborhood( CoordinateVector &dim, unsigned long padding,
			     double _radius, double norm= 2.0 )
      :	Neighborhood(), pNorm( norm ), ballRadius( _radius )
    {
      // initialization has to happen with virtual function table in place
      init( dim, padding, (unsigned long)_radius+1, norm );
    }
    
  protected:
    
    /** mask out regions outside the ball */
    virtual bool mask( CoordinateVector &pos ) const;
    
    /** p-norm to use for distance calculations */
    double pNorm;
    
    /** ball radius */
    double ballRadius;
  };
  
  
  /** \class StarNeighborhood Neighborhood.hh
      
      A neighborhood representing pixels along the major axes
      radiating out from the center.
  */
  class StarNeighborhood: public Neighborhood {
    
  public:
    
    /** constructor for a hypercube */
    inline StarNeighborhood( CoordinateVector &dim, unsigned long padding,
			     unsigned long radius )
      :	Neighborhood()
    {
      // initialization has to happen with virtual function table in place
      init( dim, padding, radius, 1.0 );
    }
    
  protected:
    
    /** mask out regions outside the ball */
    virtual bool mask( CoordinateVector &pos ) const;
  };


} /* namespace */

#endif /* ARRAY_NEIGHBORHOOD_H */

