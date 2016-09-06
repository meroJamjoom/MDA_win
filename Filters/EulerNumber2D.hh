// ==========================================================================
// $Id: EulerNumber2D.hh 372 2009-09-20 18:57:38Z heidrich $
// Local filter to compute Euler number of a binary image
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

#ifndef FILTERS_EULERNUMBER2D_H
#define FILTERS_EULERNUMBER2D_H

/*! \file  EulerNumber2D.hh
    \brief Local filter to compute Euler number of a binary image
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "HitOrMiss2D.hh"

namespace MDA {

  /** \class EulerNumber2D EulerNumber2D.hh
      Local filter to help compute the Euler number of a binary
      image. The actual Euler number is the integral of the result
      image. */
  template<class T>
  class EulerNumber2D: public HitOrMiss2D<T> {

  public:

    /** default constructor */
    EulerNumber2D( bool eightNeighbor= true )
      : HitOrMiss2D<T>()
    {
      int i, j;
      
      // initialize case table
      for( i= 0 ; i< 512 ; i++ )
      {
	// use top left 2x2 subwindow of the 3x3 HOM pattern
	j= (i&03) | ((i&030)>>1);
	HitOrMiss2D<T>::caseTable[i]=
	  eightNeighbor ? eightNeighborTable[j] : fourNeighborTable[j];
      }
    }
    
  private:
    
    /** table for 4-neighborhood */
    static const double fourNeighborTable[16];
    
    /** table for 8-neighborhood */
    static const double eightNeighborTable[16];

  };


} /* namespace */

#endif /* FILTERS_EULERNUMBER2D_H */

