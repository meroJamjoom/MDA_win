// ==========================================================================
// $Id: Prune2D.hh 374 2009-09-21 21:13:13Z heidrich $
// A 2D morphological prune filter.
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

#ifndef FILTERS_PRUNE2D_H
#define FILTERS_PRUNE2D_H

/*! \file  Prune2D.hh
    \brief A 2D morphological prune filter.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "HitOrMiss2D.hh"

namespace MDA {
  
  // forward declaration
  template<class T> class Endpoints2D;
  
  
  /** \class Prune2D Prune2D.hh
      A 2D morphological prune filter. */
  template<class T>
  class Prune2D: public HitOrMiss2D<T> {

  public:

    /** default constructor */
    Prune2D()
      : HitOrMiss2D<T>( se1, ma1, true, true, true, 0.0, 1.0 )
    {
      HitOrMiss2D<T> second( se2, ma2, true, true, true, 0.0, 1.0 );
      (*this)&= second;
    }

  private:
    
    /** structuring element 1 */
    static bool se1[9];
    
    /** mask 1 */
    static bool ma1[9];
    
    /** structuring element 2 */
    static bool se2[9];
    
    /** mask 2 */
    static bool ma2[9];
    
    /** Endpoints2D uses the same structuring elements and masks */
    friend class Endpoints2D<T>;
    
  };
  
  
  /** \class Endpoints2D Prune2D.hh
      A 2D morphological endpoint filter (selects exactly those points
      eliminated by the pruning filter). */
  template<class T>
  class Endpoints2D: public HitOrMiss2D<T> {

  public:

    /** default constructor */
    Endpoints2D()
      : HitOrMiss2D<T>( Prune2D<T>::se1, Prune2D<T>::ma1, true, true )
    {
      HitOrMiss2D<T> second( Prune2D<T>::se2, Prune2D<T>::ma2, true, true );
      (*this)|= second;
    }
    
  };

} /* namespace */

#endif /* FILTERS_PRUNE2D_H */

