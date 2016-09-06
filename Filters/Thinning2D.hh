// ==========================================================================
// $Id: Thinning2D.hh 374 2009-09-21 21:13:13Z heidrich $
// A 2D morphological thinning operator
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

#ifndef FILTERS_THINNING2D_H
#define FILTERS_THINNING2D_H

/*! \file  Thinning2D.hh
    \brief A 2D morphological thinning operator
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "HitOrMiss2D.hh"

namespace MDA {

  /** \class Thinning2D Thinning2D.hh
      A 2D morphological thinning operator */
  template<class T>
  class Thinning2D: public Filter<T> {

  public:

    /** constructor */
    Thinning2D();
    
    /** destructor */
    ~Thinning2D()
    {
      for( unsigned i= 0 ; i< 8 ; i++ )
	delete hmt[i];
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );

  private:
    
    /** the 8 Hit-or-Miss Transforms comprising a thinning operation */
    HitOrMiss2D<T> *hmt[8];
    
  };


  /** \class Thickening2D Thinning2D.hh
      A 2D morphological thickening operator (dual of thinning) */
  template<class T>
  class Thickening2D: public Filter<T> {

  public:

    /** constructor */
    Thickening2D();
    
    /** destructor */
    ~Thickening2D()
    {
      for( unsigned i= 0 ; i< 8 ; i++ )
	delete hmt[i];
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );
    
  private:
    
    /** the 8 Hit-or-Miss Transforms comprising a thinning operation */
    HitOrMiss2D<T> *hmt[8];
    
  };


} /* namespace */

#endif /* FILTERS_THINNING2D_H */

