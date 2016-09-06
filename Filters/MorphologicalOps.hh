// ==========================================================================
// $Id: MorphologicalOps.hh 247 2008-11-25 05:23:41Z heidrich $
// Morphological operators (erode, dilate, open, close etc.)
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

#ifndef FILTERS_MORPHOLOGICALOPS_H
#define FILTERS_MORPHOLOGICALOPS_H

/*! \file  MorphologicalOps.hh
    \brief Morphological operators (erode, dilate, open, close etc.)
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"
#include "Filter1D.hh"

namespace MDA {

  /** \class Erode1D MorphologicalOps.hh
      1D erosion filter */
  template<class T>
  class Erode1D: public Filter1D<T> {

  public:

    /** default constructor */
    Erode1D( unsigned _radius= 1 ):
      radius( _radius )
    {}

    /** apply filter to one channel of an array (pure virtual) */
    virtual bool apply( Array<T> &array, BoundaryMethod boundary,
                        unsigned int axis, unsigned int inChannel,
                        unsigned int outChannel );
    
  protected:
    
    /** apply filter to a single line in the array */
    virtual void apply( T *startPos, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut );
    
    /** filter radius */
    unsigned int radius;
    
  };


  /** \class Dilate1D MorphologicalOps.hh
      1D dilation filter */
  template<class T>
  class Dilate1D: public Filter1D<T> {

  public:

    /** default constructor */
    Dilate1D( unsigned _radius= 1 ):
      radius( _radius )
    {}
    
    /** apply filter to one channel of an array (pure virtual) */
    virtual bool apply( Array<T> &array, BoundaryMethod boundary,
                        unsigned int axis, unsigned int inChannel,
                        unsigned int outChannel );
    
  protected:
    
    /** apply filter to a single line in the array */
    virtual void apply( T *startPos, unsigned long incr,
			unsigned long numElements, BoundaryMethod boundary,
			T background, T *startPosOut );
    
    /** filter radius */
    unsigned int radius;
    
  };


} /* namespace */

#endif /* FILTERS_MORPHOLOGICALOPS_H */

