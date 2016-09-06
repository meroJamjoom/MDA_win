// ==========================================================================
// $Id: ImprovedHarrisCorner.hh 319 2009-05-27 21:17:01Z heidrich $
// A corner detector based on Harris', but with better localization
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

#ifndef FILTERS_IMPROVEDHARRISCORNER_H
#define FILTERS_IMPROVEDHARRISCORNER_H

/*! \file  ImprovedHarrisCorner.hh
    \brief A corner detector based on Harris', but with better localization
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"

namespace MDA {

  /** \class ImprovedHarrisCorner ImprovedHarrisCorner.hh
      A corner detector based on Harris', but with better
      localization. The corners are local maxima in the result
      channel */
  template<class T>
  class ImprovedHarrisCorner: public Filter<T> {

  public:

    /** constructor */
    ImprovedHarrisCorner( double _sigma= 1.0, double _traceWeight= 0.05,
			  double _cornerThreshold= 0.01 )
      : sigma( _sigma ), traceWeight( _traceWeight ),
	cornerThreshold( _cornerThreshold )
    {}
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes );
    
  protected:
    
    /** sigma for a Gaussian smoothing */
    double sigma;
    
    /** how to weight the trace term in the Harris corner detector */
    double traceWeight;
    
    /** threshold of original Harris value below which the
	"cornerness" is set to 0 */
    double cornerThreshold;
    
  };


} /* namespace */

#endif /* FILTERS_IMPROVEDHARRISCORNER_H */

