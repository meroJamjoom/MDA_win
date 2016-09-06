// ==========================================================================
// $Id: BilateralFilter.hh 248 2008-11-25 05:45:35Z heidrich $
// Brute-force bilateral filtering
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

#ifndef FILTERS_BILATERALFILTER_H
#define FILTERS_BILATERALFILTER_H

/*! \file  BilateralFilter.hh
    \brief Brute-force bilateral filtering
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"

namespace MDA {

  /** \class BilateralFilter BilateralFilter.hh
      Brute-force bilateral filtering */
  template<class T>
  class BilateralFilter: public Filter<T> {

  public:

    /** constructor
        \param _sigma: standard deviation
	\param _edgeStopSigma: std. dev. of edge stopping function
        \param rad: radius of filter (if -1: use 2*sigma) */
    BilateralFilter( double _sigma, double _edgeStopSigma= .1, int rad= -1 )
      : sigma( _sigma ), edgeStopSigma( _edgeStopSigma ), radius( rad )
    {}

    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes );
    
  protected:

    /** radius of the median */
    int radius;
    
    /** spatial standard deviation */
    double sigma;
    
    /** photometric / edge stop standard deviation */
    double edgeStopSigma;
    
  };


} /* namespace */

#endif /* FILTERS_BILATERALFILTER_H */

