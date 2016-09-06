// ==========================================================================
// $Id: FilterFactory.hh 252 2008-11-28 06:23:24Z heidrich $
// a factory for filters
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

#ifndef FILTERS_FILTERFACTORY_H
#define FILTERS_FILTERFACTORY_H

/*! \file  FilterFactory.hh
    \brief a factory for filters
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"
#include "FilterOption.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"


namespace MDA {

  /** \class FilterParams FilterFactory.hh
      a collection of parameters required to create any of the filters
      supported by FilterFactory */
  struct FilterParams {
    int radius; /** radius of a filter, in pixels */
    double sigma; /** standard deviation for filters involving Gaussians */
    double edgeStopSigma; /** std. dev. of edge stopping function */
    EXPR::ExpressionSequence values; /** filter values (SeparableFilter etc.) */
    bool normalize; /** whether or not to normalize the filter */
  };
  
    
  /** \class FilterFactory FilterFactory.hh
      a factory for filters */
  template<class T>
  class FilterFactory {

  public:
    
    /** default constructor */
    inline FilterFactory( FilterParams *params= NULL )
    {
      if( params== NULL )
      {
	parameters= new FilterParams;
	ownParameters= true;
      }
      else
      {
	parameters= params;
	ownParameters= false;
      }
    }
    
    /** destructor */
    inline ~FilterFactory()
    {
      if( ownParameters )
	delete parameters;
    }
    
    /** report filter praemeters */
    inline FilterParams *getFilterParameters()
    {
      return parameters;
    }
    
    /** create a new filter of filterType using previously set parameters
     *
     * \todo the MDA::FilterFactory does not currently expose all
     * options of the MDA::HitOrMiss2D operator (symmetry of patterns,
     * preservation of misses, specific hit and miss values)
     */
    Filter<T> *create( FilterType filterType );
    
    
  protected:
    
    /** filter radius in pixels */
    int radius;
    
    /** sigma for Gauss-type filters */
    double sigma;
    
    /** filter parameters */
    FilterParams *parameters;
    
    /** whether we own the filter parameters, and are responsible for
	deleting them */
    bool ownParameters;
  };


} /* namespace */

#endif /* FILTERS_FILTERFACTORY_H */

