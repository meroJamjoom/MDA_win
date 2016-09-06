// ==========================================================================
// $Id: ConnectedComponent.hh 373 2009-09-20 23:08:07Z heidrich $
// labels&counts connected components of a specific pixel value
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

#ifndef FILTERS_CONNECTEDCOMPONENT_H
#define FILTERS_CONNECTEDCOMPONENT_H

/*! \file  ConnectedComponent.hh
    \brief labels&counts connected components of a specific pixel value
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "FloodFill.hh"

namespace MDA {

  /** \class ConnectedComponent ConnectedComponent.hh

      labels&counts connected components of pixels with value 0.0 < x
      < 1.0 over a background of pixels with value 0.0. The resulting
      pixel values are integers, with 0 representing the background
      and 1...N representing membership to the Nth connected
      component. */
  template<class T>
  class ConnectedComponent: public FloodFill<T> {
    
  public:
    
    /** constructor */
    ConnectedComponent()
      : FloodFill<T>( 0.0, 1.0 ), numComp( 0 )
    {}
    
    /** find connected components in a number of channels
	\bug currently, the axis list is ignored (components are
	always connected along all dimensions) */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes );
    
    /** return number of components found during last filter
	application (if multiple channels were processed in one call,
	only the components in the last channel are reported) */
    inline unsigned long getNumComp() const
    {
      return numComp;
    }
    
  protected:
    
    /** number of components form last run */
    unsigned long numComp;
    
  };


} /* namespace */

#endif /* FILTERS_CONNECTEDCOMPONENT_H */

