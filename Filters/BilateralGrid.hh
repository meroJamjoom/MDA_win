// ==========================================================================
// $Id: BilateralGrid.hh 999 2014-05-28 15:07:31Z heidrich $
// A fast bilateral filter approximation
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

#ifndef FILTERS_BILATERALGRID_H
#define FILTERS_BILATERALGRID_H

/*! \file  BilateralGrid.hh
    \brief A fast bilateral filter approximation
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Filter.hh"


namespace MDA {

  /** \class BilateralGrid BilateralGrid.hh
      A fast bilateral filter approximation */
  template<class T>
  class BilateralGrid: public Filter<T> {

  public:

    /** constructor
        \param _sigma: standard deviation
        \param _edgeStopSigma: std. dev. of edge stopping function
        \param resize: ratio of grid spacing : sigma */
    BilateralGrid( double _sigma, double _edgeStopSigma= .1, double resize= 2 )
      : sigma( _sigma ), edgeStopSigma( _edgeStopSigma ),
	spatialScale( _sigma> resize ? _sigma/resize : 1.0 ),
	intensityScale( _edgeStopSigma/resize ),
	grid( NULL )
    {}
    
    /** destructor */
    inline ~BilateralGrid()
    {
      clear();
    }
    
    /** clean up the grid */
    inline void clear()
    {
      if( grid!= NULL )
	delete grid;
      grid= NULL;
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
                        ChannelList &channels, AxisList &axes )
    {
      construct( a, channels, channels );
      process( boundary, axes, sigma/spatialScale,
	       edgeStopSigma/intensityScale );
      slice( a, channels, channels );
      clear();
      return true;
    }
    
    /** "construct" a new grid */
    virtual bool construct( Array<T> &a, ChannelList &dataChannels,
		    ChannelList &edgeStopChannels );
    
    /** "process" an existing grid */
    bool process( BoundaryMethod boundary, AxisList &axes,
		  double sigma, double edgeStopSigma );
    
    /** "slice" an existing grid */
    bool slice( Array<T> &a, ChannelList &dataChannels,
		ChannelList &edgeStopChannels );
    
  protected:
    
    /** spatial sigma */
    double sigma;

    /** sigma in intensity dimesnions */
    double edgeStopSigma;
    
    /** spatial scaling factor */
    double spatialScale;
    
    /** intensity scaling factor */
    double intensityScale;
    
    /** the actual grid */
    Array<T> *grid;
    
    /** minimum value per channel */
    vector<T> minVal;

    /** maximum value per channel */
    vector<T> maxVal;
  };



  /** \class BilateralGridWeighted BilateralGrid.hh
      A fast bilateral filter approximation that allows an additional weight*/
  template<class T>
  class BilateralGridWeighted: public BilateralGrid<T> {

  public:

    /** constructor
        \param _sigma: standard deviation
        \param _edgeStopSigma: std. dev. of edge stopping function
        \param resize: ratio of grid spacing : sigma */
    BilateralGridWeighted( double _sigma, double _edgeStopSigma= .1, double resize= 2 )
      : BilateralGrid<T>(_sigma, _edgeStopSigma, resize)
    {}
    
    /** destructor */
    inline ~BilateralGridWeighted()
    {
      // pretty sure base class destructor is called automatically,
      // but to be safe I'll do the same thing here that happens above.
      BilateralGrid<T>::clear();
    }

    /** "construct" a new grid */
    virtual bool construct( Array<T> &a, ChannelList &dataChannels,
		    ChannelList &edgeStopChannels );

  };

} /* namespace */

#endif /* FILTERS_BILATERALGRID_H */

