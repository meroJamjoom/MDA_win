// ==========================================================================
// $Id: SampleVector.hh 257 2008-12-30 01:00:39Z heidrich $
// an array of n-dimensional samples
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

#ifndef DATAANALYSIS_SAMPLEVECTOR_H
#define DATAANALYSIS_SAMPLEVECTOR_H

/*! \file  SampleVector.hh
    \brief an array of n-dimensional samples
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>
#include "MDA/Base/ChannelList.hh"
#include "MDA/LinearAlgebra/LinAlg.hh"
#include "MDA/Array/Array.hh"

namespace MDA {

  using namespace std;
  
  /** \class KMeansClustering KMeansClustering.hh
      k-means clustering of array data */
  class SampleVector {
    
  public:
    
    /** default constructor */
    SampleVector( unsigned long numElements= 0 )
    {
      if( numElements> 0 )
	data.reserve( numElements );
      data.clear();
    }
    
    /** constructor from float Arrays */
    SampleVector( Array<float> &a, ChannelList &ch );

    /** constructor from double Arrays */
    SampleVector( Array<double> &a, ChannelList &ch );
    
    /** the actual data vector */
    vector<Vector> data;
    
  };
  
} /* namespace */

#endif /* DATAANALYSIS_SAMPLEVECTOR_H */

