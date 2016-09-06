// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: Felix Heide
// Email:   heide@informatik.uni-siegen.de
// ==========================================================================

#ifndef CONNCOMPPROPS_CLUSTERMODIFIER_H
#define CONNCOMPPROPS_CLUSTERMODIFIER_H

/*! \file  ClusterModifier.hh
    \brief A simple class for clustering, modifying n-dimensional clusters.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/Range.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "SampleVector.hh"

#include <math.h>
#include <vector>
#include <limits.h>

namespace MDA {

  using namespace std;

  /** \class ClusterModifier ClusterModifier.hh
      A simple class for clustering, modifying n-dimensional clusters. */
  
  class ClusterModifier {

  public:

    /** default constructor */
    ClusterModifier();

    /** clusters a number of samplepoints with a cutoff distance criterion. returns the number of found clusters  */
    long  clusterCutoff( double cutOffDistance, SampleVector *samples, vector<long> *clusterIndex);

    /** computes the mean of all clusters and inputs new samplepoints for means. Mapping from samples to result is done via clusterIndex. */
    void computeClusterMeans( SampleVector* samples, vector<long> *clusterIndex, SampleVector* result );

  };


} /* namespace */

#endif /* CONNCOMPPROPS_CLUSTERMODIFIER_H */

