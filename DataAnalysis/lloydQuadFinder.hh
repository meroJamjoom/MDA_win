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

#ifndef CONNCOMPPROPS_LLOYDQUADFINDER_H
#define CONNCOMPPROPS_LLOYDQUADFINDER_H

/*! \file  lloydQuadFinder.hh
    \brief A quadrilaterial line fitter for a set of 2D samplepoints.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "SampleVector.hh"
#include "MDA/GeometricPrimitives/Line.hh"
#include "MDA/GeometricPrimitives/Point.hh"
#include "GeometricPrimitiveFitting.hh"
#include "KMeansClustering.hh"

namespace MDA {

  using namespace std;

  /** \class lloydQuadFinder lloydQuadFinder.hh
      A quadrilaterial line fitter for a set of 2D samplepoints. */
  
  class lloydQuadFinder {

  // comparison function for the indexed sorting
  static bool compareIndexPairs ( Vector i, Vector j);

  public:

  /** default constructor */
  lloydQuadFinder();

  /** Does multiple lloyd iterations with a set of lines on a set of samples. Stops if converged ( no samples
      switched clusters in last iteration ) or if degenerated. */
  unsigned long lloydLines( SampleVector* samples, vector<Line*> *l, bool* degenerated);

  /** Finds quadrilateral corners of a set of samples. The samples have to be
      ordered either clockwise or counterclockwise */
  void computeCorners( SampleVector *samples, SampleVector* corners ); 

  protected:

  /** Does one lloyd iteration with a set of lines on a set of samples. ClusterIndex is a vector of 
      indices of the clusters(the line) every samplevector was in (is nearest to) in the last iteration step.
      For first iteration these values can all be set to -1. Degenerates if to at least one line no sample could
      be assigned. Returns the number of samples, which switched clusters */
  long lloydLineIter( vector<long>* clusterIndex, SampleVector* samples, vector<Line*> *l, bool* degenerated);

  /** Intersects two lines */
  void interSect( const Line &l1, const Line &l2, Vector& intersection );

  };


} /* namespace */

#endif /* CONNCOMPPROPS_LLOYDQUADFINDER_H */

