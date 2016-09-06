// ==========================================================================
// $Id:$
// k-means clustering of array data
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

#ifndef DATAANALYSIS_KMEANSCLUSTERING_H
#define DATAANALYSIS_KMEANSCLUSTERING_H

/*! \file  KMeansClustering.hh
    \brief k-means clustering of array data
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/LinearAlgebra/LinAlg.hh"
#include "MDA/Array/Array.hh"
#include "MDA/Base/ChannelList.hh"
#include "SampleVector.hh"

namespace MDA {

  /** \class KMeansClustering KMeansClustering.hh
      k-means clustering of array data */
  class KMeansClustering {

  public:

    /** \class ClusterStats KMeansClustering
	data structure for holding cluster statistics (mean,
	eigenvalues etc.) */
    class ClusterStats {

    public:
      
      /** constructor (empty stats) */
      inline ClusterStats( unsigned dim )
	: mean( dim ), cov( dim, dim ), eigenVectors( dim, dim ),
	  eigenValues( dim ), tmpMean( dim )
      {}
      
      /** copy constructor only copies mean, not rest of stats */
      inline ClusterStats( const ClusterStats &s )
      	: mean( s.mean.getSize() ), cov( s.mean.getSize(), s.mean.getSize() ),
	  eigenVectors( s.mean.getSize(), s.mean.getSize() ),
	  eigenValues( s.mean.getSize() ), tmpMean( s.mean.getSize() )
      {
	mean.copy( s.mean );
      }
      
      /** copy full content */
      inline void copy( const ClusterStats &s, unsigned stats= 1 )
      {
	mean.copy( s.mean );
	if( stats )
	{
	  memberCount= s.memberCount;
	  cov.copy( s.cov );
	  eigenValues.copy( s.eigenValues );
	  eigenVectors.copy( s.eigenVectors );
	}
      }
      
      /** mean value */
      Vector mean;
      
      /** number of data values associated with this cluster */
      unsigned memberCount;
      
      /** covariance matrix */
      Matrix cov;
      
      /** eigen vectors */
      Matrix eigenVectors;
      
      /** eigen values */
      Vector eigenValues;

      /** temp vector for computing mean */
      Vector tmpMean;
    };
    
    
    /** constructor */
    KMeansClustering( SampleVector *data, int numClusters,
		      int allocClusters= -1 );
    
    /** destructor */
    ~KMeansClustering();
    
    /** a fixed number of global iterations */
    void globalRelaxation( int numIter, double sampleSize= -1.0 );
    
    /** globally iterate until convergence */
    void globalRelaxation( double maxError, double sampleSize= -1.0 );
      
    /** split a specific cluster along its largest eigenvector */
    void splitCluster( unsigned cluster );
    
    /** iteratively split the largest cluster until a certain number of
	clusters has been reached */
    void split( unsigned targetClusters );

    /** merge two specified clusters */
    void mergeClusters( unsigned c1, unsigned c2 );
    
    /** interatively merge clusters that are closer than a certain
	multiple of their largest co-variance
	\return new number of clusters */
    unsigned merge( double closeness= 1.0 );
    
    /** return a vector of the clusters
	\param stats: which stats are required (0: mean only, 1:
	mean&eigen system, 2: mean, eigen system, and covariance matrix) */
    inline const vector<ClusterStats> *getClusters( unsigned stats= 0 )
    {
      if( clusterStatsValid< stats )
	recomputeClusterStats( stats );
      return &clusters;
    }
    
  protected:
    
    /** find the best fitting cluster for a data point */
    unsigned bestCluster( const Vector &sample );
    
    /** recompute cluster membership as well as eigen systems for each
	cluster */
    void recomputeClusterStats( unsigned stats );
    
    /** a single Lloyd iteration */
    void oneIteration( double sampleSize );
    
    /** the maximum change during the last iteration */
    double maxChange;
    
    /** number of clusters */
    int numClusters;
    
    /** actual data vector */
    SampleVector *data;
    
    /** whether or not the cluster IDs, EVs etc are valid */
    int clusterStatsValid;
    
    /** cluster statistics */
    vector<ClusterStats> clusters;
    
    /** cluster memberships */
    vector<unsigned> clusterIDs;
    
  private:
    
    /** temporary means */
    SampleVector *tmpMeans;
  };

} /* namespace */

#endif /* DATAANALYSIS_KMEANSCLUSTERING_H */

