// ==========================================================================
// $Id:$
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

#ifndef DATAANALYSIS_KMEANSCLUSTERING_C
#define DATAANALYSIS_KMEANSCLUSTERING_C

#include <math.h>

#include "MDA/LinearAlgebra/JacobiRotation.hh"
#include "KMeansClustering.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** constructor
    (KMeansClustering owns the SampleVector after this) */
KMeansClustering::KMeansClustering( SampleVector *_data, int _numClusters,
				    int allocClusters )
  : data( _data ), numClusters( _numClusters ), clusterStatsValid( 0 )
{
  unsigned dimension= data->data[0].getSize();
  unsigned long numPoints= data->data.size();
  
  // allocate memory for clusters
  // (initialize means to random values)
  allocClusters= allocClusters < numClusters ? numClusters : allocClusters;
  clusters.reserve( allocClusters );
  ClusterStats s( dimension );
  for( unsigned i= 0 ; i< numClusters ; i++ )
  {
    clusters.push_back( dimension );
    clusters[i].mean.randomize();
  }
  
  // as many cluster IDs as there are data points
  clusterIDs.resize( numPoints );
}

/** destructor */
KMeansClustering::~KMeansClustering()
{
  delete data;
}


/** a fixed number of global iterations */
void
KMeansClustering::globalRelaxation( int numIter, double sampleSize )
{
  while( numIter-- > 0 )
    oneIteration( sampleSize );
}
 

/** globally iterate until convergence */
void
KMeansClustering::globalRelaxation( double maxError, double sampleSize )
{
  do
    oneIteration( sampleSize );
  while( maxChange> maxError );
}

/** split a specific cluster along its largest eigenvector */
void
KMeansClustering::splitCluster( unsigned cluster )
{
  unsigned long i;
  
  // recompute eigen systems if necessary
  if( clusterStatsValid< 1 )
    recomputeClusterStats( 1 );
  
  // split the current mean into two along the largest eigenvector
  unsigned dimension= data->data[0].getSize();
  clusters.push_back( clusters[cluster] ); // copy mean to new cluster
  ClusterStats *c1= &(clusters[cluster]);
  ClusterStats *c2= &(clusters[numClusters]);
  c1->mean.addScalarTimesVector( c1->eigenValues[0], c1->eigenVectors[0] );
  c2->mean.addScalarTimesVector( -c1->eigenValues[0], c1->eigenVectors[0] );
  
  // update the cluster membership only for members of the origonal cluster
  unsigned long numPoints= data->data.size();  
  for( i= 0 ; i< numPoints ; i++ )
    if( clusterIDs[i]== cluster )
      if( dist( data->data[i], c2->mean )< dist( data->data[i], c1->mean ) )
	clusterIDs[i]= numClusters;
  
  // update the two means to actually be the means of the respective
  // point cluster
  c1->mean.zero();
  c2->mean.zero();
  c1->memberCount= c2->memberCount= 0;
  for( i= 0 ; i< numPoints ; i++ )
  {
    if( clusterIDs[i]== cluster )
    {
      c1->mean+= data->data[i];
      c1->memberCount++;
    }
    if( clusterIDs[i]== numClusters )
    {
      c2->mean+= data->data[i];
      c2->memberCount++;
    }
  }
  c1->mean/= (double)c1->memberCount;
  c2->mean/= (double)c2->memberCount;
  
  // update the two convariance matrices;
  Vector h( data->data[0].getSize() );
  c1->cov.zero();
  c2->cov.zero();
  for( i= 0 ; i< numPoints ; i++ )
  {
    if( clusterIDs[i]== cluster )
    {
      h.copy( data->data[i] );
      h-= c1->mean;
      c1->cov.addOuterProduct( h );
    }
    if( clusterIDs[i]== numClusters )
    {
      h.copy( data->data[i] );
      h-= c2->mean;
      c2->cov.addOuterProduct( h );
    }
  }
  c1->cov/= (double)c1->memberCount;
  c2->cov/= (double)c2->memberCount;
  
  // recompute eigen values, eigen vectors
  JacobiRotation solver;
  solver.solve( c1->cov, c1->eigenValues, c1->eigenVectors );
  solver.solve( c2->cov, c2->eigenValues, c2->eigenVectors );
  
  // we have one more cluster now
  numClusters++;
}

/** iteratively split the largest cluster until a certain number of
    clusters has been reached */
void
KMeansClustering::split( unsigned targetClusters )
{
  double maxSpread;
  unsigned largestCluster;
  unsigned i;
  unsigned dim= data->data[0].getSize();
  
  // recompute eigen systems if necessary
  if( clusterStatsValid< 1 )
    recomputeClusterStats( 1);
  
  while( numClusters< targetClusters )
  {
    // find cluster with largest eigenvalue
    maxSpread= -1.0;
    for( i= 0 ; i< numClusters ; i++ )
      if( fabs( clusters[i].eigenValues[0] )> maxSpread )
      {
	maxSpread= fabs( clusters[i].eigenValues[0] );
	largestCluster= i;
      }
    
    // split the largest cluster
    splitCluster( largestCluster );
  }
}


/** merge two specified clusters */
void
KMeansClustering::mergeClusters( unsigned c1, unsigned c2 )
{
  unsigned long i, j;
  
  // make sure c1 is the cluster with the smaller ID
  if( c1> c2 )
  {
    i= c1; c1= c2; c2=i;
  }
  /*
  cerr << "Before:\n";
  for( i= 0 ; i< numClusters ; i++ )
  {
    for( j= 0 ; j< clusters[i].mean.getSize() ; j++ )
      cerr << ' ' << clusters[i].eigenValues[j];
    cerr << "\t\t" << clusters[i].memberCount << endl;
  }
  */
  // new mean is old means weighted by memeber counts
  double w= clusters[c1].memberCount;
  clusters[c1].memberCount+= clusters[c2].memberCount;
  w/= clusters[c1].memberCount;
  clusters[c1].mean*= w;
  clusters[c1].mean.addScalarTimesVector( 1.0-w, clusters[c2].mean );
  
  // update cluster membership, covariance matrix, eigenvalues of c1
  // (only if these quantities are valid for the other clusters)
  if( clusterStatsValid> 0 )
  {
    // covariance matrix
    unsigned long numPoints= data->data.size();
    Vector h( data->data[0].getSize() );
    clusters[c1].cov.zero();
    for( i= 0 ; i< numPoints ; i++ )
    {
      // merge merge all points from c2 into c1
      if( clusterIDs[i]== c2 )
	clusterIDs[i]= c1;
      if( clusterIDs[i]== c1 )
      {
	h.copy( data->data[i] );
	h-= clusters[c1].mean;
	clusters[c1].cov.addOuterProduct( h );
      }
      // the cluster indices > c2 are now reduced by one
      if( clusterIDs[i]> c2 )
	clusterIDs[i]--;
    }
    clusters[c1].cov/= clusters[c1].memberCount;
    
    // eigenvalues and eigenvectors
    JacobiRotation solver;
    solver.solve( clusters[c1].cov, clusters[c1].eigenValues,
		  clusters[c1].eigenVectors );
  }
  
  // consolidate array
  numClusters--;
  for( i= c2 ; i< numClusters ; i++ )
    clusters[i].copy( clusters[i+1] );
  clusters.pop_back();
  /*
  cerr << "After:\n";
  for( i= 0 ; i< numClusters ; i++ )
  {
    for( j= 0 ; j< clusters[i].mean.getSize() ; j++ )
      cerr << ' ' << clusters[i].eigenValues[j];
    cerr << "\t\t" << clusters[i].memberCount << endl;
  }
  */
}


/** interatively merge clusters that are closer than a certain
    multiple of their largest co-variance */
unsigned
KMeansClustering::merge( double closeness )
{
  // recompute cluster stats if necessary
  if( clusterStatsValid< 1 )
    recomputeClusterStats( 1 );
  
  // find clusters that are closer than the extent of one of the clusters
  for( unsigned i= 0 ; i< numClusters ; i++ )
    for( unsigned j= i+1 ; j< numClusters ; j++ )
    {
      // (distance/scale)^2
      double d= dist( clusters[i].mean, clusters[j].mean ) / closeness;
      d*= d;
      // ev is variance (i.e. square of stdev) along eigenvector
      if( d< fabs( clusters[i].eigenValues[0] ) ||
	  d< fabs( clusters[j].eigenValues[0] ) )
	mergeClusters( i, j );
    }
  
  cerr << "Merged to " << numClusters << " clusters\n";
  return numClusters;
}
    

/** recompute cluster membership as well as eigen systems for each
    cluster */
void
KMeansClustering::recomputeClusterStats( unsigned stats )
{
  unsigned long i, j;

  // recompute the eigen systems if required
  if( clusterStatsValid< stats && stats> 0 )
  {
    // update cluster memberships
    unsigned long numPoints= data->data.size();
    for( i= 0 ; i< numClusters ; i++ )
      clusters[i].memberCount= 0;
    for( i= 0 ; i< numPoints ; i++ )
    {
      j= bestCluster( data->data[i] );
      clusterIDs[i]= j;
      clusters[j].memberCount++;
    }
    
    // compute covariance matrices
    Vector h( data->data[0].getSize() );
    for( i= 0 ; i< numClusters ; i++ )
      clusters[i].cov.zero();
    for( i= 0 ; i< numPoints ; i++ )
    {
      h.copy( data->data[i] );
      h-= clusters[clusterIDs[i]].mean;
      clusters[clusterIDs[i]].cov.addOuterProduct( h );
    }
    for( i= 0 ; i< numClusters ; i++ )
      clusters[i].cov/= (double)clusters[i].memberCount;
    
    // compute eigen vectors, eigen values
    JacobiRotation solver;
    for( i= 0 ; i< numClusters ; i++ )
      solver.solve( clusters[i].cov, clusters[i].eigenValues,
		    clusters[i].eigenVectors );
    
    clusterStatsValid= 1;
  }
  
  // if necessary, restore upper triangle of covariance matrices
  // (may have been destroyed by eigensolver)
  if( clusterStatsValid< 2 && stats>= 2 )
  {
    unsigned dim= data->data[0].getSize();
    for( unsigned k= 0 ; k< numClusters ; k++ )
    {
      Matrix &m= clusters[k].cov;
      for( i= 0 ; i< dim ; i++ )
	for( j= i+1 ; j< dim ; j++ )
	  m[j][i]= m[i][j];
    }
    
    clusterStatsValid= 2;
  }
}


/** a single global iteration */
void
KMeansClustering::oneIteration( double sampleSize )
{
  unsigned long i, j;
  // clear tempMeans, counts
  for( i= 0 ; i< numClusters ; i++ )
  {
    clusters[i].tmpMean.zero();
    clusters[i].memberCount= 0;
  }
  
  // perform one iteration
  unsigned long numPoints= data->data.size();  
  
  if( sampleSize<= 0.0 )
    // case 1: use all data points
    for( i= 0 ; i< numPoints ; i++ )
    {
      unsigned cluster= bestCluster( data->data[i] );
      clusters[cluster].tmpMean+= data->data[i];
      clusters[cluster].memberCount++;
    }
  else
    // case 2: use a random subset
    for( i= 0 ; i< sampleSize*numPoints ; i++ )
    {
      j= numPoints*drand48();
      unsigned cluster= bestCluster( data->data[j] );
      clusters[cluster].tmpMean+= data->data[j];
      clusters[cluster].memberCount++;
    }
  
  // update means & max change
  maxChange= 0.0;
  for( i= 0 ; i< numClusters ; i++ )
    if( clusters[i].memberCount> 0 )
    {
      clusters[i].tmpMean/= (double)clusters[i].memberCount;
      double thisChange= dist( clusters[i].tmpMean, clusters[i].mean );
      if( thisChange> maxChange )
	maxChange= thisChange;
      clusters[i].mean.copy( clusters[i].tmpMean );
    }
  /*    else
    {
      // if a mean has an empty cluster, pick a new random pixel value
      unsigned long pt= (unsigned long)(numPoints*drand48());
      clusters[i].mean.copy( data->data[pt] );
      maxChange= 1.0;
      }*/
  maxChange= sqrt( maxChange );
  
  // cluster IDs eigenvectors etc. become invalid after a global iteration
  clusterStatsValid= false;
  
  //cerr << maxChange << endl;
}


/** find the best fitting cluster for a data point */
unsigned
KMeansClustering::bestCluster( const Vector &sample )
{
  double bestDist= dist( sample, clusters[0].mean );
  double cluster= 0;
  for( unsigned i= 1 ; i< numClusters ; i++ )
  {
    double currDist= dist( sample, clusters[i].mean );
    if( currDist< bestDist )
    {
      bestDist= currDist;
      cluster= i;
    }
  }
  return cluster;
}

} /* namespace */

#endif /* DATAANALYSIS_KMEANSCLUSTERING_C */

