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

#ifndef CONNCOMPPROPS_CLUSTERMODIFIER_C
#define CONNCOMPPROPS_CLUSTERMODIFIER_C

#include "ClusterModifier.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** default constructor */

ClusterModifier::ClusterModifier()
{
}

/** clusters a number of samplepoints with a cutoff distance criterion. returns the number of found clusters  */
long ClusterModifier::clusterCutoff( double cutOffDistance, SampleVector *samples, vector<long> *clusterIndex)
{
    // Early exit if assertions not met
    if( samples == NULL || clusterIndex == NULL || (samples->data).size() == 0
			|| clusterIndex->size() != (samples->data).size() )
	return 0;
    
    long dimension = ( (samples->data).at( 0 ) ).getSize();
    
    // Initialize clusterindex
    for ( unsigned long i = 0; i < clusterIndex->size(); i++ )
	  clusterIndex->at( i ) = -1;

    long nextClusterNum = 0;
    unsigned long numSamples = (samples->data).size(); 
    for( unsigned long currSample = 0; currSample < numSamples; currSample++ )
    {
	//Assign new Clusternumber if not set
	if( clusterIndex->at( currSample ) == -1 )
	{
	    clusterIndex->at( currSample ) = nextClusterNum;
	    nextClusterNum++;
	}

	//Assign current clusternumber to all neighbors
	Vector currPoint(dimension, 0.0);
	if( ((samples->data).at( currSample ) ).getSize() != dimension )
	{
	    cerr << "Error: Dimensions must match" << endl;
	    return 0;
	}
	currPoint.assign( (samples->data).at( currSample ) );
	Vector neighPoint(dimension, 0.0);
	for( unsigned long neighbor = 0; neighbor < numSamples; neighbor++ )
	{
	    // CurrSample itself is not a neighbor of CurrSample
	    if( neighbor == currSample ) 
		continue;

	    // Continue if already in cluster
	    if( clusterIndex->at( neighbor ) == clusterIndex->at( currSample ) )
		continue;

	    // Assign current Clusternumber if in range
	    if( ((samples->data).at( neighbor ) ).getSize() != dimension )
	    {
		cerr << "Error: Dimensions must match" << endl;
		return 0;
	    }
	    neighPoint.assign( (samples->data).at( neighbor ) );
	    if( dist( currPoint, neighPoint ) <= cutOffDistance )
		clusterIndex->at( neighbor ) = clusterIndex->at( currSample );
	}
	
    }

    return nextClusterNum;
}

/** computes the mean of all clusters and inputs new samplepoints for means. Mapping from samples to result is done via clusterIndex. */
void ClusterModifier::computeClusterMeans( SampleVector* samples, vector<long> *clusterIndex, SampleVector* result )
{
    // Early exit if assertions not met
    if( samples == NULL || clusterIndex == NULL || result == NULL || (samples->data).size() == 0 
			|| clusterIndex->size() != (samples->data).size() )
	return;

    long dimension = ( (samples->data).at( 0 ) ).getSize();
    (result->data).clear();

    // Find number of clusters
    unsigned long clusterNum = 0;
    for( unsigned long i = 0; i < clusterIndex->size(); i++ )
	clusterNum = ( clusterIndex->at( i ) > clusterNum ) ? clusterIndex->at( i ) : clusterNum ;

    // Iterate over all clusters
    Vector mean;
    Vector zeroVec( dimension, 0.0 );
    unsigned long clusterSize;
    for( unsigned long currCluster = 0; currCluster <= clusterNum; currCluster++ )
    {
	mean.copy( zeroVec );
	clusterSize = 0;
	for( unsigned long i = 0; i < clusterIndex->size(); i++ )
	{
	    if( clusterIndex->at( i ) == currCluster )
	    {
		mean += (samples->data).at( i );
		clusterSize++;
	    }
	}
	if( clusterSize > 0 )
	    mean /= (double) clusterSize;
	
	(result->data).push_back( mean );
    }

    return;
}



} /* namespace */

#endif /* CONNCOMPPROPS_CLUSTERMODIFIER_C */

