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

#ifndef CONNCOMPPROPS_LLOYDQUADFINDER_C
#define CONNCOMPPROPS_LLOYDQUADFINDER_C

#include "lloydQuadFinder.hh"
#include <algorithm>

#define MAXITERATIONS 100

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  // comparison function for the indexed sorting
  bool lloydQuadFinder::compareIndexPairs ( Vector i, Vector j)
  { 
      return ( i[0] < j[0] ); 
  }

  /** default constructor */
  lloydQuadFinder::lloydQuadFinder(){}

  /** Does one lloyd iteration with a set of lines on a set of samples. ClusterIndex is a vector of 
      indices of the clusters(the line) every samplevector was in (is nearest to) in the last iteration step.
      For first iteration these values can all be set to -1. Degenerates if to at least one line no sample could
      be assigned. Returns the number of samples, which switched clusters */
  long lloydQuadFinder::lloydLineIter( vector<long>* clusterIndex, SampleVector* samples, vector<Line*> *l, bool* degenerated)
  {
      // Test for sanity
      if( clusterIndex == NULL || l == NULL || samples == NULL ||
	  clusterIndex->size() != (samples->data).size() )
	  return -1;

      vector<SampleVector> newClusters(l->size());

      vector<bool> degenerate( l->size(), true );

      long switchedSamples = 0;
      Vector currPoint( 2, 0.0 );
      for( unsigned long i = 0; i < (samples->data).size(); i++)
      {
	  currPoint.copy( (samples->data).at(i) );
    
	  //Assign current Point to nearest line
	  int lowIndex = 0;
	  double currDist;
	  double lowDist = dist( *(l->at(0)), currPoint );
	  for( int j = 0; j < l->size(); j++ )
	  {
	      currDist = dist( *(l->at(j)), currPoint );
	      if( currDist < lowDist )
	      {
		  lowIndex = j;
		  lowDist = currDist;
	      }
	  }

	  //Check if cluster switched
	  if( lowIndex != clusterIndex->at(i) )
	      switchedSamples++;

	  //Update new cluster Index
	  clusterIndex->at(i) = lowIndex;
	  newClusters.at(lowIndex).data.push_back( currPoint );
	  degenerate.at(lowIndex) = false;
      }

      for( int i = 0; i < l->size(); i++ )
	  fitToDataPoints( *(l->at(i)), newClusters.at(i) );

      //Check for degeneration ( no points assigned to at least one line )
      *degenerated = false;
      for( unsigned long i = 0; i < degenerate.size(); i++)
	  if( degenerate.at(i) )	
	      *degenerated = true;

      return switchedSamples;
  }


  /** Does multiple lloyd iterations with a set of lines on a set of samples. Stops if converged ( no samples
      switched clusters in last iteration ) or if degenerated. */
  unsigned long lloydQuadFinder::lloydLines( SampleVector* samples, vector<Line*> *l, bool* degenerated)
  {
      vector<long>* clusterIndex = new vector<long>( (samples->data).size(), -1 );

      bool degenerate = false;

      //Do iterations until points do not longer switch clusters (or upper limit of iterations achieved )
      unsigned long iterations;
      for( iterations  = 0; iterations < MAXITERATIONS; iterations++ )
      {
	  // Do an iteration for the first cluster if not finished
	  if( lloydLineIter( clusterIndex, samples, l, &degenerate ) == 0 )
	      break;

	  if( degenerate )
	      break;
      }

      //Check for degeneration
      if( iterations >= MAXITERATIONS || degenerate )
	  *degenerated = true;
      else	
	  *degenerated = false;

      return iterations;
  }

  /** Intersects two lines */
  void lloydQuadFinder::interSect( const Line &l1, const Line &l2, Vector& intersection )
  {
      double u = (l2.p2[0] - l2.p1[0])*(l1.p1[1] - l2.p1[1]) - (l2.p2[1] - l2.p1[1])*(l1.p1[0] - l2.p1[0]);
      u /= (l2.p2[1] - l2.p1[1])*(l1.p2[0] - l1.p1[0]) - (l2.p2[0] - l2.p1[0])*(l1.p2[1] - l1.p1[1]);

      if( u == 0.0 )
      {
	  intersection.set(0, 0.0);
	  intersection.set(1, 0.0);
      }else
      {
	  intersection.set( 0, l1.p1[0] + u * ( l1.p2[0] - l1.p1[0] ) );
	  intersection.set( 1, l1.p1[1] + u * ( l1.p2[1] - l1.p1[1] ) );
      }
      
      return;
  }

  /** Finds quadrilateral corners of a set of samples. The samples have to be
      ordered either clockwise or counterclockwise */
  void lloydQuadFinder::computeCorners( SampleVector *samples, SampleVector* corners)
  {
      //Exit if wrong size
      if( samples == NULL || (samples->data).size() == 0 )
	  return;

      // Clear result vector
      (corners->data).clear();

      //Compute derivation of ordered sample sequence
      SampleVector* deriv = new SampleVector();
      Vector derivVec(2, 0.0);
      for( long i = 0; i < (samples->data).size(); i++ )
      {
	  if( (samples->data).at((i+1) % (samples->data).size() ).getSize() != 2 
	    || (samples->data).at((i-1 + (samples->data).size()) % (samples->data).size() ).getSize() != 2 )
	  {
	      cerr << "Supports only 2D sample points." << endl;
	      return;
	  }else
	  {
	      derivVec.copy( (samples->data).at( (i+1) % (samples->data).size() ) );
	      derivVec -=  (samples->data).at( (i - 1 + (samples->data).size()) % (samples->data).size() );
	      derivVec /= 2.0;
	      (deriv->data).push_back( derivVec );
	  }
      }

      //Compute the samplesize factor so that at least 4 random samples are used
      double sampleSize =  5.0 /(double)( (deriv->data).size() );

      //Do a 4-Means Clustering on the derivation data to find the 4 main line orientations.
      KMeansClustering kmeansSolver(deriv, 4);
      kmeansSolver.globalRelaxation( 0.001, sampleSize  );
      const vector<KMeansClustering::ClusterStats>* clusters = kmeansSolver.getClusters(2);

      // Return if not exactly 4 clusters were found
      if( clusters->size() < 4 )
	  return;

      //Debug Output
      //cerr << "Clusters found: " << clusters->size() << endl;

      //Extract the samplepoints of the 4 main line orientations:
      vector<SampleVector> initialClusters(clusters->size());
      Vector currPoint(2, 0.0);
      for( unsigned long i = 0; i < (samples->data).size(); i++ )
      {
	  currPoint.copy( (samples->data).at(i) );
    
	  //Assign current Point to nearest cluster
	  int lowIndex = 0;
	  double currDist;
	  double lowDist = dist( (clusters->at(0)).mean , (deriv->data).at(i)  );
	  for( int j = 1; j < clusters->size(); j++ )
	  {
	      currDist = dist( (clusters->at(j)).mean , (deriv->data).at(i)  );
	      if( currDist < lowDist )
	      {
		  lowIndex = j;
		  lowDist = currDist;
	      }
	  }

	  //Update new cluster means
	  initialClusters.at(lowIndex).data.push_back( currPoint );
      }

      //Fit lines to resulting initial clusters
      //Copy lines to restore if lloyd degenerates
      vector<Line*> l;
      vector<Line*> lCopy;
      for( int i = 0; i < clusters->size(); i++ )
      {
	  if( initialClusters.at(i).data.size() == 0 )
	      return;
	  
	  l.push_back( new Line(2) );
	  fitToDataPoints( *(l.back()), initialClusters.at(i) );
	  lCopy.push_back( new Line( l.back()->p1, l.back()->p2 ) );
      }

      //Finally do lloyd iterations on initial lines.
      bool degenerated;
      long lloydIter = lloydLines( samples, &l, &degenerated );

      //If lloyd degenerated restore original initial guesses
      if( degenerated )
	  l = lCopy;

      //Debug Output
      //cerr << "Degnerated: " << degenerated << endl;

      // Find the two most parallel lines
      vector<Vector> angles;
      Vector zeroVec(2, 0.0);
      for( int i = 0; i < clusters->size(); i++ )
      {
	  currPoint.copy(zeroVec);
	  currPoint.set(0, atan2( ((l.at(i))->p2)[1] - ((l.at(i))->p1)[1],  
				  fabs( ((l.at(i))->p2)[0] - ((l.at(i))->p1)[0] ) ) );
	  currPoint.set(1, (double)i );
	  angles.push_back( currPoint );
      }

      std::sort( angles.begin(), angles.end(), lloydQuadFinder::compareIndexPairs );

      //Order lines
      vector<Line*> lOrdered;
      for( int i = 0; i < clusters->size(); i++)
	  lOrdered.push_back( l.at( (unsigned int)(angles.at(i)[1]) ) );

      //Now search for most parallel lines
      int lowestDiff = 0;
      for( int i = 0; i < clusters->size() - 1; i++ )
      {
	  if( fabs( angles.at(i)[0] - angles.at( i + 1 )[0] ) < 
	      fabs( angles.at(lowestDiff)[0] - angles.at( lowestDiff + 1)[0] ) )
	      lowestDiff = i;
      }

      //Compute the four intersection Points and save them
      for( int j = 0; j < 2; j++ )
      {
	  for( int i = 0; i < clusters->size(); i++ )
	  {
	      if( i != lowestDiff && i != lowestDiff + 1 )
	      {
		  currPoint.copy(zeroVec);
		  interSect( *(lOrdered.at(i)), *(lOrdered.at(lowestDiff + j)), currPoint );
		  (corners->data).push_back( currPoint);
	      }
	  }
      }

      return;
  }


} /* namespace */

#endif /* CONNCOMPPROPS_LLOYDQUADFINDER_C */

