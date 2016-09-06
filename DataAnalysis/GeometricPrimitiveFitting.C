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

#ifndef GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_C
#define GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_C

#include "GeometricPrimitiveFitting.hh"

#include "MDA/LinearAlgebra/Matrix.hh"
#include "MDA/LinearAlgebra/JacobiRotation.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  /** fit a point to a SampleVector of data-points */
  void fitToDataPoints( Point &p, const SampleVector &dataPoints )
  {
      //Exit if dataPoints empty
      if( dataPoints.data.size() == 0 )
        return;

      Vector result(p.pos.getSize(), 0.0);
      for( unsigned long i = 0; i < dataPoints.data.size(); i++ )
      {
          //Exit if wrong size
          if( dataPoints.data.at(i).getSize() != p.pos.getSize() )
          {
              cerr << "Vector dimensions must match! " << std::endl;
              return;
          }

          result += dataPoints.data.at(i);
      }
      
      result /= (double)dataPoints.data.size();
      p.pos.assign( result );

      return;
  }

  /** fit a line to a SampleVector of data-points */
  void fitToDataPoints( Line &l, const SampleVector &dataPoints )
  {
      // The least squares line goes through the mean of all points
      // and the direction is the first eigenvector of the co-variance 
      // matrix of the points

      //Exit if datapoints empty
      if( dataPoints.data.size() == 0 )
        return;

      // First, compute arithmetic mean
      unsigned long dimension = l.p1.getSize();
      Vector mean(dimension, 0.0);
      for( unsigned long i = 0; i < dataPoints.data.size(); i++ )
      {
          //Exit if wrong size
          if( dataPoints.data.at(i).getSize() != l.p1.getSize() )
          {
              cerr << "Vector dimensions must match! " << std::endl;
              return;
          }

          mean += dataPoints.data.at(i);
      }
      
      mean /= (double)dataPoints.data.size();

      // Compute covariance-matrix
      Matrix cov( dimension, dimension, true );
      Vector covEntry(dimension, 0.0);
      for( unsigned long i = 0 ; i< dataPoints.data.size() ; i++ )
      {
          covEntry.copy( dataPoints.data.at(i) );
          covEntry -= mean;
          cov.addOuterProduct( covEntry );
      }

      // Solve Eigenvalue problem
      JacobiRotation solver;
      Vector eigValues( dimension, 0.0 );
      Matrix eigVectors( dimension, dimension, true);
      solver.solve( cov, eigValues, eigVectors );
      
      // Find the eigenVector with the larges eigenValue
      unsigned long largestEigenVector = 0;
      for( unsigned long i = 0; i < dimension; i++ )
      {
          if( eigValues[i] > eigValues[largestEigenVector] )
              largestEigenVector = i;
      }

      // Direction of the line is then
      Vector direction(dimension, 0.0);
      direction.assign( eigVectors.getRowVector(largestEigenVector) );

      // Store in 2-point representation
      l.p1.assign( mean );
      l.p2.assign( mean );
      l.p2 += direction;
      
      return;
  }


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_C */

