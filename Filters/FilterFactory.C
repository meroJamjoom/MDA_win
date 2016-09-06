// ==========================================================================
// $Id: FilterFactory.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef FILTERS_FILTERFACTORY_C
#define FILTERS_FILTERFACTORY_C

#include "MDA/Base/Errors.hh"

#include "FilterFactory.hh"
#include "Filter1D.hh"
#include "Linear1DFilter.hh"
#include "SeparableFilter.hh"
#include "UnsharpMasking.hh"
#include "SobelFilter.hh"
#include "SimpleEdgeFilter.hh"
#include "FloodFill.hh"
#include "MedianFilter.hh"
#include "MorphologicalOps.hh"
#include "BilateralFilter.hh"
#include "BilateralFilterMasked.hh"
#include "BilateralGrid.hh"
#include "Thinning2D.hh"
#include "ImprovedHarrisCorner.hh"
#include "ExtremaDetector.hh"
#include "ConnectedComponent.hh"
#include "Prune2D.hh"
#include "EulerNumber2D.hh"
#include "DistanceTransform.hh"
#include "Thinning3D.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;



/** create a new filter of filterType using previously set parameters */
template<class T>
Filter<T> *
FilterFactory<T>::create( FilterType filterType )
{
  int i;
  Filter1D<T> *f1d= NULL;
  int radius= parameters->radius;
  double sigma= parameters->sigma;
  double edgeStopSigma= parameters->edgeStopSigma;
  
  switch( filterType )
  {
  case BilateralFiltering: // brute force bilateral filtering
    return new BilateralFilter<T>( sigma, edgeStopSigma,
				   radius> 0 ? radius : (int)(2.0*sigma+.5) );
  case BilateralFilteringMasked: // brute force bilateral filtering (masked)
    return new BilateralFilterMasked<T>( sigma, edgeStopSigma,
				   radius> 0 ? radius : (int)(2.0*sigma+.5) );
  case BilateralGridFiltering: // fast bilateral filter
    return new BilateralGrid<T>( sigma, edgeStopSigma );
  case BilateralGridWeightedFiltering: // fast bilateral filter with weights
    return new BilateralGridWeighted<T>( sigma, edgeStopSigma );
  case BoxFiltering: // box
    f1d= new BoxFilter1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case ConnectedComponentFiltering: // connected components
    return new ConnectedComponent<T>;
  case CornerFiltering: // Harris corner detector
    return new ImprovedHarrisCorner<T>( sigma );
  case DilateFiltering: // Morphological dilate
    f1d= new Dilate1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case DistanceFiltering: // Distance transform
    return new DistanceTransform<T>;
  case EdgeFiltering: // Simple edge filter
    return new SimpleEdgeFilter<T>( radius> 0 ? radius : 1 );
  case ErodeFiltering: // Morphological erode
    f1d= new Erode1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case EulerNumberFiltering: // filter for Euler number
    return new EulerNumber2D<T>;
  case ExtremaFiltering: // find and label local extrema
    return new ExtremaDetector<T>( radius> 0 ? radius : 1 );
  case FastGaussFiltering: // fastgauss
    f1d= new FastGaussian1D<T>( sigma );
    return new SeparableFilter<T>( f1d );
  case FirstDerivativeFiltering: // Derivative
    f1d= new FirstDerivative1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case FloodfillFiltering: // Flood fill
    return new FloodFill<T>;
  case GaussFiltering: // gauss
    f1d= new Gaussian1D<T>( sigma,
			    radius> 0 ? radius : (int)(2.0*sigma+.5) );
    return new SeparableFilter<T>( f1d );
  case HatFiltering: // hat
    f1d= new HatFilter1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case HitOrMissFiltering: // morphological hit-or-miss operator
    {
      errorCond( parameters->values.size()== 9,
		 "Hit-or-miss operator requres exactly 9 filter values\n" );
      bool structuringElement[9];
      bool mask[9];
      for( i= 0 ; i< 9 ; i++ )
      {
	double h= parameters->values[i]->eval( NULL, NULL );
	structuringElement[i]= (h> 0.0);
	mask[i]= (h>= 0.0);
      }
      return new HitOrMiss2D<T>( structuringElement, mask );
    }
  case LoGFiltering: // Laplacian of Gaussian
    f1d= new LaplacianOfGaussian1D<T>( sigma,
					  radius> 0 ? radius :
					  (int)(4.0*sigma+.5) );
    return new SeparableFilter<T>( f1d );
  case MedianFiltering: // Median
    return new MedianFilter<T>( radius> 1 ? radius : 1 );
  case MedianMaskFiltering: // Median
    return new MedianFilterMasked<T>( radius> 1 ? radius : 1 );
  case PruneFiltering: // pruning
    return new Prune2D<T>;
  case SecondDerivativeFiltering: // Second derivative
    f1d= new SecondDerivative1D<T>( radius> 0 ? radius : 1 );
    return new SeparableFilter<T>( f1d );
  case SeparableFiltering: // Custom separable
    {
      double *valArray;
      double distance;
      double norm= 0.0;
      
      // if the expression sequence has only one entry, we use the
      // radius to determine the filter width, and apply the same
      // formula to all entries
      // if more than one expression is provided, then we determine
      // the new radius from the number of expressions, and evaluate
      // a different expression for each filter entry
      if( parameters->values.size()== 1 )
      {
	valArray= new double[2*radius+1];
	for( i= 0 ; i<= 2*radius ; i++ )
	{
	    distance= i-radius;
	    norm+= valArray[i]= parameters->values[0]->eval( &distance, NULL );
	}
      }
      else
      {
	radius= parameters->values.size()%2 ?
	  (parameters->values.size()-1)/2 : parameters->values.size()/2;
	valArray= new double[2*radius+1];
	for( i= 0 ; i<= 2*radius ; i++ )
	  if( i>= parameters->values.size() )
	    valArray[i]= 0.0;
	  else
	  {
	    distance= i-radius;
	    norm+= valArray[i]= parameters->values[i]->eval( &distance, NULL );
	  }
      }
      
      // normalize the filter if the used desires it
      if( parameters->normalize )
	for( i= 0 ; i<= 2*radius ; i++ )
	  valArray[i]/= norm;
      
      // create the filter
      f1d= new Linear1DFilter<T>( radius );
      ((Linear1DFilter<T> *)f1d)->setFilter( radius, valArray );
      delete [] valArray;
      return new SeparableFilter<T>( f1d );
    }
    break;
  case SobelFiltering: // Sobel
    return new SobelFilter<T>;
  case ThinningFiltering: // Thinning
    return new Thinning2D<T>();
  case UnsharpMaskFiltering: // unsharp masking
    return new UnsharpMasking<T>( sigma );
  case Thinning3DFiltering: // Thinning 3D
    return new Thinning3D<T>();
  default: // UndefinedFiltering
    warning( "  unsupported filter" );
  }

  // if we get here, the filter was unknown...
  return NULL;
}


template class FilterFactory<float>;  
template class FilterFactory<double>;  
  
} /* namespace */
  
#endif /* FILTERS_FILTERFACTORY_C */

