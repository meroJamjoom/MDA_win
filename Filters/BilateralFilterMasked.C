// ==========================================================================
// $Id: BilateralFilter.C 415 2009-10-15 10:51:23Z heidrich $
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

#ifndef FILTERS_BILATERALFILTERMASKED_C
#define FILTERS_BILATERALFILTERMASKED_C

#include <MDA/Array/Neighborhood.hh>
#include <MDA/Array/PaddedChannelTraverser.hh>
#include "BilateralFilterMasked.hh"

#include <math.h>
#include <fstream>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** apply the filter to a number of dimensions and channels */
template<class T>
bool
BilateralFilterMasked<T>::apply( Array<T> &a, BoundaryMethod boundary,
			   ChannelList &channels, AxisList &axes )
{
  unsigned long i, j;
  
  CoordinateVector dim= a.getDimension();
  unsigned dimension= dim.vec.size();
  unsigned numChannels= channels.vec.size();
  
  //
  // unlike most filters, the bilateral filter needs to process all
  // channels at once, hence they all need to be padded at the same
  // time
  //
  
  // compute number of pixels in both padded and unpadded arrays
  unsigned long numPixels= dim.vec[0];
  unsigned long numScanlines= scanlines( dim );
  CoordinateVector paddedDim= pad( dim, radius );
  unsigned long numPaddedPixels= size( paddedDim );
  
  // create padded channels and cache pointers to unpadded channel data
  T **oChannels= new T*[numChannels];
  T **paddedChannels= new T*[numChannels];
  for( i= 0 ; i< numChannels ; i++ )
  {
    oChannels[i]= &((*a[channels.vec[i]])[0]);
    paddedChannels[i]= new T[numPaddedPixels];
    fetchChannel( paddedChannels[i], oChannels[i], dim, radius,
		  boundary, a[channels.vec[i]]->getBackground() );
  }
  
  // create a neighborhood of the specified radius along specified axes
  Neighborhood n( dim, radius, axes, radius );
  unsigned neighborhoodSize= n.getNumPixels();
  
  // precompute spatial weights for each pixel in the neighborhood as
  // a Gaussian of their distance
  vector<double> spatialWeights( neighborhoodSize );
  for( i=0,n.begin() ; !n.isAtEnd() ; ++n,i++ ){
    spatialWeights[i]= exp( -(*n).dist*(*n).dist / (2.0*sigma*sigma) );
    //cerr<<(*n).dist<<" "<<spatialWeights[i]<<endl;
  }
  //
  // traverse all pixels and apply filter
  //
  
  double weight, weightSum;
  double edgeStopMult= -1.0 / (2.0 * edgeStopSigma * edgeStopSigma);
  double *accum= new double[numChannels];
  unsigned long uOff, pOff;
  
  // process all pixels in sequence (all channels at once)
  PaddedChannelTraverser traverser( dim, radius );

  for( traverser.begin() ; !traverser.isAtEnd() ; ++traverser )
  {
    // current positions within padded and original channels
    uOff= traverser.uOffset();
    pOff= traverser.pOffset();
    
    //excluding zero mask pixels
    if (oChannels[numChannels-1][uOff]==0) continue;

    // for each pixel, first clear the accumulator for each channel
    for( i= 0 ; i< numChannels ; i++ )
      accum[i]= 0.0;
    weightSum= 0.0;
    
    // compute filtered value for one pixel by iterating over its neighborhood
    for( j= 0, n.begin() ; !n.isAtEnd() ; ++n, j++ )
    {

      //excluding zero pixels
      if (paddedChannels[numChannels-1][pOff+(*n).pOff]==0)  continue;
      
      // get photometric distance squared
      weight= 0.0;
      for( i= 0 ; i< numChannels-1 ; i++ )
      {
	double tmp=
	  oChannels[i][uOff] - paddedChannels[i][pOff+(*n).pOff];
	weight+= tmp*tmp;
      }
      // compute full weight formula (Gaussian of photo distance)
      weight= exp( weight * edgeStopMult )*spatialWeights[j];
      weightSum+= weight;

      // accumulate weighted values
      for( i= 0 ; i< numChannels-1 ; i++ )
	accum[i]+= weight * paddedChannels[i][pOff+(*n).pOff];
    }
    
    // finally, normalize pixel result
    for( i= 0 ; i< numChannels-1 ; i++ )
    	oChannels[i][uOff]= accum[i] / weightSum;
  }
  //
  // clean up
  //
  
  delete [] accum;
  for( i= 0 ; i< numChannels ; i++ )
    delete [] paddedChannels[i];
  delete [] paddedChannels;
  delete [] oChannels;

  return true;
}
  
  
// template instantiation code
template class BilateralFilterMasked<float>;
template class BilateralFilterMasked<double>;

} /* namespace */

#endif /* FILTERS_BILATERALFILTER_C */

