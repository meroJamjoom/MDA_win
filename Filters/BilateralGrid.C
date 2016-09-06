// ==========================================================================
// $Id: BilateralGrid.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef FILTERS_BILATERALGRID_C
#define FILTERS_BILATERALGRID_C

#include <stdio.h>

#include "MDA/Resampling/PointSampler.hh"
#include "Linear1DFilter.hh"
#include "BilateralGrid.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** "construct" a new grid */
template<class T>
bool
BilateralGrid<T>::construct( Array<T> &a, ChannelList &dataChannels,
			     ChannelList &edgeStopChannels )
{
  unsigned long i, j, k;
  
  // source dimensions
  CoordinateVector srcDim= a.getDimension();
  unsigned srcDimension= srcDim.vec.size();
  unsigned long srcChannelsize= 1;
  for( i= 0 ; i< srcDimension ; i++ )
    srcChannelsize*= srcDim.vec[i];
  
  // grid dimensions
  // first, the spatial dimensions
  CoordinateVector gridDim;
  minVal.clear();
  maxVal.clear();
  for( i= 0 ; i< srcDim.vec.size() ; i++ )
    gridDim.vec.push_back( (int)(srcDim.vec[i]/spatialScale)+1 );
  // then, then intensity dimensions
  for( i= 0 ; i< edgeStopChannels.vec.size() ; i++ )
  {
    // compute min/max for the channel first
    minVal.push_back( (*a[edgeStopChannels.vec[i]])[0] );
    maxVal.push_back( (*a[edgeStopChannels.vec[i]])[0] );
    for( j= 1 ; j< srcChannelsize ; j++ )
    {
      T val= (*a[edgeStopChannels.vec[i]])[j];
      if( minVal[i]> val )
	minVal[i]= val;
      if( maxVal[i]< val )
	maxVal[i]= val;
    }
    gridDim.vec.push_back( (int)((maxVal[i]-minVal[i]) / intensityScale) + 1 );
  }
  unsigned gridDimension= gridDim.vec.size();
  unsigned long gridChannelsize= 1;
  for( i= 0 ; i< gridDimension ; i++ )
    gridChannelsize*= gridDim.vec[i];
  unsigned numGridChannels=  dataChannels.vec.size()+1;
  
  // warning for large grids
  warnMem( gridChannelsize*(dataChannels.vec.size()+1)*sizeof(T),
	   MEGA_BYTES( 100 ) );
  
  // allocate grid
  if( grid!= NULL )
    clear();
  grid= new Array<T>( gridDim );
  // allocate and initialize n+1 channels, the last one for the weights
  vector<T *> gridChannelPtrs;
  for( i= 0 ; i< numGridChannels ; i++ )
  {
    grid->addChannel();
    gridChannelPtrs.push_back( &((*(*grid)[i])[0]) );
    bzero( gridChannelPtrs[i], gridChannelsize );
  }
  
  // finally, actually fill the grid with image data
  CoordinateVector srcPos( srcDimension );
  Vector gridPos( gridDimension );
  vector<T *> srcChannelPtrs;
  for( i= 0 ; i< dataChannels.vec.size() ; i++ )
    srcChannelPtrs.push_back( &((*a[dataChannels.vec[i]])[0]) );
  NearestNeighborPointSampler<T> splatter( gridDim );
  for( i= 0 ; i< srcChannelsize ; i++ )
  {
    // set spatial dimensions of grid pos
    for( j= 0 ; j< srcDimension ; j++ )
      gridPos[j]= srcPos.vec[j] / spatialScale;
    // set intesity dimensions of grid pos
    for( k= 0 ; j< gridDimension ; j++, k++ )
      gridPos[j]=
	(srcChannelPtrs[k][i] - minVal[k]) / intensityScale;
    
    // set sample location
    splatter.setSampleLocation( gridPos, Clamp );
    
    // splat each channel and the weights
    for( j= numGridChannels-1 ; j> 0 ; )
    {
      j--;
      splatter.splatSample( srcChannelPtrs[j][i], gridChannelPtrs[j] );
    }
    splatter.splatSample( 1.0, gridChannelPtrs[numGridChannels-1] );
    
    // update position in source array
    for( j= 0 ; j< srcDimension ; j++ )
      if( ++srcPos.vec[j]>= srcDim.vec[j] )
	srcPos.vec[j]= 0;
      else
	break;
  }
  
  return true;
}


/** "process" an existing grid */
template<class T>
bool
BilateralGrid<T>::process( BoundaryMethod boundary, AxisList &axes,
			   double sigma, double edgeStopSigma )
{
  if( grid== NULL )
    return false;
  
  unsigned i, j;

  // filter spatial dimensions
  FastGaussian1D<T> spatialFilter( sigma );
  for( i= 0 ; i< axes.vec.size() ; i++ )
    for( j= 0 ; j< grid->getNumChannels() ; j++ )
      spatialFilter.apply( *grid, boundary, i, j, j );
  
  // filter intensity dimensions
  FastGaussian1D<T> intensityFilter( edgeStopSigma );
  for( ; i< grid->getDimension().vec.size() ; i++ )
    for( j= 0 ; j< grid->getNumChannels() ; j++ )
      intensityFilter.apply( *grid, boundary, i, j, j );
  
  return true;
}


/** "slice" an existing grid */
template<class T>
bool
BilateralGrid<T>::slice( Array<T> &a, ChannelList &dataChannels,
			 ChannelList &edgeStopChannels )
{
  unsigned long i, j, k;
  
  if( grid== NULL )
    return false;

  // destination info
  CoordinateVector dstDim= a.getDimension();
  unsigned dstDimension= dstDim.vec.size();
  unsigned dstChannelsize= 1;
  for( i= 0 ; i< dstDimension ; i++ )
    dstChannelsize*= dstDim.vec[i];
  
  // grid info
  CoordinateVector gridDim= grid->getDimension();
  unsigned gridDimension= gridDim.vec.size();
  unsigned numGridChannels= grid->getNumChannels();
  vector<T *> gridChannelPtrs;
  for( i= 0 ; i< numGridChannels ; i++ )
    gridChannelPtrs.push_back( &((*(*grid)[i])[0]) );
  
  // protect against some basic internal errors
  errorCond( dstDimension+dataChannels.vec.size() == gridDimension,
	     "  inconsistent grid dimensions" );
  errorCond( dataChannels.vec.size() == edgeStopChannels.vec.size(),
	     "  inconsistent #channels" );
  
  // finally, actually fill the image data with grid info
  CoordinateVector dstPos( dstDimension );
  Vector gridPos( gridDimension );
  vector<T *> dstChannelPtrs;
  vector<T *> edgeStopChannelPtrs;
  for( i= 0 ; i< dataChannels.vec.size() ; i++ )
  {
    dstChannelPtrs.push_back( &((*a[dataChannels.vec[i]])[0]) );
    edgeStopChannelPtrs.push_back( &((*a[edgeStopChannels.vec[i]])[0]) );
  }
  LinearPointSampler<T> reconstructor( gridDim );
  for( i= 0 ; i< dstChannelsize ; i++ )
  {
    // set spatial dimensions of grid pos
    for( j= 0 ; j< dstDimension ; j++ )
      gridPos[j]= dstPos.vec[j] / spatialScale;
    // set intesity dimensions of grid pos
    for( k= 0 ; j< gridDimension ; j++, k++ )
      gridPos[j]=
	(edgeStopChannelPtrs[k][i] - minVal[k]) / intensityScale;
    
    // set sample location
    reconstructor.setSampleLocation( gridPos, Clamp );
    
    // splat each channel and the weights
    double weight= reconstructor.getSample(gridChannelPtrs[numGridChannels-1]);
    for( j= numGridChannels-1 ; j> 0 ; )
    {
      j--;
      dstChannelPtrs[j][i]=
	(weight==0.0 ? 0 : reconstructor.getSample( gridChannelPtrs[j] ) / weight);
    }
    
    // update position in source array
    for( j= 0 ; j< dstDimension ; j++ )
      if( ++dstPos.vec[j]>= dstDim.vec[j] )
	dstPos.vec[j]= 0;
      else
	break;
  }
  return true;
}
 

/*********************************************************************/

/** "construct" a new grid */
template<class T>
bool
BilateralGridWeighted<T>::construct( Array<T> &a, ChannelList &dataChannels,
			     ChannelList &edgeStopChannels )
{
  unsigned long i, j, k;
  
  // source dimensions
  CoordinateVector srcDim= a.getDimension();
  unsigned srcDimension= srcDim.vec.size();
  unsigned long srcChannelsize= 1;
  for( i= 0 ; i< srcDimension ; i++ )
    srcChannelsize*= srcDim.vec[i];
  
  // grid dimensions
  // first, the spatial dimensions
  CoordinateVector gridDim;
  BilateralGrid<T>::minVal.clear();
  BilateralGrid<T>::maxVal.clear();
  for( i= 0 ; i< srcDim.vec.size() ; i++ )
    gridDim.vec.push_back( (int)(srcDim.vec[i]/BilateralGrid<T>::spatialScale)+1 );
  // then, then intensity dimensions
  //for( i= 0 ; i< edgeStopChannels.vec.size() ; i++ )
  for( i= 0 ; i< edgeStopChannels.vec.size()-1 ; i++ )
  {
    // compute min/max for the channel first
    BilateralGrid<T>::minVal.push_back( (*a[edgeStopChannels.vec[i]])[0] );
    BilateralGrid<T>::maxVal.push_back( (*a[edgeStopChannels.vec[i]])[0] );
    for( j= 1 ; j< srcChannelsize ; j++ )
    {
      T val= (*a[edgeStopChannels.vec[i]])[j];
      if( BilateralGrid<T>::minVal[i]> val )
	BilateralGrid<T>::minVal[i]= val;
      if( BilateralGrid<T>::maxVal[i]< val )
	BilateralGrid<T>::maxVal[i]= val;
    }
    gridDim.vec.push_back( (int)((BilateralGrid<T>::maxVal[i]-BilateralGrid<T>::minVal[i]) / BilateralGrid<T>::intensityScale) + 1 );
  }
  unsigned gridDimension= gridDim.vec.size();
  unsigned long gridChannelsize= 1;
  for( i= 0 ; i< gridDimension ; i++ )
    gridChannelsize*= gridDim.vec[i];
  //unsigned numGridChannels=  dataChannels.vec.size()+1;
  unsigned numGridChannels=  dataChannels.vec.size(); // last one is weights
  
  // warning for large grids
  warnMem( gridChannelsize*(dataChannels.vec.size()+1)*sizeof(T),
	   MEGA_BYTES( 100 ) );
  
  // allocate grid
  if( BilateralGrid<T>::grid!= NULL )
    BilateralGrid<T>::clear();
  BilateralGrid<T>::grid= new Array<T>( gridDim );
  // allocate and initialize n+1 channels, the last one for the weights
  vector<T *> gridChannelPtrs;
  for( i= 0 ; i< numGridChannels ; i++ )
  {
    BilateralGrid<T>::grid->addChannel();
    gridChannelPtrs.push_back( &((*(*BilateralGrid<T>::grid)[i])[0]) );
    bzero( gridChannelPtrs[i], gridChannelsize );
  }
  
  // finally, actually fill the grid with image data
  CoordinateVector srcPos( srcDimension );
  Vector gridPos( gridDimension );
  vector<T *> srcChannelPtrs;

  //assume last channel is weights, so don't use it for anything else
  unsigned long weightChannel =  dataChannels.vec.size()-1;
  //for( i= 0 ; i< dataChannels.vec.size() ; i++ )
  for( i= 0 ; i< dataChannels.vec.size()-1 ; i++ )
    srcChannelPtrs.push_back( &((*a[dataChannels.vec[i]])[0]) );
  NearestNeighborPointSampler<T> splatter( gridDim );
  for( i= 0 ; i< srcChannelsize ; i++ )
  {
    // set spatial dimensions of grid pos
    for( j= 0 ; j< srcDimension ; j++ )
      gridPos[j]= srcPos.vec[j] / BilateralGrid<T>::spatialScale;
    // set intesity dimensions of grid pos
    for( k= 0 ; j< gridDimension ; j++, k++ )
      gridPos[j]=
	(srcChannelPtrs[k][i] - BilateralGrid<T>::minVal[k]) / BilateralGrid<T>::intensityScale;
    
    // set sample location
    splatter.setSampleLocation( gridPos, Clamp );
    
    // splat each channel and the weights
    for( j= numGridChannels-1 ; j> 0 ; )
    {
      j--;
      //splatter.splatSample( srcChannelPtrs[j][i], gridChannelPtrs[j] );
      splatter.splatSample( srcChannelPtrs[j][i]*(*a[dataChannels.vec[weightChannel]])[i], gridChannelPtrs[j] );
    }
    //splatter.splatSample( 1.0, gridChannelPtrs[numGridChannels-1] );
    splatter.splatSample( (*a[dataChannels.vec[weightChannel]])[i], gridChannelPtrs[numGridChannels-1] );
    
    // update position in source array
    for( j= 0 ; j< srcDimension ; j++ )
      if( ++srcPos.vec[j]>= srcDim.vec[j] )
	srcPos.vec[j]= 0;
      else
	break;
  }

  // now remove the last channel as it was only weights
  dataChannels.vec.pop_back();
  
  return true;
}


// explicit template instation code

template class BilateralGrid<float>;
template class BilateralGrid<double>;
template class BilateralGridWeighted<float>;
template class BilateralGridWeighted<double>;


} /* namespace */

#endif /* FILTERS_BILATERALGRID_C */

