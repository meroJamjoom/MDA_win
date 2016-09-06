// ==========================================================================
// $Id: ConnectedComponent.C 999 2014-05-28 15:07:31Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef FILTERS_CONNECTEDCOMPONENT_C
#define FILTERS_CONNECTEDCOMPONENT_C

#include "ConnectedComponent.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** find connected components in a number of channels */
template <class T>
bool
ConnectedComponent<T>::apply( Array<T> &a, BoundaryMethod boundary,
			      ChannelList &channels, AxisList &axes )
{
  unsigned long i, j;
  unsigned long lineIndex;
  CoordinateVector pos;
  
  list<typename FloodFill<T>::ScanSegment> segment;
  CoordinateVector dim= a.getDimension();
  unsigned dimension= dim.vec.size();
  
  // compute per-dimension increments
  CoordinateVector incr;
  incr.vec.push_back( 1 );
  for( i= 1 ; i< dimension ; i++ )
    incr.vec.push_back( incr.vec[i-1]*dim.vec[i-1] );
  
  // process all channels independently
  for( i= 0 ; i< channels.vec.size() ; i++ )
  {
    // reset fill value to 0
    // (indicating no components found in current channel)
    numComp= 0;
    
    // get channel
    typename Array<T>::Channel *channel= a[channels.vec[i]];
    if( channel== NULL )
      continue;
    
    // position of first scanline
    pos.vec.resize( dimension, 0 );
    for( j= 0 ; j< dimension ; j++ )
      pos.vec[i]= 0;
    lineIndex= 0;
    
    // traverse all scanlines
    while( pos.vec[dimension-1]< dim.vec[dimension-1] )
    {
      // pointer to start of line
      T *linePtr= &((*channel)[lineIndex]);
      
      // find first pixel that is neither background, nor already
      // labelled
      for( j= 0 ; j< dim.vec[0] ; j++ )
	if( linePtr[j]> 0.0 && linePtr[j]< 1.0 )
	{
	  segment.push_back( typename FloodFill<T>::ScanSegment( pos,
								 lineIndex,
								 j, j ) );
	  FloodFill<T>::fillVal= (T)(++numComp);
	  this->growSeed( &((*channel)[0]), &((*channel)[0]),
			  segment, dim, axes, incr );
	}
      
      // update scanline position
      lineIndex+= incr.vec[1];
      for( j= 1 ; j< dimension ; j++ )
      {
	if( ++pos.vec[j]< dim.vec[j] )
	  break;
	else
	  if( j< dimension-1 )
	    pos.vec[j]= 0;
      }
    }
    
  }
  return true;
}



/* we use explicit instantiation for this class */
template class ConnectedComponent<float>;
template class ConnectedComponent<double>;


} /* namespace */

#endif /* FILTERS_CONNECTEDCOMPONENT_C */

