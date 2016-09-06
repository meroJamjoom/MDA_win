// ==========================================================================
// $Id: FloodFill.C 301 2009-03-25 04:16:15Z heidrich $
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

#ifndef FILTERS_FLOODFILL_C
#define FILTERS_FLOODFILL_C

#include "MDA/Base/Errors.hh"

#include "FloodFill.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


// begin template definitions

/** apply flood fill to a number of dimensions and channels */
template <class T>
bool
FloodFill<T>::apply( Array<T> &a, BoundaryMethod boundary,
		     ChannelList &channels, AxisList &axes )
{
  unsigned long i, j;
  unsigned long minX, maxX;
  unsigned long lineIndex;
  CoordinateVector pos;
  
  list<ScanSegment> segments;
  CoordinateVector dim= a.getDimension();
  unsigned dimension= dim.vec.size();

  // get mask channel (and remove it from regular channel list)
  unsigned maskChannel= channels.vec[channels.vec.size()-1];
  channels.vec.pop_back();
  typename Array<T>::Channel *mask= a[maskChannel];
  if( !warnCond( mask!= NULL, "  mask channel out of range" ) )
    return false;
  

  // check if axis 0 is part of the axis list
  int k;
  for( k= axes.vec.size()-1 ; k>= 0 ; k-- )
    if( axes.vec[k]== 0 )
      break;
  warnCond( k>= 0,
	    "  FloodFill requires axis 0 to be part of axis list\n"
	    "  (continuing as if it was...)" );
  
  // compute per-dimension increments
  CoordinateVector incr;
  incr.vec.push_back( 1 );
  for( i= 1 ; i< dimension ; i++ )
    incr.vec.push_back( incr.vec[i-1]*dim.vec[i-1] );
  
  // process all channels independently (they could have different seeds...)
  for( i= 0 ; i< channels.vec.size() ; i++ )
  {
    // get channel
    typename Array<T>::Channel *channel= a[channels.vec[i]];
    if( channel== NULL )
      continue;

    // position of first scanline
    pos.vec.clear();
    for( j= 0 ; j< dimension ; j++ )
      pos.vec.push_back( 0 );
    lineIndex= 0;
    
    // initialize list with all seed segments:
    // traverse all scanlines
    while( pos.vec[dimension-1]< dim.vec[dimension-1] )
    {
      // add seed segments for this scanline
      for( j= 0 ; j< dim.vec[0] ; j++ )
	if( (*channel)[lineIndex+j]== fillVal )
	{
	  minX= j;
	  for( ; j< dim.vec[0]-1 && (*mask)[lineIndex+j+1]!= boundaryVal ; j++ )
	    ;
	  segments.push_back( ScanSegment( pos, lineIndex, minX, j ) );
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
    
    growSeed( &((*channel)[0]), &((*mask)[0]), segments, dim, axes, incr );
  }
  
  return true;
}

/** flood fill a single channel from a given (set of) seed segment(s) */
template <class T>
void
FloodFill<T>::growSeed( T *channel, T *mask,
			list<ScanSegment> &segments, CoordinateVector &dim,
			AxisList &axes, CoordinateVector &incr )
{
  unsigned long i, j;
  unsigned long minX, maxX;
  unsigned long lineIndex;
  CoordinateVector pos;
  
  // now process all segments until none is left
  while( !segments.empty() )
  {
    // get segment info
    minX= segments.front().minX;
    maxX= segments.front().maxX;
    lineIndex= segments.front().lineIndex;
    pos= segments.front().pos;
    segments.pop_front();
    
    // extend segment in negative x direction
    if( mask[lineIndex+minX]!= boundaryVal && minX> 0)
      while( mask[lineIndex+minX-1]!= boundaryVal && minX> 0 )
	minX--;
    // extend segment in positive x direction
    if( mask[lineIndex+maxX]!= boundaryVal && maxX< dim.vec[0]-1)
      while( mask[lineIndex+maxX+1]!= boundaryVal && maxX< dim.vec[0]-1 )
	maxX++;
    
    // traverse the segment
    while( minX<= maxX )
    {
      // skip boundary pixels
      while( minX<= maxX && mask[lineIndex+minX]== boundaryVal )
	minX++;
      
      if( minX> maxX )
	break;
	
      // search for right boundary (end of sub-segment)
      bool pixelSet= false;
      for( j= minX ; j<= maxX ; j++ )
      {
	// exit loop if boundary found
	if( mask[lineIndex+j]== boundaryVal )
	  break;
	// otherwise, see if we have an unset interior pixel
	if( channel[lineIndex+j]!= fillVal )
	{
	  channel[lineIndex+j]= fillVal;
	  pixelSet= true;
	}
      }
      // sub-segment ranges from minX to j-1
      
      
      // if any action on this segment, add all spatially adjacent
      // segments to the back of the list
      if( pixelSet )
      {
	for( i= 0 ; i< axes.vec.size() ; i++ )
	{
	  int dir= axes.vec[i];
	  if( dir== 0 )
	    continue;
	  if( pos.vec[dir]> 0 )
	  {
	    segments.push_front( ScanSegment( pos, lineIndex-incr.vec[dir],
					      minX, j-1 ) );
	    segments.front().pos.vec[dir]--;
	  }
	  if( pos.vec[dir]< dim.vec[dir]-1 )
	  {
	    segments.push_front( ScanSegment( pos, lineIndex+incr.vec[dir],
					      minX, j-1 ) );
	    segments.front().pos.vec[dir]++;
	  }
	}
      }
      
      // search for next sub-segment
      minX= j;
    }
  }
}
    



/* we use explicit instantiation for this class */
template class FloodFill<float>;
template class FloodFill<double>;

  
} /* namespace */

#endif /* FILTERS_FLOODFILL_C */

