// ==========================================================================
// $Id: SampleVector.C 224 2008-10-25 04:52:07Z heidrich $
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

#ifndef DATAANALYSIS_SAMPLEVECTOR_C
#define DATAANALYSIS_SAMPLEVECTOR_C

#include "SampleVector.hh"

namespace MDA {

using namespace std;


/** create a SampleVector from Array Channels */
SampleVector::SampleVector( Array<float> &a, ChannelList &ch )
{
  unsigned numChannels= ch.vec.size();
  vector<float*> channelPtrs( numChannels );
  unsigned long i, j;
  for( i= 0 ; i< numChannels ; i++ )
    channelPtrs[i]= &((*a[ch.vec[i]])[0]);
  unsigned long numElem= 1;
  CoordinateVector dim= a.getDimension();
  for( i= 0 ; i< dim.vec.size() ; i++ )
    numElem*= dim.vec[i];
  
  data.reserve( numElem );
  data.clear();
  for( i= 0 ; i< numElem ; i++ )
  {
    data.push_back( Vector( numChannels ) );
    for( j= 0 ; j< numChannels ; j++ )
      data[i][j]= channelPtrs[j][i];
  }
}

/** create a SampleVector from Array Channels */
SampleVector::SampleVector( Array<double> &a, ChannelList &ch )
{
  unsigned numChannels= ch.vec.size();
  unsigned long numElem= 1;
  vector<double *> channelPtrs( numChannels );
  unsigned long i, j;
  for( i= 0 ; i< numChannels ; i++ )
  {
    channelPtrs[i]= &((*a[ch.vec[i]])[0]);
    numElem*= ch.vec[i];
  }
  
  data.reserve( numElem );
  data.clear();
  for( i= 0 ; i< numElem ; i++ )
  {
    data.push_back( Vector( numChannels ) );
    for( j= 0 ; j< numChannels ; j++ )
      data[i][j]= channelPtrs[j][i];
  }
}


} /* namespace */

#endif /* DATAANALYSIS_SAMPLEVECTOR_C */

