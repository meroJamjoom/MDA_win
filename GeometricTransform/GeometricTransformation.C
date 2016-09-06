// ==========================================================================
// $Id: GeometricTransformation.C 712 2010-04-23 09:16:34Z martinle $
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

#ifndef GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_C
#define GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_C

#include "MDA/Base/Errors.hh"

#include "GeometricTransformation.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  
/** apply the transform to all channels in the source array and write
    the result to new channels in the destination array. Optionally
    delete the source channels as soon as they are processed. */
template<class T>
bool
GeometricTransformation<T>::apply( Array<T> &srcArray, Array <T> &dstArray,
				   bool deleteSource )
{
  bool success= true;
  int numChannels= srcArray.getNumChannels();
  
  for( int i= 0 ; i< numChannels ; i++ )
  {
    if( deleteSource )
    {
      // processing&deleting from front, so next channel always has index 0
      success&= apply( srcArray, 0, dstArray, dstArray.addChannel() );
      srcArray.deleteChannel( 0 );
    }
    else
      success&= apply( srcArray, i, dstArray, dstArray.addChannel() );
  }
  return success;
}

 
/** apply the transform to a list of channels, and write to a list
    of output channels (possibly in different array) */
template<class T>
bool
GeometricTransformation<T>::apply( Array<T> &srcArray, ChannelList srcChannels,
				   Array<T> &dstArray, ChannelList dstChannels)
{
  int numChannels= srcChannels.vec.size();
  if( !warnCond( numChannels== dstChannels.vec.size(),
		 "  #destination channels does not match #source channels" ) )
    numChannels= numChannels< dstChannels.vec.size() ?
      numChannels : dstChannels.vec.size();
  bool success= true;
  for( int i= 0 ; i< numChannels ; i++ )
    success&=
      apply( srcArray, srcChannels.vec[i], dstArray, dstChannels.vec[i] );
  
  return success;
}


//
// template instantiations
//

template class GeometricTransformation<float>;
template class GeometricTransformation<double>;
template class GeometricTransformation<unsigned char>;

} /* namespace */

#endif /* GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_C */

