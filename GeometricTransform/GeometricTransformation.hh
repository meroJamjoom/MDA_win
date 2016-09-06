// ==========================================================================
// $Id: GeometricTransformation.hh 999 2014-05-28 15:07:31Z heidrich $
// Bas class for transformations of arrays
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

#ifndef GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_H
#define GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_H

/*! \file  Transformation.hh
    \brief Base class for transformations of arrays
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/Array.hh"
#include "MDA/Resampling/Resampling.hh"

namespace MDA {

  /** \class Transformation Transformation.hh
      Bas class for transformations of arrays */
  template<class T>
  class GeometricTransformation {

  public:

    /** default constructor */
    inline GeometricTransformation( Resampler<T> *_resampler= NULL )
    {
      if( _resampler== NULL )
      {
	resampler= new CubicResampler<T>;
	resampler->setBoundaryMethod( Clamp );
      }
      else
	resampler= _resampler;
    }
    
    /** destructor */
    virtual ~GeometricTransformation()
    {
      if( resampler!= NULL )
	delete resampler;
    }
    
    /** get a pointer to the current resampler */
    inline Resampler<T> *getResampler() const
    {
      return resampler;
    }
    
    /** set a new resampler */
    inline void setResampler( Resampler<T> *_resampler )
    {
      delete resampler;
      resampler= _resampler;
    }
    
    /** apply the transform to a single channel in one array , and
	write the result to another channel in a (possibly different)
	array. */
    virtual bool apply( Array<T> &srcArray, int srcChannel,
			Array<T> &dstArray, int dstChannel )= 0;
    
    /** apply the transform to all channels in the source array and
	write the result to new channels in the destination
	array. Optionally delete the source channels as soon as they
	are processed. */
    virtual bool apply( Array<T> &srcArray, Array <T> &dstArray,
			bool deleteSource= false );
    
    /** apply the transform to a list of channels, and write to a list
	of output channels (possibly in a different array) */
    virtual bool apply( Array<T> &srcArray, ChannelList srcChannels,
			Array<T> &dstArray, ChannelList dstChannels );
    
  protected:
    
    /** the resampler used for all transformations */
    Resampler<T> *resampler;
    
  };


} /* namespace */


#endif /* GEOMETRICTRANSFORM_GEOMETRICTRANSFORMATION_H */

