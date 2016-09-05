// ==========================================================================
// $Id: Boundary.hh 423 2009-10-16 18:18:15Z skoch $
// handling of boundaries
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

#ifndef ARRAY_BOUNDARY_H
#define ARRAY_BOUNDARY_H

/*! \file  Boundary.hh
    \brief handling of boundaries in arrays
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CoordinateVector.hh"
#include "MDA/Base/CommandlineParser.hh"


namespace MDA {

  /** BoundaryMethod: handling of boundaries in arrays */
  enum BoundaryMethod
  {
    Background,
    Clamp,
    Cyclic,
    Mirror,
    Renormalize
  };
  
  /** write boundary mode name to istream */
  ostream & operator<<( ostream &os, BoundaryMethod boundary );

  /** read boundary mode name from istream */
  istream & operator>>( istream &is, BoundaryMethod &boundary );
  
  
  /** \class BoundaryOption Boundary.hh
      parser for a boundary mode */
  class BoundaryOption: public CommandlineOption {
    
  public:
    
    /** constructor: optional paramters are help string and option texts */
    BoundaryOption( BoundaryMethod &b,
		    const char *msg= "\tboundary mode: one of "
		    "background, clamp, cyclic, mirror, renormalize\n",
		    const char *lTxt= "--boundary-mode", const char *sTxt= "-b")
      : CommandlineOption( msg, lTxt, sTxt ), boundary( b )
    {}

    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );

    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** current boundary value */
    BoundaryMethod &boundary;
  };
  
  
  /* copy an axis-aligned line from an array into a separate buffer.
     The destination buffer will also contain replicated data at the
     boundaries, according to the boundary mode */
  template<class T> inline void fetchLine( T *dstBuf, T *srcBuf,
                                          unsigned long srcStride, 
                                          unsigned long numElements,
                                          unsigned long radius,
                                          BoundaryMethod boundary, 
                                          T background ){
    fetchLine( dstBuf,srcBuf, srcStride, numElements, radius, radius,
                     boundary, background);
    return;
  };
  
  /* copy an axis-aligned line from an array into a separate buffer.
     more general than fetchLine, padding can be different at 
     front and back */
  
  template<class T> void fetchLine( T *dstBuf, T *srcBuf,
                                   unsigned long srcStride, 
                                   unsigned long numElements,
                                   unsigned long paddingleft, 
                                   unsigned long paddingright,
                                   BoundaryMethod boundary, T background );

  
  /* copy a channel into a separate buffer.
     The destination buffer will also contain replicated data at the
     boundaries, according to the boundary mode */
  template<class T> inline void fetchChannel( T *dstBuf, T *srcBuf,
                                             CoordinateVector srcRes,
                                             unsigned long radius,
                                             BoundaryMethod boundary,
                                             T background ){
    
    CoordinateVector padding;
    for(int i = 0; i < srcRes.vec.size()*2; i++)
      padding.vec.push_back(radius);
    
    fetchChannel(dstBuf, srcBuf, srcRes, padding, boundary,background );
    return;
    
  };
  
  /* copy a channel into a separate buffer.
     more general than fetchChannel, padding can be different in different
     dimensions */
  
  template<class T> void fetchChannel( T *dstBuf, T *srcBuf,
                                      CoordinateVector srcRes,
                                      CoordinateVector padding,
                                      BoundaryMethod boundary,
                                      T background );

   
  

} /* namespace */


#endif /* ARRAY_BOUNDARY_H */

