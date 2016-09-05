// ==========================================================================
// $Id: Boundary.C 423 2009-10-16 18:18:15Z skoch $
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

#ifndef ARRAY_BOUNDARY_C
#define ARRAY_BOUNDARY_C

#include <string.h>

#include "MDA/Base/Errors.hh"
#include "Boundary.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

static const int maxBoundaryNameLength= 12;
static const char *boundaryModeNames[5]=
{
  "background", "clamp", "cyclic", "mirror", "renormalize"
};


/** write boundary mode name to ostream */
ostream &
operator<<( ostream &os, BoundaryMethod boundary )
{
  return os << boundaryModeNames[boundary];
}

/** read boundary mode name from istream */
istream &
operator>>( istream &is, BoundaryMethod &boundary )
{
  int i;
  char buffer[maxBoundaryNameLength];

  is >> ws;
  for( i= 0 ; i< maxBoundaryNameLength && isalpha( is.peek() ) ; i++ )
    buffer[i]= is.get();
  buffer[i]= '\0';

  for( i= 0 ; i< 5 ; i++ )
    if( !strncasecmp( buffer, boundaryModeNames[i], maxBoundaryNameLength ) )
    {
      boundary= (BoundaryMethod)i;
      return is;
    }

  // if we get here, things have gone wrong
  warning( "  could not read boundary mode\n" );
  is.setstate( ios::failbit );
  return is;
}


  //
  // BoundaryOption methods
  //

/** the actual parsing function */
bool
BoundaryOption::parse( int &index, int argc, char *argv[] )
{
  if( index>  argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> boundary;
  return !optStr.fail();
}

/** output usage string */
void
BoundaryOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <mode>\n"
     << helpTxt
     << "\tCurrently: " << boundary << "\n\n";
}


//
// functions for getting boundary padded arrays and lines
//

  /** pad data line according to boundary mode */
  template<class T>
  void
  fetchLine( T *dstBuf, T *srcBuf, unsigned long srcStride, 
            unsigned long numElements, unsigned long left,unsigned long right,
            BoundaryMethod boundary, T background )
  {
    unsigned long i, j;
    T *tmp;
    // copy the interior data first
    if( srcStride== 1 )
      memcpy( dstBuf+left, srcBuf, sizeof(T)*numElements );
    else
      for( i= numElements, tmp= dstBuf+right ;
          i> 0 ;
          i--, srcBuf+= srcStride )
        *tmp++= *srcBuf;

    // then handle the boundaries
    switch( boundary )
    {
      case Background:
        // all values outside the domain are zero
        for( i= 0 ; i< left ; i++ )
          dstBuf[i]= background;
        for( i= 0 ; i< right ; i++ )
          dstBuf[numElements+left+i] = background;
        break;
      case Clamp:
        // outside values are clamped to the value of the boundary 
        for( i= 0 ; i< left ; i++ )
          dstBuf[i]= dstBuf[left];
        for( i= 0 ; i< right ; i++ )
          dstBuf[numElements+left+i]= dstBuf[numElements+left-1];
        break;
      case Cyclic:
        // periodic tiling
        
        for( i= 0 ; i< left ; i++ )
        {
          // Note that the order in which the entries are written accounts
          // for multiple wrap-arounds
          dstBuf[left-1-i]= dstBuf[numElements+left-1-i];
        }
        
        for (i=0; i < right; i++) 
          dstBuf[numElements+left+i] = dstBuf[left+i];
        
        break;
      case Mirror:
        // the outside is a mirror reflection of the inside
        for( i= 0 ; i< left ; i++ )
        {
          // Note that the order in which the entries are written accounts
          // for multiple wrap-arounds
          dstBuf[left-1-i]= dstBuf[left+i];
        }
        for( i= 0 ; i< right ; i++ )
          dstBuf[numElements+left+i]= dstBuf[numElements+left-1-i];
        
        break;
      case Renormalize:
        // renormalize involves altering a posible filter function, so
        // nothing can be done here.
        break;
    }
  }
  
inline bool
inRange( const CoordinateVector &pos, const LongRangeList &box )
{
  for( int i= pos.vec.size()-1 ; i> 0 ; i-- )
    if( pos.vec[i]< box.vec[i].val.first ||
        pos.vec[i]> box.vec[i].val.second )
      return false;
  return true;
}
  
  
  /* copy a channel into a separate buffer.
   more general than fetchChannel, padding can be different in different
   dimensions */
  
template<class T> 
void 
fetchChannel( T *dstBuf, T *srcBuf, CoordinateVector srcRes,
                     CoordinateVector padding, BoundaryMethod boundary,
                     T background ){
  unsigned long i;
  
  // dimension and array increments for the different axes
  unsigned dimension= srcRes.vec.size();
  
  LongRangeList srcBox( dimension );
  LongRangeList dstBox( dimension );
  
  for( i= 0 ; i< dimension ; i++ )
  {
    srcBox.vec[i].val.first= 0l;
    srcBox.vec[i].val.second= srcRes.vec[i]-1;
    dstBox.vec[i].val.first= -padding.vec[i*2];
    dstBox.vec[i].val.second= srcRes.vec[i]+padding.vec[i*2+1]-1;
  }
  
  CoordinateVector srcIncr( dimension, 1 );
  CoordinateVector dstIncr( dimension, 1 );
  for( i= 1 ; i< dimension ; i++ )
  {
    srcIncr.vec[i]= srcIncr.vec[i-1] * srcRes.vec[i-1];
    dstIncr.vec[i]= dstIncr.vec[i-1] * (srcRes.vec[i-1]
                                        +padding.vec[(i-1)*2]
                                        +padding.vec[(i-1)*2+1]);
  }
  
  //
  // first, copy original scanlines, already boundary padded in the
  // direction of axis 0
  //  
  T *srcPtr= srcBuf;
  T *dstPtr= dstBuf;
  int step = 0;
  CoordinateVectorIter iter( dstBox );
  for( iter.begin() ; !iter.isAtEnd() ; iter.incrComp( 1 ) )
  {
    if( inRange( iter.getPos(), srcBox ) )
    {
      fetchLine( dstPtr, srcPtr, 1, srcRes.vec[0], padding.vec[0], 
                       padding.vec[1], boundary, background );
      srcPtr+= srcIncr.vec[1];
    }
    dstPtr+= dstIncr.vec[1];
  }
  
  //
  // second, pad the boundary scanlines
  //
  T *backgroundLine= NULL;
  if( boundary== Background )
  {
    backgroundLine= new T[dstIncr.vec[1]];
    for( i= 0 ; i< dstIncr.vec[1] ; i++ )
      backgroundLine[i]= background;
  }
  for( dstPtr= dstBuf, iter.begin() ; !iter.isAtEnd() ; iter.incrComp( 1 ) )
  {
    CoordinateVector pos= iter.getPos();
    if( !inRange( pos, srcBox ) )
      switch( boundary )
    {
      case Background:
        // pad with background color
        memcpy( dstPtr, backgroundLine, dstIncr.vec[1]*sizeof(T) );
        break;
      case Clamp:
        // find closest scanline on boundary
        srcPtr= dstBuf;
        for( i= 1 ; i< dimension ; i++ )
          if( pos.vec[i]>= srcRes.vec[i] )
            srcPtr+= (srcRes.vec[i]+padding.vec[i*2]-1) * dstIncr.vec[i];
          else if( pos.vec[i]>= 0 )
            srcPtr+= (pos.vec[i]+padding.vec[i*2]) * dstIncr.vec[i];
          else
            srcPtr+= padding.vec[i*2] * dstIncr.vec[i];
        memcpy( dstPtr, srcPtr, dstIncr.vec[1]*sizeof(T) );
        break;
      case Cyclic:
        // wrap around cyclically in all directions
        srcPtr= dstBuf;
        for( i= 1 ; i< dimension ; i++ )
        {
          int tmp= pos.vec[i] % srcRes.vec[i];
          if( tmp< 0 )
            tmp+= srcRes.vec[i];
          srcPtr+= (tmp+padding.vec[i*2]) * dstIncr.vec[i];
        }
        memcpy( dstPtr, srcPtr, dstIncr.vec[1]*sizeof(T) );
        break;
      case Mirror:
        // reflect along all boundaries
        // wrap around cyclically in all directions
        srcPtr= dstBuf;
        for( i= 1 ; i< dimension ; i++ )
        {
          int tmp= pos.vec[i] % (2*srcRes.vec[i]);
          if( tmp>= srcRes.vec[i] )
            tmp= 2*srcRes.vec[i] - tmp-1;
          else if( tmp< -srcRes.vec[i] )
            tmp+= 2*srcRes.vec[i];
          else if( tmp< 0 )
            tmp= -tmp-1;
          srcPtr+= (tmp+padding.vec[i*2]) * dstIncr.vec[i];
        }
        memcpy( dstPtr, srcPtr, dstIncr.vec[1]*sizeof(T) );
        break;
      default: // can't do anything for renormalize mode
        break;
    }
    
    dstPtr+= dstIncr.vec[1];
  }
  
  //
  // clean up
  //
  if( backgroundLine!= NULL )
    delete [] backgroundLine;
  
}
  
  



//
// template instantiation code
//
template void fetchLine<float>( float *dstBuf, float *srcBuf,
				unsigned long srcStride,
				unsigned long numElements,
				unsigned long radius,
				BoundaryMethod boundary, float background );
template void fetchLine<double>( double *dstBuf, double *srcBuf,
				 unsigned long srcStride,
				 unsigned long numElements,
				 unsigned long radius,
				 BoundaryMethod boundary, double background );
template void fetchChannel<float>( float *dstBuf, float *srcBuf,
				   CoordinateVector srcRes,
				   unsigned long radius,
				   BoundaryMethod boundary,
				   float background );
template void fetchChannel<double>( double *dstBuf, double *srcBuf,
				    CoordinateVector srcRes,
				    unsigned long radius,
				    BoundaryMethod boundary,
				    double background );

template void fetchLine<float>( float *dstBuf, float *srcBuf,
                               unsigned long srcStride,
                               unsigned long numElements,
                               unsigned long paddingleft, 
                               unsigned long paddingright,
                               BoundaryMethod boundary, float background );
template void fetchLine<double>( double *dstBuf, double *srcBuf,
                                unsigned long srcStride,
                                unsigned long numElements,
                                unsigned long paddingleft, 
                                unsigned long paddingright,
                                BoundaryMethod boundary, double background );
template void fetchChannel<float>( float *dstBuf, float *srcBuf,
                                  CoordinateVector srcRes,
                                  CoordinateVector padding,
                                  BoundaryMethod boundary,
                                  float background );
template void fetchChannel<double>( double *dstBuf, double *srcBuf,
                                   CoordinateVector srcRes,
                                   CoordinateVector padding,
                                   BoundaryMethod boundary,
                                   double background );
  
  

} /* namespace */

#endif /* ARRAY_BOUNDARY_C */

