// ==========================================================================
// $Id: HitOrMiss2D.C 372 2009-09-20 18:57:38Z heidrich $
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

#ifndef FILTERS_HITORMISS2D_C
#define FILTERS_HITORMISS2D_C

#include "MDA/Base/Errors.hh"
#include "MDA/Threading/SMPJobManager.hh"

#include "HitOrMiss2D.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/* the indices in the structuring element and mask are provided by
 * the usual column-first order used everywhere in MDA arrays:
 *
 *	0 1 2
 *	3 4 5
 *	6 7 8
 *
 * However, the bits for indices into the case table are allocated to
 * pixels in a neighborhood according to the transpose of the above
 * pattern:
 *
 *	0 3 6
 *	1 4 7
 *	2 5 8
 *
 * This is convenient for updating the neighborhood along a scanline
 * from left to right: if index[n] is the index for pixel n, then
 * index[n+1]= (index[n]>>3 | newBits), where newBits are then new
 * bits 6-8, representing the next column over.
 */
unsigned bitMapping[9] = {0,3,6,1,4,7,2,5,8};
    

/** constructor from a more traditional structural element (&mask) */
template <class T>
HitOrMiss2D<T>::HitOrMiss2D( bool structureElem[9], bool mask[9],
			     bool rot, bool refl, bool _preserveMisses,
			     T _hitValue, T _missValue )
  : preserveMisses( _preserveMisses ),
    hitValue( _hitValue ), missValue( _missValue )
{
  caseTable= new double[512];
  unsigned i, j, k;
  
  unsigned maskBits= neighborhoodToIndex( mask );
  unsigned structureBits= neighborhoodToIndex( structureElem ) & maskBits;
  
  for( i= 0 ; i< 512 ; i++ )
  {
    caseTable[i]= missValue;
    
    // check all rotations and reflection for a match with the
    // structuring element
    for( j= i, k= 0 ; k< (rot ? 4 : 1) ; k++ )
      if( (j & maskBits)== structureBits ||
	  (refl && ((reflect( j ) & maskBits)== structureBits)) )
      {
	caseTable[i]= hitValue;
	break;
      }
      else
	j= rotate( j );
  }
}


/** apply the filter to a number of dimensions and channels */
template <class T>
bool
HitOrMiss2D<T>::apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )
{
  unsigned long i, j, k;
  
  CoordinateVector dim= a.getDimension();
  if( !warnCond( dim.vec.size()== 2 && axes.vec.size()== 2,
		 "  only works for 2D arrays with both axes active\n" ) )
    return false;
  
  // data for one line in original data type, and a full channel in bool
  unsigned long w= dim.vec[0];
  unsigned long h= dim.vec[1];
  unsigned long wp= w+2;
  unsigned long hp= h+2;
  bool *boolData= new bool[wp*hp];
  T *lineData= new T[wp];
  
  // process all channels
  for( k= 0 ; k< channels.vec.size() ; k++ )
  {
    // grab the full channel, and store it as a padded bool array
    T *data= &((*a[channels.vec[k]])[0]);
    T bg= a[k]->getBackground();
    for( i= 0 ; i< h ; i++ )
    {
      // grab one line, translate to bool
      fetchLine<T>( lineData, data+i*w, 1, w, 1, boundary, bg );
      for( j= 0 ; j< wp ; j++ )
	boolData[(i+1)*wp+j]= (lineData[j]> 0.0);
    }
    // take care of top and bottom lines
    memcpy( boolData, boolData+wp, wp*sizeof( bool ) );
    memcpy( boolData+(hp-1)*wp, boolData+(hp-2)*wp, wp*sizeof( bool ) );
    
    // execute line jobs in parallel
    SMPJobList jobs;
    for( i= 0 ; i< h ; i++ )
      jobs.push_back( new HitOrMissLineJob<T>( (HitOrMiss2D<T> *)this,
					       boolData+(i+1)*wp,
					       data+i*w, w ) );
    SMPJobManager::getJobManager()->batch( jobs );
  }
  
  return true;
}


/** apply to a single scanline */
template <class T>
void
HitOrMiss2D<T>::apply( bool *currLine, T* dst, unsigned long numElements )
{
  bool *prevLine= currLine-(numElements+2);
  bool *nextLine= currLine+(numElements+2);
  
  // initialize table index to first two columns
  unsigned index= 0;
  index|= (*(prevLine++) ?   8 : 0);
  index|= (*(currLine++) ?  16 : 0);
  index|= (*(nextLine++) ?  32 : 0);
  index|= (*(prevLine++) ?  64 : 0);
  index|= (*(currLine++) ? 128 : 0);
  index|= (*(nextLine++) ? 256 : 0);
  
  // go over the rest of the line, update index, and determine hit or miss
  for( unsigned long i= 0 ; i< numElements ; i++ )
  {
    // update index
    index>>= 3;
    index|= (*(prevLine++) ?  64 : 0);
    index|= (*(currLine++) ? 128 : 0);
    index|= (*(nextLine++) ? 256 : 0);
    
    // match result
    // if result equals the missValue, and we preserve misses, then just
    // keep the current pixel value
    if( !preserveMisses || caseTable[index]!= missValue )
      dst[i]= caseTable[index];
  }
}

/** convert a neighborhood bit vector to a table index */
template <class T>
unsigned
HitOrMiss2D<T>::neighborhoodToIndex( const bool neighborhood[9] )
{
  unsigned index= 0;
  for( unsigned i= 0 ; i< 9 ; i++ )
    index|= (neighborhood[i] << bitMapping[i]);
  return index;
}
  
/** convert a table index to a neighborhood bit vector */
template <class T>
void
HitOrMiss2D<T>::indexToNeighborhood( unsigned index, bool neighborhood[9] )
{
  for( unsigned i= 0 ; i< 9 ; i++ )
    neighborhood[i]= index & (1 << bitMapping[i]);
}   

/** rotate a table index by 90 degrees ccw */
template <class T>
unsigned
HitOrMiss2D<T>::rotate( unsigned ind )
{
  bool inNeighbor[9];
  bool outNeighbor[9];
  
  indexToNeighborhood( ind, inNeighbor );
  outNeighbor[0]= inNeighbor[2];
  outNeighbor[1]= inNeighbor[5];
  outNeighbor[2]= inNeighbor[8];
  outNeighbor[3]= inNeighbor[1];
  outNeighbor[4]= inNeighbor[4];
  outNeighbor[5]= inNeighbor[7];
  outNeighbor[6]= inNeighbor[0];
  outNeighbor[7]= inNeighbor[3];
  outNeighbor[8]= inNeighbor[6];
  
  return neighborhoodToIndex( outNeighbor );
}
  
/** reflect a table index horizontally */
template <class T>
unsigned
HitOrMiss2D<T>::reflect( unsigned ind )
{
  unsigned left= ind & 7;
  unsigned center= ind & 56;
  unsigned right= ind & 448;
  
  return (left << 6) | center | (right >> 6);
}
    



/** execute line job */
template <class T>
void
HitOrMissLineJob<T>::execute( int jobID )
{
  filter->apply( currLine, dst, numElements );
}


    
// explicit template instation code

template class HitOrMiss2D<float>;
template class HitOrMiss2D<double>;
template class HitOrMissLineJob<float>;
template class HitOrMissLineJob<double>;


} /* namespace */

#endif /* FILTERS_HITORMISS2D_C */

