// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2012-, UBC
// 
// Creator: krim ()
// Email:   krim@cs.ubc.ca
// ==========================================================================

#ifndef FILTERS_THINNING3D_C
#define FILTERS_THINNING3D_C

#include "Thinning3D.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/Array.hh"

#include <cassert>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  /*
   * Copies a block of a 3d array (src), into destination array (dst) 
   * setting all out of bounds points to background (bg).
   * Destination array should be of size xnum * ynum * znum.
   */
template<typename T>
void copySubArray3D( T *dst, const T *src,
                     int width, int height, int depth,
                     int xs, int ys, int zs,
                     int xnum, int ynum, int znum,
                     T bg = T(0) )
{
  int i = 0;
  for ( int z = zs; z < zs + znum; z++ )
  {
    bool outer_z = z < 0 || z >= depth;
    for ( int y = ys; y < ys + ynum; y++ )
    {
      bool outer_y = outer_z || y < 0 || y >= height;
      for ( int x = xs; x < xs + xnum; x++ )
      {
        bool outer = outer_y || x < 0 || x >= width;
        dst[ i++ ] = outer ? bg : src[ x + width * ( y + z * height ) ];
      }
    }
  }
}


                   
//
// Thinning3D code
//

/** constructor */
template <class T>
Thinning3D<T>::Thinning3D()
{
}

/** apply the filter to a number of dimensions and channels */
template <class T>
bool
Thinning3D<T>::apply( Array<T> &a, BoundaryMethod boundary,
		      ChannelList &channels, AxisList &axes )
{
  unsigned long i, j, k;

  CoordinateVector dim = a.getDimension();
  if ( !warnCond( dim.vec.size() == 3 && axes.vec.size() == 3,
                  "  only works for 3D arrays with both axes active\n" ) )
    return false;

  unsigned long w= dim.vec[0];
  unsigned long h= dim.vec[1];
  unsigned long d= dim.vec[2];

  T neighbors[ 5*5*5 ];

  // process all channels
  for( int ch= 0 ; ch< channels.vec.size() ; ch++ )
  {
    // grab the full channel, and store it as a padded bool array
    T *data= &((*a[channels.vec[ch]])[0]);
    T bg= a[ch]->getBackground();
    T *pt = data;  // current point

    long deletedPoints = 0;
    do
    {
      pt = data;
      deletedPoints = 0;
      for ( k=0; k<d; k++ )
      {
        for ( j=0; j<h; j++ )
        {
          for ( i=0; i<w; i++, pt++ )
          {
            if ( *pt != bg )  // is non background point?
            {
              copySubArray3D( neighbors, data,
                              w, h, d,
                              i-2, j-2, k-2,
                              5, 5, 5,
                              bg );
              Neighborhood hood( neighbors, bg );
              if ( !hood.isTailPoint() && hood.isDeleteTemplate() )
              {
                *pt = bg;   //data[ i + w * ( j + h * k ) ] = bg;
                deletedPoints++;
              }
            }
          }
        }
      }
    } while ( deletedPoints > 0 );

  } // for ch

  return true;
}


//
// Thickening3D code
//

// TO BE IMPLEMENTED



// explicit template instation code
template class Thinning3D<unsigned char>;
template class Thinning3D<char>;
template class Thinning3D<unsigned short>;
template class Thinning3D<short>;
template class Thinning3D<unsigned int>;
template class Thinning3D<int>;
template class Thinning3D<float>;
template class Thinning3D<double>;
//template class Thickening3D<float>;
//template class Thickening3D<double>;



} /* namespace */

#endif /* FILTERS_THINNING3D_C */

