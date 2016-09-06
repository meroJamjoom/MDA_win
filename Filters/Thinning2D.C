// ==========================================================================
// $Id: Thinning2D.C 374 2009-09-21 21:13:13Z heidrich $
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

#ifndef FILTERS_THINNING2D_C
#define FILTERS_THINNING2D_C

#include "Thinning2D.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

bool strucElem[8][9]= {
  {0,0,0,0,1,0,1,1,1},
  {0,0,0,1,1,0,0,1,0},
  {1,1,1,0,1,0,0,0,0},
  {0,1,0,1,1,0,0,0,0},
  {0,0,1,0,1,1,0,0,1},
  {0,1,0,0,1,1,0,0,0},
  {1,0,0,1,1,0,1,0,0},
  {0,0,0,0,1,1,0,1,0}
};

bool mask[8][9]= {
  {1,1,1,0,1,0,1,1,1},
  {0,1,1,1,1,1,0,1,0},
  {1,1,1,0,1,0,1,1,1},
  {0,1,0,1,1,1,0,1,1},
  {1,0,1,1,1,1,1,0,1},
  {0,1,0,1,1,1,1,1,0},
  {1,0,1,1,1,1,1,0,1},
  {1,1,0,1,1,1,0,1,0}
};

//
// Thinning2D code
//

/** constructor */
template <class T>
Thinning2D<T>::Thinning2D()
{
  for( unsigned i= 0 ; i< 8 ; i++ )
    hmt[i]= new HitOrMiss2D<T>( strucElem[i], mask[i],
				false, false, true, 0.0, 1.0 );
}


/** apply the filter to a number of dimensions and channels */
template <class T>
bool
Thinning2D<T>::apply( Array<T> &a, BoundaryMethod boundary,
		      ChannelList &channels, AxisList &axes )
{
  for( unsigned i= 0 ; i< 8 ; i++ )
    hmt[i]->apply( a, boundary, channels, axes );
  
  return true;
}


//
// Thickening2D code
//

/** constructor */
template <class T>
Thickening2D<T>::Thickening2D()
{
  bool sElem[9];
  
  for( unsigned i= 0 ; i< 8 ; i++ )
  {
    // negate mask
    for( unsigned j= 0 ; j< 9 ; j++ )
      sElem[j]= !(strucElem[i][j]);
    hmt[i]= new HitOrMiss2D<T>( sElem, mask[i], false, false, true );
  }
}


/** apply the filter to a number of dimensions and channels */
template <class T>
bool
Thickening2D<T>::apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )
{
  for( unsigned i= 0 ; i< 8 ; i++ )
    hmt[i]->apply( a, boundary, channels, axes );
  
  return true;
}



// explicit template instation code
template class Thinning2D<float>;
template class Thinning2D<double>;
template class Thickening2D<float>;
template class Thickening2D<double>;



} /* namespace */

#endif /* FILTERS_THINNING2D_C */

