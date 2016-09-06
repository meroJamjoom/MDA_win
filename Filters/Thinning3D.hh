// ==========================================================================
// $Id:$
// A 3D morphological thinning operator
// Derived from the 2D thinning filter
// Based on the papers:
// * A Fully 3D Thinning Algorithm and Its Applications by Ma and Sonka
// * A note on 'A Fully 3D Thinning Algorithm and Its Applications' by Wang and Basu
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

#ifndef FILTERS_THINNING3D_H
#define FILTERS_THINNING3D_H

/*! \file  Thinning3D.hh
    \brief A 3D morphological thinning operator
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <set>
#include <cstring>

#include "MDA/Array/Array.hh"
#include "Filter.hh"

/* TODO:
 * Implement algorithm 2.5 which first marks boundary points and only then performs thinng.
 * Implement parallel version with jobs
 * Fix D template according to 'notes on...' - done in commented lines in D templates
*/
namespace MDA {

  //  typedef unsigned char ubyte;
    
  /** \class Thinning3D Thinning3D.hh
      A 3D morphological thinning operator */
  template<class T>
  class Thinning3D: public Filter<T> {

  public:

    /** constructor */
    Thinning3D();
    
    /** destructor */
    ~Thinning3D()
    {
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes );

  private:
    /*
      z   y
      |   /
      |  /
      | /
      |/_______x


      Order of the points in v array:
                                       
  o - stands for outer
  c - stands for center

              UNW--------UN---------UNE
            / |        / |        / |
          UW---------U----------UE  |
        / |   |    / |   |    / |   |
      USW--------US---------USE |   |
      |   |   NW-|---|---N--|---|---NE
      |   | / |  |   | / |  |   | / |
  ow -|---W------|---C------|---E----------oe
      | / |   |  | / |   |  | / |   |      |
      SW---------S----------SE  |   |      |
      |   |  DNW/|---|---DN-|---|---DNE    |
      |   | / os |   | /    |   | /        |
      |   DW-----|---D------|---DE---------doe
      | /        | / |      | / |          |
      DSW--------DS---------DSE |          |
                     |          |          |
		     |          |          |
                     od---------ode--------odoe

  the outer points above appear in other directions as well 
  (odow, ouw, osow, odos, ouos)


  OD plane (z=0)

      100------101------102------103------104
      |        |        |        |        |
      |        |        |        |        |
      |        |        |        |        |
      75-------76-------77-------78-------79
      |        |        |        |        |
      |        |        |        |        |
  odow|     odw|      od|     ode|    odoe|
      50-------51-------52-------53-------54
      |        |        |        |        |
      |        |        |        |        |
      |        |     ods|        |        |
      25-------26-------27-------28-------29
      |        |        |        |        |
      |        |        |        |        |
      |        |    odos|        |        |
      0--------1--------2--------3--------4


  D plane (z=1)

      105------106------107------108------109
      |        |        |        |        |
      |        |        |        |        |
      |     dnw|      dn|     dne|        |
      80-------81-------82-------83-------84
      |        |        |        |        |
      |        |        |        |        |
   dow|      dw|       d|      de|        |
      55-------56-------57-------58-------59
      |        |        |        |        |
      |        |        |        |        |
      |     dsw|      ds|     dse|        |
      30-------31-------32-------33-------34
      |        |        |        |        |
      |        |        |        |        |
      |        |     dos|        |        |
      5--------6--------7--------8--------9


  C plane (z=2)

      110------111------112------113------114
      |        |        |        |        |
      |        |        |        |        |
      |      nw|       n|      ne|        |
      85-------86-------87-------88-------89
      |        |        |        |        |
      |        |        |        |        |
    ow|       w|       c|       e|      oe|
      60-------61-------62-------63-------64
      |        |        |        |        |
      |        |        |        |        |
   sow|      sw|       s|      se|     soe|
      35-------36-------37-------38-------39
      |        |        |        |        |
      |        |        |        |        |
  osow|     osw|      os|     ose|    osoe|
      10-------11-------12-------13-------14


  U plane (z=3)

      115------116------117------118------119
      |        |        |        |        |
      |        |        |        |        |
      |     unw|      un|     une|        |
      90-------91-------92-------93-------94
      |        |        |        |        |
      |        |        |        |        |
   uow|      uw|       u|      ue|        |
      65-------66-------67-------68-------69
      |        |        |        |        |
      |        |        |        |        |
      |     usw|      us|     use|        |
      40-------41-------42-------43-------44
      |        |        |        |        |
      |        |        |        |        |
      |        |     uos|        |        |
      15-------16-------17-------18-------19


  OU plane (z=4)

      120------121------122------123------124
      |        |        |        |        |
      |        |        |        |        |
      |        |        |        |        |
      95-------96-------97-------98-------99
      |        |        |        |        |
      |        |        |        |        |
  ouow|     ouw|      ou|     oue|    ouoe|
      70-------71-------72-------73-------74
      |        |        |        |        |
      |        |        |        |        |
      |        |     ous|        |        |
      45-------46-------47-------48-------49
      |        |        |        |        |
      |        |        |        |        |
      |        |    ouos|        |        |
      20-------21-------22-------23-------24

    */

    struct Neighborhood
    {
      const T *p;    // original point reference
      T bg;          // background value
      Neighborhood( const T *_p, T _bg = T(0) ) : p( _p ), bg( _bg )
      {
      }

      // point aliases
      static const int Z = -1;  // none
      static const int C = 62;  // centre
      // 6-adjacent aliases
      static const int S = 37;  // south
      static const int N = 87;  // north
      static const int W = 61;  // west
      static const int E = 63;  // east
      static const int U = 67;  // up
      static const int D = 57;  // down
      // 18-adjacent aliases
      static const int SW = 36;
      static const int SE = 38;
      static const int NW = 86;
      static const int NE = 88;
      static const int DW = 56;
      static const int DE = 58;
      static const int UW = 66;
      static const int UE = 68;
      static const int DS = 32;
      static const int DN = 82;
      static const int US = 42;
      static const int UN = 92;
      // 26-adjacent aliases
      static const int DSW = 31;
      static const int DSE = 33;
      static const int DNW = 81;
      static const int DNE = 83;
      static const int USW = 41;
      static const int USE = 43;
      static const int UNW = 91;
      static const int UNE = 93;
      // outer aliases
      static const int OW = 60;  // outer west
      static const int OS = 12;  // outer south
      static const int OD = 52;  // outer down

      static const int OE = 64;
      static const int OU = 72;
      static const int OSE = 13;
      static const int SOE = 39;
      static const int OSOE = 14;
      static const int OSW = 11;
      static const int SOW = 35;
      static const int OSOW = 10;
      static const int OUS = 47;
      static const int UOS = 17;
      static const int OUOS = 22;
      static const int ODS = 27;
      static const int DOS = 7;
      static const int ODOS = 2;

      static const int OUW = 71;
      static const int UOW = 65;
      static const int OUOW = 70;
      static const int ODW = 51;
      static const int DOW = 55;
      static const int ODOW = 50;

      // true if point at index is an object point
      inline bool isObjectPoint( int index ) const     { return p[ index ] != bg; }
      // true if point at index is a background point
      inline bool isBackgroundPoint( int index ) const { return p[ index ] == bg; }

      int countNeighborObjectPoints( int stopafter = 26 ) const
      {
        static const int points[] = { DSW, DS, DSE, DW, D, DE, DNW, DN, DNE, 
                                       SW,  S,  SE,  W,     E,  NW,  N,  NE, 
                                      USW, US, USE, UW, U, UE, UNW, UN, UNE, Z };
        int count = 0;
        for ( int i = 0; points[ i ] != Z && count < stopafter; ++i )
        {
          if ( isObjectPoint( points[ i ] ) )
            count++;
        }
        return count;
      }

      // returns true if point p is 26-adjacent to a single object point
      // Rule 2.2 in paper
      bool isLineEndPoint() const
      {
        return isObjectPoint( C ) && countNeighborObjectPoints( 2 ) == 1;  // centre point and another one
      }

      // returns true if point p is 26-adjacent to two object points
      // Rule 2.2 in paper
      bool isNearLineEndPoint() const
      {
        if ( !isObjectPoint( C ) || countNeighborObjectPoints( 3 ) != 2 )  // centre point and two others
          return false;

        // no need to xor as that is taken by the above if
        return ( isObjectPoint( S ) && ( isObjectPoint( E ) || isObjectPoint( U ) ) ) ||
               ( isObjectPoint( W ) && ( isObjectPoint( N ) || isObjectPoint( U ) ) ) ||
               ( isObjectPoint( D ) && ( isObjectPoint( N ) || isObjectPoint( E ) ) );
      }

      // used for algorithm 2.5
      // currently not in use
      inline bool is26AdjacentToBackground() const
      {
        static const int points[] = { DSW, DS, DSE, DW, D, DE, DNW, DN, DNE, 
                                       SW,  S,  SE,  W,     E,  NW,  N,  NE, 
                                      USW, US, USE, UW, U, UE, UNW, UN, UNE, Z };
        return anyBackground( points );
      }

      inline bool isTailPoint() const
      {
        return isLineEndPoint() || isNearLineEndPoint();
      }

      // does pattern match any of the deletion templates?
      inline bool isDeleteTemplate() const
      {
         return isObjectPoint( C ) && ( isATemplate() || isBTemplate() || isCTemplate() || isDTemplate() );
      }

      inline bool isATemplate() const 
      {
	return matchA1() || matchA2() || matchA3() || matchA4() || matchA5() || matchA6();
      }
      
      inline bool isBTemplate() const 
      {
	return matchB1() || matchB2() || matchB3() || matchB4() || matchB5() || matchB6() ||
               matchB7() || matchB8() || matchB9() || matchB10() || matchB11() || matchB12();
      }

      inline bool isCTemplate() const 
      {
	return matchC1() || matchC2() || matchC3() || matchC4() || 
               matchC5() || matchC6() || matchC7() || matchC8();
      }
      
      inline bool isDTemplate() const 
      {
	return matchD1() || matchD2() || matchD3() || matchD4() || matchD5() || matchD6() ||
               matchD7() || matchD8() || matchD9() || matchD10() || matchD11() || matchD12();
      }

      // do all indices contain background points?
      bool allBackground( const int *pos ) const
      {
	while ( *pos != Z )
	{
	  if ( isBackgroundPoint( *pos++ ) == false )
	    return false;
	}
	return true;
      }

      // do any of the indices contain background points?
      bool anyBackground( const int *pos ) const
      {
	while ( *pos != Z )
	{
	  if ( isBackgroundPoint( *pos++ ) )
            return true;
        }
        return false;
      }

      // used by D template family for examining outer object points
      inline bool anyObjectPoint( int i1, int i2, int i3, int i4, int i5 ) const
      {
        return 
          isObjectPoint( i1 ) || 
          isObjectPoint( i2 ) || 
          isObjectPoint( i3 ) || 
          isObjectPoint( i4 ) || 
          isObjectPoint( i5 );
      }
 
      // A cases
      bool matchA1() const
      {
	static const int background[] = { DSW, SW, USW, DW, W, UW, DNW, NW, UNW, Z }; // west
	return isObjectPoint( E ) && 
          allBackground( background );
      }

      bool matchA2() const
      {
	static const int background[] = { DSE, SE, USE, DE, E, UE, DNE, NE, UNE, Z }; // east
	return isObjectPoint( W ) && 
          allBackground( background ) &&
          isObjectPoint( OW );
      }

      bool matchA3() const
      {
	static const int background[] = { DSW, DS, DSE, SW, S, SE, USW, US, USE, Z }; // south
	return isObjectPoint( N ) && 
          allBackground( background );
      }

      bool matchA4() const
      {
	static const int background[] = { DNW, DN, DNE, NW, N, NE, UNW, UN, UNE, Z }; // north
	return isObjectPoint( S ) && 
          allBackground( background ) &&
          isObjectPoint( OS );
      }

      bool matchA5() const
      {
	static const int background[] = { DSW, DS, DSE, DW, D, DE, DNW, DN, DNE ,Z }; // down
	return isObjectPoint( U ) && 
          allBackground( background );
      }

      bool matchA6() const
      {
	static const int background[] = { USW, US, USE, UW, U, UE, UNW, UN, UNE, Z }; // up
	return isObjectPoint( D ) && 
          allBackground( background ) &&
          isObjectPoint( OD );
      }

      // B cases
      bool matchB1() const
      {
	static const int background[] = { DW, W, UW, DNW, DN, NW, N, UNW, UN, Z }; // north-west
	return isObjectPoint( S ) && isObjectPoint( E ) && 
          allBackground( background ) &&
          isObjectPoint( OS );
      }

      bool matchB2() const
      {
	static const int background[] = { DE, E, UE, DN, DNE, N, NE, UN, UNE, Z }; // north-east
	return isObjectPoint( S ) && isObjectPoint( W ) && 
          allBackground( background ) &&
          isObjectPoint( OS ) && isObjectPoint( OW );
      }

      bool matchB3() const
      {
	static const int background[] = { DS, DSE, S, SE, US, USE, DE, E, UE, Z }; // south-east
	return isObjectPoint( W ) && isObjectPoint( N ) && 
          allBackground( background ) &&
          isObjectPoint( OW );
      }

      bool matchB4() const
      {
	static const int background[] = { DSW, DS, SW, S, USW, US, DW, W, UW, Z }; // south-west
	return isObjectPoint( E ) && isObjectPoint( N ) && 
          allBackground( background );
      }

      bool matchB5() const
      {
	static const int background[] = { DW, D, DE, DNW, DN, DNE, NW, N, NE, Z }; // down-north
	return isObjectPoint( S ) && isObjectPoint( U ) && 
          allBackground( background ) &&
          isObjectPoint( OS );
      }

      bool matchB6() const
      {
	static const int background[] = { DSW, DS, DSE, SW, S, SE, DW, D, DE, Z }; // down-south
	return isObjectPoint( U ) && isObjectPoint( N ) && 
          allBackground( background );
      }

      bool matchB7() const
      {
	static const int background[] = { UW, U, UE, NW, N, NE, UNW, UN, UNE, Z }; // up-north
	return isObjectPoint( S ) && isObjectPoint( D ) && 
          allBackground( background ) &&
          isObjectPoint( OS ) && isObjectPoint( OD );
      }

      bool matchB8() const
      {
	static const int background[] = { SW, S, SE, USW, US, USE, UW, U, UE, Z }; // up-south
	return isObjectPoint( D ) && isObjectPoint( N ) && 
          allBackground( background ) &&
          isObjectPoint( OD );
      }

      bool matchB9() const
      {
	static const int background[] = { DS, DSE, SE, D, DE, E, DN, DNE, NE, Z }; // down-east
	return isObjectPoint( W ) && isObjectPoint( U ) && 
          allBackground( background ) &&
          isObjectPoint( OW );
      }

      bool matchB10() const
      {
	static const int background[] = { DSW, DS, SW, DW, D, W, DNW, DN, NW, Z }; // down-west
	return isObjectPoint( E ) && isObjectPoint( U ) && 
          allBackground( background );
      }

      bool matchB11() const
      {
	static const int background[] = { SE, US, USE, E, U, UE, NE, UN, UNE, Z }; // up-east
	return isObjectPoint( W ) && isObjectPoint( D ) && 
          allBackground( background ) &&
          isObjectPoint( OW ) && isObjectPoint( OD );
      }

      bool matchB12() const
      {
	static const int background[] = { SW, USW, US, W, UW, U, NW, UNW, UN, Z }; // up-west
	return isObjectPoint( E ) && isObjectPoint( D ) && 
          allBackground( background ) &&
          isObjectPoint( OD );
      }

      // C cases
      bool matchC1() const
      {
	static const int background[] = { DS, DSE, S, SE, D, DE, E, Z }; // down-south-east
	return isObjectPoint( W ) && isObjectPoint( N ) && isObjectPoint( U ) &&
          allBackground( background ) &&
          isObjectPoint( OW );
      }

      bool matchC2() const
      {
	static const int background[] = { DSW, DS, SW, S, DW, D, W, Z }; // down-south-west
	return isObjectPoint( E ) && isObjectPoint( N ) && isObjectPoint( U ) &&
          allBackground( background );
      }

      bool matchC3() const
      {
	static const int background[] = { DW, D, W, DNW, DN, NW, N, Z }; // down-north-west
	return isObjectPoint( E ) && isObjectPoint( S ) && isObjectPoint( U ) &&
          allBackground( background ) &&
          isObjectPoint( OS );
      }

      bool matchC4() const
      {
	static const int background[] = { D, DE, E, DN, DNE, N, NE, Z }; // down-north-east
	return isObjectPoint( W ) && isObjectPoint( S ) && isObjectPoint( U ) &&
          allBackground( background ) &&
          isObjectPoint( OW ) && isObjectPoint( OS );
      }

      bool matchC5() const
      {
	static const int background[] = { S, SE, US, USE, E, U, UE, Z }; // up-south-east
	return isObjectPoint( W ) && isObjectPoint( N ) && isObjectPoint( D ) &&
          allBackground( background ) &&
          isObjectPoint( OW ) && isObjectPoint( OD );
      }

      bool matchC6() const
      {
	static const int background[] = { SW, S, USW, US, W, UW, U, Z }; // up-south-west
	return isObjectPoint( E ) && isObjectPoint( N ) && isObjectPoint( D ) &&
          allBackground( background ) &&
          isObjectPoint( OD );
      }

      bool matchC7() const
      {
	static const int background[] = { W, UW, U, NW, N, UNW, UN, Z }; // up-north-west
	return isObjectPoint( E ) && isObjectPoint( S ) && isObjectPoint( D ) &&
          allBackground( background ) &&
          isObjectPoint( OS ) && isObjectPoint( OD );
      }

      bool matchC8() const
      {
	static const int background[] = { E, U, UE, N, NE, UN, UNE, Z }; // up-north-east
	return isObjectPoint( W ) && isObjectPoint( S ) && isObjectPoint( D ) &&
          allBackground( background ) &&
          isObjectPoint( OW ) && isObjectPoint( OS ) && isObjectPoint( OD );
      }

      // D cases 
      // commented out code has the original paper implementation
      // the non commented line before the comment is from the 'notes on...' paper
      bool matchD1() const
      {
	static const int background[] = { DSW, DS, DSE, S, DW, D, DE, U, DNW, DN, DNE, NW, N, NE, UNW, UN, UNE, Z }; // down-north
	return isObjectPoint( US ) && !( isObjectPoint( W ) && isObjectPoint( E ) ) &&
          //return isObjectPoint( US ) && 
          allBackground( background ) &&
          anyObjectPoint( OU, OUS, OS, UOS, OUOS );
      }

      bool matchD2() const
      {
	static const int background[] = { S, USW, US, USE, D, UW, U, UE, DNW, DN, DNE, NW, N, NE, UNW, UN, UNE, Z }; // up-north
	return isObjectPoint( DS ) && !( isObjectPoint( W ) && isObjectPoint( E ) ) &&
          //return isObjectPoint( DS ) && 
          allBackground( background ) &&
          anyObjectPoint( OD, ODS, OS, DOS, ODOS );
      }

      bool matchD3() const
      {
	static const int background[] = { S, DSE, SE, USE, W, DE, E, UE, DNW, DN, DNE, NW, N, NE, UNW, UN, UNE, Z }; // north-east
	return isObjectPoint( SW ) && !( isObjectPoint( D ) && isObjectPoint( U ) ) &&
          //return isObjectPoint( SW ) && 
          allBackground( background ) &&
          anyObjectPoint( OS, OSW, OW, SOW, OSOW );
      }

      bool matchD4() const
      {
	static const int background[] = { DSW, SW, S, USW, DW, W, E, UW, DNW, DN, DNE, NW, N, NE, UNW, UN, UNE, Z }; // north-west
	return isObjectPoint( SE ) && !( isObjectPoint( D ) && isObjectPoint( U ) ) &&
          //return isObjectPoint( SE ) && 
          allBackground( background ) &&
          anyObjectPoint( OS, OSE, OE, SOE, OSOE );
      }

      bool matchD5() const
      {
	static const int background[] = { DSW, DS, DSE, SE, USE, DW, D, DE, W, E, U, UE, DNW, DN, DNE, NE, UNE, Z }; // down-east
	return isObjectPoint( UW ) && !( isObjectPoint( S ) && isObjectPoint( N ) ) &&
          //return isObjectPoint( UW ) && 
          allBackground( background ) &&
          anyObjectPoint( OW, UOW, OU, OUW, OUOW );
      }

      bool matchD6() const
      {
	static const int background[] = { DSE, SE, USW, US, USE, D, DE, W, E, UW, U, UE, DNE, NE, UNW, UN, UNE, Z }; // up-east
	return isObjectPoint( DW ) && !( isObjectPoint( S ) && isObjectPoint( N ) ) &&
          //return isObjectPoint( DW ) && 
          allBackground( background ) &&
          anyObjectPoint( OW, DOW, OD, ODW, ODOW );
      }

        // THESE ARE INCOMPLETE
      bool matchD7() const
      {
	static const int background[] = { SW, S, SE, W, E, NW, N, Z }; // 
	return isObjectPoint( NE ) && ( matchD7_1() || matchD7_2() || matchD7_3() ) &&
          //return isObjectPoint( NE ) && isObjectPoint( D ) && isObjectPoint( U ) &&
          allBackground( background );
      }

      bool matchD7_1() const
      {
        static const int background[] = { DSW, DS, DSE, DW, D, DNW, USW, US, USE, UW, U, UNW, Z }; // 
        return 
          allBackground( background );
      }

      bool matchD7_2() const
      {
        static const int background[] = { USW, US, USE, UW, U, UNW, Z }; // 
        return isObjectPoint( D ) &&
          allBackground( background );
      }

      bool matchD7_3() const
      {
        static const int background[] = { DSW, DS, DSE, DW, D, DNW, Z }; // 
        return isObjectPoint( U ) &&
          allBackground( background );
      }

      bool matchD8() const
      {
	static const int background[] = { SW, S, SE, W, E, N, NE, Z }; // 
	return isObjectPoint( NW ) && ( matchD8_1() || matchD8_2() || matchD8_3() ) &&
          //return isObjectPoint( NW ) && isObjectPoint( D ) && isObjectPoint( U ) &&
          allBackground( background );
      }

      bool matchD8_1() const
      {
        static const int background[] = { DSW, DS, DSE, D, DE, DNE, USW, US, USE, U, UE, UNE, Z }; // 
        return 
          allBackground( background );
      }

      bool matchD8_2() const
      {
        static const int background[] = { USW, US, USE, U, UE, UNE, Z }; // 
        return isObjectPoint( D ) && 
          allBackground( background );
      }

      bool matchD8_3() const
      {
        static const int background[] = { DSW, DS, DSE, D, DE, DNE, Z }; // 
        return isObjectPoint( U ) && 
          allBackground( background );
      }

      bool matchD9() const
      {
	static const int background[] = { DS, S, US, D, U, DN, N, Z }; //
	return isObjectPoint( UN ) && ( matchD9_1() || matchD9_2() || matchD9_3() ) &&
          //return isObjectPoint( UN ) && isObjectPoint( W ) && isObjectPoint( E ) &&
          allBackground( background );
      }

      bool matchD9_1() const
      {
        static const int background[] = { DSW, DW, DNW, SW, W, USW, DSE, DE, DNE, SE, E, USE, Z }; //
	return 
          allBackground( background );
      }

      bool matchD9_2() const
      {
        static const int background[] = { DSW, DW, DNW, SW, W, USW, Z }; //
	return isObjectPoint( E ) && 
          allBackground( background );
      }

      bool matchD9_3() const
      {
        static const int background[] = { DSE, DE, DNE, SE, E, USE, Z }; //
	return isObjectPoint( W ) && 
          allBackground( background );
      }

      bool matchD10() const
      {
	static const int background[] = { DS, S, US, D, U, N, UN, Z }; //
	return isObjectPoint( DN ) && ( matchD10_1() || matchD10_2() || matchD10_3() ) &&
          //return isObjectPoint( DN ) && isObjectPoint( W ) && isObjectPoint( E ) &&
          allBackground( background );
      }

      bool matchD10_1() const
      {
        static const int background[] = { DSW, SW, W, USW, UW, UNW, DSE, SE, E, USE, UE, UNE, Z }; //
	return
          allBackground( background );
      }

      bool matchD10_2() const
      {
        static const int background[] = { DSW, SW, W, USW, UW, UNW, Z }; //
	return isObjectPoint( E ) &&
          allBackground( background );
      }

      bool matchD10_3() const
      {
        static const int background[] = { DSE, SE, E, USE, UE, UNE, Z }; //
	return isObjectPoint( W ) &&
          allBackground( background );
      }

      bool matchD11() const
      {
	static const int background[] = { DW, D, W, E, UW, U, UE, Z }; //
	return isObjectPoint( DE ) && ( matchD11_1() || matchD11_2() || matchD11_3() ) &&
          //return isObjectPoint( DE ) && isObjectPoint( S ) && isObjectPoint( N ) &&
          allBackground( background );
      }

      bool matchD11_1() const
      {
        static const int background[] = { DSW, SW, S, USW, US, USE, DNW, NW, N, UNW, UN, UNE, Z }; //
	return
          allBackground( background );
      }

      bool matchD11_2() const
      {
        static const int background[] = { DNW, NW, N, UNW, UN, UNE, Z }; //
	return isObjectPoint( S ) &&
          allBackground( background );
      }

      bool matchD11_3() const
      {
        static const int background[] = { DSW, SW, S, USW, US, USE, Z }; //
	return isObjectPoint( N ) &&
          allBackground( background );
      }

      bool matchD12() const
      {
	static const int background[] = { DW, D, DE, W, E, UW, U, Z }; //
        return isObjectPoint( UE ) && ( matchD12_1() || matchD12_2() || matchD12_3() ) &&
          //return isObjectPoint( UE ) && isObjectPoint( S ) && isObjectPoint( N ) &&
          allBackground( background );
      }

      bool matchD12_1() const
      {
        static const int background[] = { DSW, DS, DSE, SW, S, USW, DNW, DN, DNE, NW, N, UNW, Z }; //
        return 
          allBackground( background );
      }

      bool matchD12_2() const
      {
        static const int background[] = { DNW, DN, DNE, NW, N, UNW, Z }; //
        return isObjectPoint( S ) && 
          allBackground( background );
      }

      bool matchD12_3() const
      {
        static const int background[] = { DSW, DS, DSE, SW, S, USW, Z }; //
        return isObjectPoint( N ) && 
          allBackground( background );
      }

    }; //neighborhood

  }; //thinning


  /** \class Thickening2D Thinning2D.hh
      A 3D morphological thickening operator (dual of thinning) */
  template<class T>
  class Thickening3D: public Filter<T> {

  public:

    /** constructor */
    Thickening3D()
    {
    }
    
    /** destructor */
    ~Thickening3D()
    {
    }
    
    /** apply the filter to a number of dimensions and channels */
    virtual bool apply( Array<T> &a, BoundaryMethod boundary,
			ChannelList &channels, AxisList &axes )
    {
      return false;
    }

    
  private:

  };


} /* namespace */

#endif /* FILTERS_THINNING3D_H */

