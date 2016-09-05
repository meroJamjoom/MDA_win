// ==========================================================================
// $Id: CoordinateVector.hh 413 2009-10-15 10:47:40Z heidrich $
// Index vector into a multi-dimensional array
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2006-, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef BASE_COORDINATEVECTOR_H
#define BASE_COORDINATEVECTOR_H

/*! \file  CoordinateVector.hh
    \brief Index vector into a multi-dimensional array
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>

#include "Range.hh"
#include "CommandlineParser.hh"


namespace MDA {
  
  using namespace std;
  
  
  /** \class CoordinateVector CoordinateVector.hh
      a coordinate index vector into a multi-dimensional array */
  class CoordinateVector {
  public:
    /** constructor for empty coordinate vector */
    inline CoordinateVector() {}
    /** constructor for coordinate vector of a certain dimension */
    inline CoordinateVector( unsigned dimension, long initValue= 0l )
      : vec( dimension, initValue )
    {}
    /** the actual data */
    vector<long> vec;
  };
  
  /** output a coordinate vector to an ostream */
  ostream &operator<<( ostream &os, const CoordinateVector &index );
  
  /** input a coordinate vector from an istream
      (note that the dimension needs to be fixed upfront) */
  istream &operator>>( istream &is, CoordinateVector &index );
  

  
  /** \class CoordinateVectorIter CoordinateVector.hh
      Sequential iterator for indices in an array */
  class CoordinateVectorIter {
  public:
    /** constructor (optionally providing array dimensions) */
    inline CoordinateVectorIter( const CoordinateVector &dim=
				 CoordinateVector() )
    {
      setDimension( dim );
    }
    
    /** constructor from range list specifying a traversal box */
    inline CoordinateVectorIter( const LongRangeList &newBox )
    {
      setBox( newBox );
    }
    
    /** set dimension (create box from 0 to dim-1 on each dimension */
    inline void setDimension( const CoordinateVector &dim )
    {
      box.vec.resize( dim.vec.size() );
      for( int i= 0 ; i< dim.vec.size() ; i++ )
      {
	box.vec[i].val.first= 0l;
	box.vec[i].val.second= dim.vec[i]-1;
      }
      begin();
    }
    
    /** set traversal box */
    inline void setBox( const LongRangeList &newBox )
    {
      box= newBox;
      begin();
    }
    
    /** get traversal box */
    inline const LongRangeList &getBox() const
    {
      return box;
    }
    
    /** move iterator to beginning of array */
    inline void begin()
    {
      int numEntries= box.vec.size();
      pos.vec.clear();
      for( int i= 0 ; i< numEntries ; i++ )
	pos.vec.push_back( box.vec[i].val.first );
      atEnd= false;
    }
    
    /** check if we have reached the end */
    inline bool isAtEnd()
    {
      return atEnd;
    }
    
    /** increment position */
    inline CoordinateVectorIter &operator++()
    {
      return incrComp( 0 );
    }
    
    /** increment position by a whole hyperplane */
    CoordinateVectorIter &incrComp( unsigned int component );
    
    /** get position */
    inline CoordinateVector &getPos()
    {
      return pos;
    }
    
  private:
    
    /** the dimmension of the array */
    LongRangeList box;
    
    /** the current position */
    CoordinateVector pos;
    
    /** for sepeed, store if we have reached the end already */
    bool atEnd;
  };
  
  
  
  
  
  /** \class CoordinateOption CoordinateVector.hh
      parser for sequence of coordinates */
  class CoordinateOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to CoordinateVector object */
    CoordinateOption( CoordinateVector &d,
		      const char *msg= "\tcomma separated coordinate vector\n"
		      "\t(from fastest to slowest changing dimension)\n",
		      const char *longOpt= "--coordinates",
		      const char *shortOpt= "-co" )
      : CommandlineOption( msg, longOpt, shortOpt ), dim( d )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the DataType object */
    CoordinateVector &dim;
  };
  


} /* namespace */



#endif /* BASE_COORDINATEVECTOR_H */

