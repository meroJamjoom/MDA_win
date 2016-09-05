// ==========================================================================
// $Id: CoordinateVector.C 676 2010-04-04 00:03:47Z heidrich $
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

#ifndef ARRAY_ARRAYINDEX_C
#define ARRAY_ARRAYINDEX_C

#include <string.h>
#include <iostream>

#include "Errors.hh"
#include "CoordinateVector.hh"


namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inlcusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;



/** output a coordinate vector to an ostream */
ostream &
operator<<( ostream &os, const CoordinateVector &index )
{
  int numElem= index.vec.size();
  if( numElem> 0 )
  {
    os << index.vec[0];
    for( int i= 1 ; i< numElem ; i++ )
      os << ',' << index.vec[i];
    os << ' ';
  }
  return os;
}

 
/** input a coordinate vector from an istream */
istream &
operator>>( istream &is, CoordinateVector &index )
{
  unsigned long coord;
  
  index.vec.clear();
  do
  {
    // read one coordinate and append it to the vector
    is >> coord;
    if( is.fail() )
      break;
    index.vec.push_back( coord );
    
    // stop parsing if we no longer get commas
    if( is.eof() )
      break;
    if( is.peek()== ',' )
      is.get();
    else
      break;
  }
  while( !is.fail() );
  
  warnCond( !is.fail(), "  could not read CoordinateVector!" );
  return is;
}




  //
  // CoordinateVectorIter functions
  // 


/** increment position */
CoordinateVectorIter &
CoordinateVectorIter::incrComp( unsigned int component )
{
  if( atEnd )
    return *this;
  
  int numElem= box.vec.size();
  for( int i= component ; i< numElem ; i++ )
  {
    pos.vec[i]++;
    if( pos.vec[i]<= box.vec[i].val.second )
      return *this;
    else
      pos.vec[i]= box.vec[i].val.first;
  }
  
  // if we reach this, we are through with the array
  atEnd= true;
  return *this;
}



  //
  // CoordinateOption functions
  // 

/** the actual parsing function */
bool
CoordinateOption::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> dim;
  return !optStr.fail();
}
    
/** output usage string */
void
CoordinateOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <value,value,...>\n"
     << helpTxt
     << "\tCurrently: " << dim << "\n\n";
}


} /* namespace */

#endif /* ARRAY_ARRAYINDEX_C */

