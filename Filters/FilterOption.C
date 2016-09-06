// ==========================================================================
// $Id: FilterOption.C 954 2012-06-05 20:37:36Z krim $
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

#ifndef FILTERS_FILTEROPTION_C
#define FILTERS_FILTEROPTION_C

#include <string.h>

#include "MDA/Base/Errors.hh"
#include "FilterOption.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;
  
const int maxFilterNameLength= 25;
const char *filterNames[] =
{
  "bilateral",
  "bilateralmasked",
  "bilateralgrid",
  "bilateralgridweighted",
  "box",
  "connectedcomponent",
  "corner",
  "dilate",
  "distance",
  "edge",
  "erode",
  "eulernumber",
  "extrema",
  "fastgauss",
  "firstderiv",
  "floodfill",
  "gauss",
  "hat",
  "hitmiss",
  "LoG",
  "median",
  "medianmask",
  "prune",
  "secondderiv",
  "separable",
  "sobel",
  "thin",
  "unsharpmask",
  "thinvoxel",
  "undefined"
};
  

/** write filter name to ostream */
ostream &
operator<<( ostream &os, FilterType filter )
{
  if( filter< UndefinedFiltering )
    os << filterNames[filter] << ' ';
  else
    os << filterNames[UndefinedFiltering] << ' ';
  return os;
}


/** read filter name from istream */
istream &
operator>>( istream &is, FilterType &filter )
{
  int i;
  char buffer[maxFilterNameLength];

  is >> ws;
  for( i= 0 ; i< maxFilterNameLength-1 && isalpha( is.peek() ) ; i++ )
    buffer[i]= is.get();
  buffer[i]= '\0';

  for( int i= 0 ; i< (int)UndefinedFiltering ; i++ )
    if( !strncasecmp( buffer, filterNames[i], maxFilterNameLength ) )
    {
      filter= (FilterType)i;
      return is;
    }
  
  // if we get here, there was trouble parsing the option
  warning( "  could not read filter name\n" );
  is.setstate( ios::failbit );
  return is;
}


//
// FilterOption members
//

/** the actual parsing function */
bool
FilterOption::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> filterType;
  return !optStr.fail();
}

/** output usage string */
void
FilterOption::usage( ostream &os )
{
  int i, length;
  
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <filter name>\n"
     << helpTxt
     << "\tPossible values:";
  for( i= 0 ; i< UndefinedFiltering ; i++ )
  {
    length= strlen( filterNames[i] );
    os << "\n\t  " << filterNames[i++];
    while( i< UndefinedFiltering )
    {
      if( (length+= strlen( filterNames[i] ))> 50 )
	break;
      os << ", " << filterNames[i++];
    }
  }
  os << "\n\tCurrently: " << filterType << "\n\n";
}




} /* namespace */

#endif /* FILTERS_FILTEROPTION_C */

