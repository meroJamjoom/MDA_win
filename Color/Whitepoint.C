// ==========================================================================
// $Id:$
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

#ifndef COLOR_WHITEPOINT_C
#define COLOR_WHITEPOINT_C

#include "Whitepoint.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

//
// declarations and data local to this file
//

class NamedWhitepoint {
public:
  const char *name;
  double x;
  double y;
};
    

NamedWhitepoint namedWhitepoints[]= {
  { "D65",	0.3127,		0.3290 }, // first entry is the default
  { "A",	0.45,		0.41 },
  { "B",	0.35,		0.35 },
  { "C",	0.3101,		0.3161 },
  { "D50",	0.3457,		0.3584 },
  { "D55",	0.33,		0.35 },
  { "D75",	0.30,		0.32 },
  { "E",	0.3333,		0.3333 },
  { "" } // end of list
};


//
// whitepoint member data
//

/** set whitepoint acording to a name in the list */
void
Whitepoint::setWhitepoint( const char *name, double relativeScale )
{
  unsigned i= 0;
  // find whitepoint standard
  while( namedWhitepoints[i].name[0]!= '\0' )
    if( !strcmp( name, namedWhitepoints[i].name ) )
      break;
    else
      i++;
  if( namedWhitepoints[i].name[0]== '\0' )
  {
    i= 0;
    cerr << "Whitepoint \"" << name << "\" unknown - using \""
	 << namedWhitepoints[i].name  << '\"' << endl;
  }
  
  // compute white point in XYZ
  XYZ[1]= relativeScale;
  XYZ[0]= relativeScale * namedWhitepoints[i].x / namedWhitepoints[i].y;
  XYZ[2]= relativeScale * (1.0-namedWhitepoints[i].x-namedWhitepoints[i].y) /
    namedWhitepoints[i].y;
}


/** return a list of the names of known standard gammas
    (names are appended to end of existing list) */
void
Whitepoint::getStandardNames( list<const char *> &names )
{
  unsigned i= 0;
  while( namedWhitepoints[i].name[0]!= '\0' )
    names.push_back( namedWhitepoints[i++].name );
}

//
// I/O methods
//

/** write Whitepoint to an ostream */
ostream &
operator<<( ostream &os, const Whitepoint &wp )
{
  Vector XYZ( wp.getXYZ() );
  double sum= XYZ[0]+XYZ[1]+XYZ[2];
  os << XYZ[0]/sum << ',' << XYZ[1]/sum;
  return os;
}

/** read Whitepoint from an istream */
istream &
operator>>( istream &is, Whitepoint &wp )
{
  double x, y;
  char std[16];
  
  is.clear();
  is >> ws >> x;
  if( is.fail() )
  {
    is.clear();
    is >> ws >> std;
    wp.setWhitepoint( std );
  }
  else
    if( is.peek()== ',' )
    {
      is.get();
      is >> y;
      wp.setWhitepoint( x, y );
    }
    else
      is.setstate( ios::failbit );
  return is;
}


//
// whitepoint parser
//

/** the actual parsing function */
bool
WhitepointOption::parse( int &index, int argc, char *argv[] )
{
  if( index>  argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> wp;
  return !optStr.fail();
}


/** output usage string */
void
WhitepointOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " (<x>,<y>)|<standard>\n"
     << helpTxt
     << "\tRecognized standards:\n\t ";
  
  list<const char *> standards;
  Whitepoint::getStandardNames( standards );
  for( list<const char *>::iterator iter= standards.begin() ;
       iter!= standards.end() ;
       iter++ )
    os << ' ' << *iter;
  os << "\n\tCurrently: " << wp << "\n\n";
}


} /* namespace */

#endif /* COLOR_WHITEPOINT_C */

