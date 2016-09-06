// ==========================================================================
// $Id: ResamplerOption.C 319 2009-05-27 21:17:01Z heidrich $
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

#ifndef RESAMPLING_RESAMPLEOPTION_C
#define RESAMPLING_RESAMPLEOPTION_C

#include <string.h>
#include "MDA/Base/Errors.hh"
#include "Resampling.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


const int maxResamplerNameLength= 15;
const char *resamplerNames[] =
{
  "nearest",
  "box",
  "linear",
  "cubic",
  "gauss",
  "min",
  "max",
  "undefined"
};


/** write resampler name to ostream */
ostream &
operator<<( ostream &os, ResampleMode resampler )
{
  if( resampler< UndefinedResampling )
    os << resamplerNames[resampler] << ' ';
  else
    os << resamplerNames[UndefinedResampling] << ' ';
  return os;
}


/** read resampler name from istream */
istream &
operator>>( istream &is, ResampleMode &resampler )
{
  int i;
  char buffer[maxResamplerNameLength];

  is >> ws;
  for( i= 0 ; i< maxResamplerNameLength-1 && isalpha( is.peek() ) ; i++ )
    buffer[i]= is.get();
  buffer[i]= '\0';
  
  for( int i= 0 ; i< (int)UndefinedResampling ; i++ )
    if( !strncasecmp( buffer, resamplerNames[i], maxResamplerNameLength ) )
    {
      resampler= (ResampleMode)i;
      return is;
    }
  
  // if we get here, things have gone wrong...
  warning( "  could not read resampler name!" );
  is.setstate( ios::failbit );
  return is;
}


//
// ResamplerOption members
//

/** the actual parsing function */
bool
ResamplerOption::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> resampleMode;
  return !optStr.fail();
}


/** output usage string */
void
ResamplerOption::usage( ostream &os )
{
  int i, length;
  
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <resampling filter>\n"
     << helpTxt
     << "\tPossible values:";
  for( i= 0 ; i< UndefinedResampling ; i++ )
  {
    length= strlen( resamplerNames[i] );
    os << "\n\t  " << resamplerNames[i++];
    while( i< UndefinedResampling )
    {
      if( (length+= strlen( resamplerNames[i] ))> 50 )
	break;
      os << ", " << resamplerNames[i++];
    }
  }
  os << "\n\tCurrently: " << resampleMode << "\n\n";
}

} /* namespace */

#endif /* RESAMPLING_RESAMPLEOPTION_C */

