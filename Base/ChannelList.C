// ==========================================================================
// $Id: ChannelList.C 163 2007-08-31 22:20:52Z heidrich $
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

#ifndef BASE_CHANNELLIST_C
#define BASE_CHANNELLIST_C

#include <ctype.h>

#include "ChannelList.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** write ChannelList to an ostream */
ostream
&operator<<( ostream &os, const ChannelList &channels )
{
  if( !channels.vec.empty() )
  {
    int numElem= channels.vec.size();
    
    os << channels.vec[0];
    for( int i= 1 ; i < numElem ; i++ )
      os << ',' << channels.vec[i];
  }
  
  return os;
}
  

/** read ChannelList from an istream */
istream
&operator>>( istream &is, ChannelList &channels )
{
  unsigned int current, rangeEnd;
  
  channels.vec.clear();
  do
  {
    // read one channel ID, append it to channel list
    is >> current;
    channels.vec.push_back( current );
    
    // if we now read a colon, the user has provided a range...
    if( is.peek()== ':' )
    {
      is.get();
      is >> rangeEnd;
      // error if we didn't read a second number
      if( is.fail() )
	break;
      // append the whole range to the channel list
      for( current++ ; current<= rangeEnd ; current++ )
	channels.vec.push_back( current );
    }
    
    // stop parsing if we no longer get commas
    if( is.eof() )
      break;
    if( is.peek()== ',' )
      is.get();
    else
      break;
  }
  while( !is.fail() );
  
  return is;
}
  


  //
  // ChannelListOption functions
  // 

/** the actual parsing function */
bool
ChannelListOption::parse( int &index, int argc, char *argv[] )
{
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> channels;
  return true; //!!!optStr.fail();
}

/** output usage string */
void
ChannelListOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <index,...>\n"
     << helpTxt
     << "\tCurrently: " << channels << "\n\n";
}



} /* namespace */

#endif /* BASE_CHANNELLIST_C */

