// ==========================================================================
// $Id: ChannelList.hh 413 2009-10-15 10:47:40Z heidrich $
// array of channel numbers
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

#ifndef BASE_CHANNELLIST_H
#define BASE_CHANNELLIST_H

/*! \file  ChannelList.hh
    \brief array of channel numbers
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>
#include <iostream>

#include "CommandlineParser.hh"


namespace MDA {

  using namespace std;
  
  /** \class ChannelList ChannelList.hh
      array of channel numbers */
  class ChannelList{
  public:
    /** constructor with min allocation, but initially empty */
    inline ChannelList( unsigned length= 0 )
    {
      if( length> 0 )
	vec.reserve( length );
    }
    
    /** the actual data */
    vector<unsigned int> vec;
  };
  
  /** write ChannelList to an ostream */
  ostream &operator<<( ostream &os, const ChannelList &channels );
  
  /** read ChannelList from an istream */
  istream &operator>>( istream &is, ChannelList &channels );
  
  
  

  /** \class ChannelListOption ChannelList.hh
      parser for list of channel IDs */
  class ChannelListOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to channel list object */
    ChannelListOption( ChannelList &c,
		       const char *msg= "\tlist of channels numbers\n"
		       "\t(comma separated scalars and/or ranges)\n"
		       "\te.g. \"1,2,3,4,5,9\" is the same as \"1:5,9\"\n",
		       const char *longOpt= "--channellist",
		       const char *shortOpt= "-cl" )
      : CommandlineOption( msg, longOpt, shortOpt ), channels( c )
    {
    }
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the DataType object */
    ChannelList &channels;
  };


  /** a list of coordinate axes */
  typedef ChannelList AxisList;
  
  
} /* namespace */



#endif /* BASE_CHANNELLIST_H */

