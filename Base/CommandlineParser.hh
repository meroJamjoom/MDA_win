// ==========================================================================
// $Id: CommandlineParser.hh 277 2009-02-15 01:10:47Z heidrich $
// parser for commandline options
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

#ifndef BASE_COMNMANDLINEPARSER_H
#define BASE_COMNMANDLINEPARSER_H

/*! \file  CommandlineParser.hh
    \brief parser for commandline options
 */

#if defined(_WIN32) || defined(_WIN64)
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>

#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

#include <iostream>
#include <sstream>
#include <list>
#include <limits.h>

#include "Errors.hh"

namespace MDA {

  using namespace std;

  
  
  /** \class CommandlineOption CommandlineParser.hh
      baseclass for parsing a specific commandline option and its arguments
      (the baseclass implements the -h and --help options) */
  class CommandlineOption {
    
  public:
    
    /** constructor: optional paramters are help string and option texts */
    CommandlineOption( const char *msg= "\toutput this help message\n",
		       const char *lTxt= "--help", const char *sTxt= "-h" )
      : helpTxt( msg ), longTxt( lTxt ), shortTxt( sTxt )
    {}
    
    /** determine if this parser is responsible for the provided option */
    virtual bool isResponsible( const char *option );
    
    /** the actual parsing function (pure virtual):
	a return value of false indicates an error */
#if defined (_WIN32) || defined (_WIN64)
	 virtual bool parse( int &index, int argc, char *argv[] );
#else
    virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** the short version of the option string */
    const char *shortTxt;
    
    /** the long version of the option string */
    const char *longTxt;
    
    /** a help text (indented by tab) */
    const char *helpTxt;
  };
  
  
  
  
  /** \class CommandlineParser CommandlineParser.hh
      parser for commandline options */
  class CommandlineParser {

  public:

    /** default constructor */
    CommandlineParser();
    
    /** registering an option */
    void registerOption( CommandlineOption *option )
    {
      options.push_back( option );
    }
    
    /** parse the commandline options */
    bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */ 
    void usage( const char *progName, const char *paramTxt= NULL,
		ostream &os= cerr );
    
    
  protected:
    
    /** a list of all registered options the parser should search for */
    list<CommandlineOption *>	options;

    /** every parser has at least a help option */
    CommandlineOption		helpOption;
  };
  
  
  
  
  /** \class BoolOption CommandlineParser.hh
      parser for a single boolean flag */
  class BoolOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to DataType object */
    BoolOption( bool &val, const char *msg, const char *longTrueOpt,
                const char *shortTrueOpt, const char *longFalseOpt,
		const char *shortFalseOpt )
      : CommandlineOption( msg, longTrueOpt, shortTrueOpt ), value( val ),
	shortFalseTxt( shortFalseOpt ), longFalseTxt( longFalseOpt )
    {}
    
    /** determine if this parser is responsible for the provided option */
    virtual bool isResponsible( const char *option );
    
    /** the actual parsing function */
#if defined (_WIN32) || defined (_WIN64)
	bool parse( int &index, int argc, char *argv[] );
#else
    virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the DataType object */
    bool &value;
    
    /** the short version of the option string for a value of false */
    const char *shortFalseTxt;
    
    /** the long version of the option string for a value of false */
    const char *longFalseTxt;
  };
  
  
  
  
  /** \class ScalarOption CommandlineParser.hh
      parser for a single scalar parameter of any type */
  template <class T> class ScalarOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to scalar object */
    ScalarOption( T &val, const char *msg, const char *longOpt,
                  const char *shortOpt, T _min= INT_MIN, T _max= INT_MAX )
      : CommandlineOption( msg, longOpt, shortOpt ),
#if defined (_WIN32) || defined (_WIN64)
	  value( val ), min_s( _min ), max_s( _max )
#else
	value( val ), min( _min ), max( _max )
#endif
    {}
    
    /** the actual parsing function */
#if defined (_WIN32) || defined (_WIN64)
    virtual bool __forceinline parse( int &index, int argc, char *argv[] );
#else
	virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the DataType object */
    T &value;
    
    /** minimum acceptable value */
#if defined (_WIN32) || defined (_WIN64)
	T min_s;
#else
    T min;
#endif
    /** maximum acceptable value */
#if defined (_WIN32) || defined (_WIN64)
	T max_s;
#else
    T max;
#endif

  };
  
  
  /** parser for a single integer scalar */
  typedef ScalarOption<int> IntOption;
 
  /** parser for a single long integer scalar */
  typedef ScalarOption<long> LongOption;
 
  /** parser for a single floating point scalar */
  typedef ScalarOption<float> FloatOption;
  
  /** parser for a single double precision floating point scalar */
  typedef ScalarOption<double> DoubleOption;
  
  
  
  /** \class FileOption CommandlineParser.hh
      parser for a file name parameter */
  class FileOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to DataType object */
    FileOption( char *&name, const char *msg, const char *longOpt,
		const char *shortOpt )
      : CommandlineOption( msg, longOpt, shortOpt ), fileName( name )
    {}
    
    /** the actual parsing function */
#if defined (_WIN32) || defined (_WIN64)
	 bool parse( int &index, int argc, char *argv[] );
#else
    virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the file name */
    char *&fileName;
  };



  /** \class SelectionOption CommandlineParser.hh
      parser for reading one of a selection of options */
  class SelectionOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to int */
    SelectionOption( int &which, const char *msg, list<const char *>opts )
      : CommandlineOption( msg, NULL, NULL ), options( opts ), selection( which)
    {}
    
    /** constructor with empty option list */
    SelectionOption( int &which, const char *msg )
      : CommandlineOption( msg, NULL, NULL ), selection( which )
    {}
    
    /** add options after object creation */
    inline void addOptions( list<const char *>opts )
    {
      options.splice( options.end(), opts );
    }
    
    /** check if an argument matches one of the possible options */
    virtual bool isResponsible( const char *option );
    
    /** the actual parsing function */
#if defined (_WIN32) || defined (_WIN64)
	 bool parse( int &index, int argc, char *argv[] );
#else
    virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the file name */
    list<const char *> options;
    
    /** which option was selected */
    int &selection;
    
  };



  /** \class StringSelectorOption CommandlineParser.hh
      parser for getting an argument corresponding to a string from a
      given selection */
  class StringSelectorOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to int */
    StringSelectorOption( const char *&str, const char *msg,
			  list<const char *>opts,
			  const char *lTxt, const char *sTxt= NULL )
      : CommandlineOption( msg, lTxt, sTxt), options( opts ), selection( str )
    {}
    
    /** constructor with empty option list */
    StringSelectorOption( const char *&str, const char *msg,
			  const char *lTxt, const char *sTxt= NULL )
      : CommandlineOption( msg, lTxt, sTxt ), selection( str )
    {}
    
    /** add options after object creation */
    inline void addOptions( list<const char *>opts )
    {
      options.splice( options.end(), opts );
    }
    
    /** the actual parsing function */
#if defined (_WIN32) || defined (_WIN64)
	 bool parse( int &index, int argc, char *argv[] );
#else
    virtual bool parse( int &index, int argc, char *argv[] );
#endif
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the file name */
    list<const char *> options;
    
    /** which option was selected */
    const char *&selection;
    
  };



  //
  // ScalarOption template functions
  // 

/** the actual parsing function */

template <class T> bool
 
ScalarOption<T>::parse( int &index, int argc, char *argv[] )

{
  if( index> argc-1 )
    return false;
  istringstream scalar(  argv[index++] );
  scalar >> value;
#if defined (_WIN32) || defined (_WIN64)
  if( !warnCond( value>= min_s && value<= max_s, "  value out of range" ) )
	  return false;
#else
  if( !warnCond( value>= min && value<= max, "  value out of range" ) )
	  return false;
#endif
    
  return true;
}

/** output usage string */
template <class T> void
ScalarOption<T>::usage( ostream &os )
{
  os << "  ";
  if( shortTxt )
    os << shortTxt << " | ";
  os << longTxt << " <num>:\n"
     << helpTxt
     << "\tcurrently: " << value << "\n\n";
}

/*
template<> bool

ScalarOption<char>::parse(int &index, int argc, char *argv[])

{
	if (index > argc - 1)
		return false;

	value = argv[index++][0];
#if defined (_WIN32) || defined (_WIN64)
	if (value< min_s || value> max_s)
#else
	if (value< min || value> max)
#endif
	{
#ifdef DEBUG
		cerr << "Value " << value << " out of range " << min << "..." << max
			<< " for option " << argv[index - 2] << endl;
#endif
		return false;
	}

	return true;

}
*/
} /* namespace */


#endif /* BASE_COMNMANDLINEPARSER_H */

