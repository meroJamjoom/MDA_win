// ==========================================================================
// $Id: CommandlineParser.C 661 2010-03-22 01:19:14Z heidrich $
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

#ifndef BASE_COMMANDLINEPARSER_C
#define BASE_COMMANDLINEPARSER_C

#include <string.h>
#include <sstream>

#include "CommandlineParser.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


  //
  // CommandlineParser functions
  // 

/** default constructor */
CommandlineParser::CommandlineParser()
{
  options.clear();
  options.push_back( &helpOption );
}


/** parse the commandline options */
bool
CommandlineParser::parse( int &index, int argc, char *argv[] )
{
  while( index< argc )
  {
    list<CommandlineOption *>::iterator iter;
    for( iter= options.begin() ; iter!= options.end() ; iter++ )
    {
      if( (*iter)->isResponsible( argv[index] ) )
      {
	if( !(*iter)->parse( ++index, argc, argv ) )
	  // did an error occur?
	  return false;
	break;
      }
    }
    
    // abort if no option feels responsible for the current argument
    // the remaining parameters can be parsed directly by the application
    if( iter == options.end() )
      return true;
  }
  
  return true;
}
    

/** output usage string */ 
void
CommandlineParser::usage( const char *progname, const char *paramTxt,
			  ostream &os )
{
  os << "\nUsage: " << progname << ' '
     << (paramTxt == NULL ? "<options>" : paramTxt) << endl << endl
     << "Valid options:\n";
  for( list<CommandlineOption *>::iterator iter= options.begin() ;
       iter!= options.end() ; iter++ )
    (*iter)->usage( os );
}



  //
  // CommandlineOption functions
  // 


/** determine if this parser is responsible for the provided option */
bool
CommandlineOption::isResponsible( const char *option )
{
  return (shortTxt!= NULL && !strcmp( option, shortTxt ))
    || !strcmp( option, longTxt );
}


/** the actual parsing function (pure virtual):
    a return value of false indicates an error */
bool
CommandlineOption::parse( int &index, int argc, char *argv[] )
{
  // trigger the output of the help text
  return false;
}
  
/** output usage string */
void
CommandlineOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << endl
     << helpTxt << endl;
}




  //
  // BoolOption functions
  // 

/** determine if this parser is responsible for the provided option */
bool
BoolOption::isResponsible( const char *option )
{
  return
    (shortTxt!= NULL && !strcmp( option, shortTxt ))
    || !strcmp( option, longTxt )
    || (shortFalseTxt!= NULL && !strcmp( option, shortFalseTxt ))
    || !strcmp( option, longFalseTxt );
}


/** the actual parsing function */
bool
BoolOption::parse( int &index, int argc, char *argv[] )
{
  if( (shortTxt!= NULL && !strcmp( argv[index-1], shortTxt ))
      || !strcmp( argv[index-1], longTxt ) )
    value= true;
  else
    value= false;
  
  return true;
}

/** output usage string */
void
BoolOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt )
    os << '(' << shortTxt << " | ";
  os << longTxt;
  if( shortTxt )
    os << ')';
  os << " | ";
  if( shortFalseTxt )
    os << '(' << shortFalseTxt << " | ";
  os << longFalseTxt;
  if( shortFalseTxt )
    os << ')';
  os << endl
     << helpTxt
     << "\tCurrently: " << (value ? longTxt : longFalseTxt) << "\n\n";
}




  //
  // FileOption functions
  // 

/** the actual parsing function */
bool
FileOption::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  fileName= argv[index++];
  return true;
}

/** output usage string */
void
FileOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt )
    os << shortTxt << " | ";
  os << longTxt << " <file name>:\n"
     << helpTxt
     << "\tCurrently: " << (fileName== NULL ? "(null)" : fileName) << "\n\n";
}


  //
  // SelectionOption functions
  // 

/** check if an argument matches one of the possible options */
bool
SelectionOption::isResponsible( const char *option )
{
  list<const char *>::iterator iter;
  for( iter= options.begin() ; iter!= options.end() ; iter++ )
    if( !strcmp( *iter, option ) )
      return true;
  return false;
}
    
/** the actual parsing function */
bool
SelectionOption::parse( int &index, int argc, char *argv[] )
{
  list<const char *>::iterator iter;
  int i;
  for( i= 0, iter= options.begin() ; iter!= options.end() ; iter++, i++ )
    if( !strcmp( *iter, argv[index-1] ) )
    {
      selection= i;
      return true;
    }
  
  i= -1;
  return false;
}

/** output usage string */
void
SelectionOption::usage( ostream &os )
{
  const char *which= NULL;
  int i;
  list<const char *>::iterator iter;
  for( i= 0, iter= options.begin() ; iter!= options.end() ; iter++, i++ )
  {
    os << "  " << *iter << endl;
    if( i== selection )
      which= *iter;
  }
  os << helpTxt
     << "\tCurrently: " << (which== NULL ? "(undefined)" : which) << "\n\n";
}



  //
  // StringSelectorOption functions
  // 

/** the actual parsing function */
bool
StringSelectorOption::parse( int &index, int argc, char *argv[] )
{
  list<const char *>::iterator iter;
  
  if( index< argc )
    for( iter= options.begin() ; iter!= options.end() ; iter++ )
      if( !strcasecmp( *iter, argv[index] ) )
      {
	selection= *iter;
	index++;
	return true;
      }
  
  return false;
}

/** output usage string */
void
StringSelectorOption::usage( ostream &os )
{
  list<const char *>::iterator iter;

  os << "  " << longTxt;
  if( shortTxt!= NULL )
    os << " | " << shortTxt;
  os << " <flag>\n\twith <flag> being one of:\n";
  for( iter= options.begin() ; iter!= options.end() ; iter++ )
    os << "\t  " << *iter << endl;
  os << helpTxt
     << "\tCurrently: " << (selection== NULL ? "(undefined)" : selection)
     << "\n\n";
}






/** the parsing function, specialized for char */
template<> bool
ScalarOption<char>::parse( int &index, int argc, char *argv[] )
{
  if( index> argc-1 )
    return false;
  value= argv[index++][0];

#if defined (_WIN32) || defined (_WIN64)
  if( !warnCond( value>= min_s && value<= max_s, "  value out of range" ) )
	  return false;
#else
  if( !warnCond( value>= min && value<= max, "  value out of range" ) )
	  return false;
#endif
   
  return true;
}


template class ScalarOption<char>;


} /* namespace */

#endif /* BASE_COMMANDLINEPARSER_C */

