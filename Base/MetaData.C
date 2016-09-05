// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: mmt (Matthew Trentacoste)
// Email:   mmt@cs.ubc.ca
// ==========================================================================

#ifndef METADATA_METADATA_C
#define METADATA_METADATA_C

#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "MDA/Config.hh"

#include "Errors.hh"
#include "ChannelList.hh"
#include "CoordinateVector.hh"
#include "Range.hh"

#include "MetaData.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** default constructor */

/** default constructor - creates empty XML document */
MetaData::MetaData()
{
  TiXmlDeclaration *decl= new TiXmlDeclaration( "1.0", "", "" );
  document.LinkEndChild( decl );
}

/** read metadata from a separate file */
bool 
MetaData::read( char const *fileName )
{
  document= TiXmlDocument( fileName );
    
  bool loaded = document.LoadFile();
  warnCond( loaded, "MetaData failed to load file" );
  
  return loaded;
}

/** write to file */
bool
MetaData::write( const char *fileName )
{
  return document.SaveFile( fileName );
}
  

//
// Basic data accessors
//


/** returns the specified node as a string */
bool
MetaData::get( const char *path, string &val )
{
  TiXmlElement *node;
  if( !get( path, &node ) )
    return false;
  
  // Chop off the attribute string
  string mpath;
  string attrib;
  chopAttribString( path, mpath, attrib );
  
  // If path string contains an attribute, grab that, otherwise the value
  if( attrib== "" )
  {
    // Get the child final with the actual data, checking whether it
    // exists or not
    TiXmlNode* tmp= node->FirstChild();
    if( tmp )
#if defined(_WIN32) || defined(_WIN64)
	  val= tmp->Value();
#else
      val= tmp->ValueStr();
#endif
    else
      return false;
    }
  else
  {
    // Get the attribute, checking whether it exists or not
    const char *tmpval = node->Attribute( attrib.c_str() );
    if( !tmpval )
      return false;
    val = tmpval;
  }
  
  return true;
}

/** sets the specified node as a string, creating paths as needed */
void
MetaData::set( const char *path, const string &val )
{
  // get node, creating it if necessary
  TiXmlElement *node= make( path );
  
  // Chop off the attribute string
  string mpath;
  string attrib;
  chopAttribString( path, mpath, attrib );
  
  // If path string contains an attribute, set that, otherwise the value
  if( attrib== "" )
  {
    // Get the child final with the actual data, checking whether it
    // exists. If not, create a new one
    TiXmlNode* tmp = node->FirstChild();
    if( tmp )
    {
#if defined(_WIN32) || defined(_WIN64)
		tmp->SetValue( val.c_str() );
#else
      tmp->SetValue( val );
#endif
    }
    else
    {
#if defined(_WIN32) || defined(_WIN64)
	  TiXmlText *text= new TiXmlText( val.c_str() );
#else
      TiXmlText *text= new TiXmlText( val );
#endif
      node->LinkEndChild( text );
    }
  }
  else
  {
    // Set the attribute
#if defined(_WIN32) || defined(_WIN64)
	  node->SetAttribute( attrib.c_str(), val.c_str() );
#else
    node->SetAttribute( attrib, val );
#endif
  }
}


//
// protected members
//


/** make a new path or parts of a new path, as necessary */
TiXmlElement *
MetaData::make( const char *path )
{
  // Chop off the attribute string
  string mpath;
  string attrib;
  chopAttribString( path, mpath, attrib );
  
  // Divide up into chunks
  vector<string> parts;
  chopString( mpath, parts );
  
  // Get the top-level handle and associated element
  TiXmlHandle prevHandle( &document );
  TiXmlElement *prevElem= NULL;
  
  // loop over elements, getting the next sub-element
  string name;
  unsigned int index;
  for( unsigned int i = 0; i < parts.size() ; i++ )
  {
    chopArrayString( parts[i], name, index );
#if defined(_WIN32) || defined(_WIN64)
    TiXmlHandle thisHandle= prevHandle.Child( name.c_str(), index );
#else
	TiXmlHandle thisHandle= prevHandle.Child( name, index );
#endif
    TiXmlElement *thisElem= thisHandle.ToElement();
    
    // if the element does not yet exit, we create it (repeatedly, if
    // necessary, to get to the right index count
    while( !(thisElem) )
    {
#if defined(_WIN32) || defined(_WIN64)
      thisElem= new TiXmlElement( name.c_str() );
#else
	  thisElem = new TiXmlElement(name);
#endif
      // new root node if necessary...
      if( prevElem )
	prevElem->LinkEndChild( thisElem );
      else
	document.LinkEndChild( thisElem );
      
      // query again, to see if the right count of "name" elements is
      // available now
#if defined(_WIN32) || defined(_WIN64)
	  thisHandle = prevHandle.Child(name.c_str(), index);
#else
      thisHandle= prevHandle.Child( name, index );
#endif
      thisElem= thisHandle.ToElement();
    }
    
    prevHandle= thisHandle;
    prevElem= thisElem;
  }
  
  return prevElem;
}


/** returns the specified node as XML element */
bool
MetaData::get( const char *path, TiXmlElement **node )
{
  // Chop off the attribute string
  string mpath;
  string attrib;
  chopAttribString( path, mpath, attrib );
  
  // Divide up into chunks
  vector<string> parts;
  chopString( mpath, parts );
  
  // Get the top-level handle
  TiXmlHandle handle( &document );
  
  // loop over elements, getting the next sub-element
  for( unsigned int i = 0; i < parts.size(); i++ )
  {
    string name;
    unsigned int index;
    chopArrayString( parts[i], name, index );
#if defined(_WIN32) || defined(_WIN64)
	handle = handle.Child(name.c_str(), index);
#else
    handle= handle.Child( name, index );
#endif
  }
  
  // Convert back to element
  *node = (handle.ToElement());
  
  // Return correct error code
  if( *node )
    return true;
  else
    return false;
}

/** chop a path into tokens at a separator character */
void 
MetaData::chopString( string const &path, vector<string> &result,
		      char const sep )
{
  size_t oldpos = 0;
  size_t newpos = 0;
  
  result.clear();
  
  // Loop over each period found, grab the string to that point, add, move along
  while( newpos != string::npos )
  {
    newpos = path.find( sep, oldpos );
    
    result.push_back( path.substr( oldpos, newpos-oldpos ) );
    oldpos = newpos + 1;
  }
}

/** chop a string into array syntax */
bool
MetaData::chopArrayString( string const &element, string &node,
			   unsigned int &index )
{
  string tnode;
  int tindex = 0;
  
  // Find our braces
  size_t obrace = element.find( '[' );
  size_t cbrace = element.find( ']' );
  
  // If not present, return 0 without error
  if( !obrace && !cbrace )
    return true;
  
  // If one present, return error
  if ( ( obrace && !cbrace ) || ( !obrace && cbrace ) )
    return false;
  
  // Get the value
  tnode   = element.substr( 0, obrace );
  tindex  = strtol( element.substr( obrace+1 ).c_str(), NULL, 0 );
  
  // Some range checking
  if( ( tindex < 0) || ( tindex > 100000 ) )
    return false;
  
  node  = tnode;
  index = tindex;
  return true;
}

/** extract the attribute portion of a path */
bool 
MetaData::chopAttribString( string const &element, string &node,
			    string &attrib )
{
  // Find the colon
  size_t colon = element.find( ':' );
  
  // If no colon/attribute
  if( colon == string::npos )
  {
    node = element;
    attrib = "";
    return true;
  }
  
  // If there is a colon/attribute 
  node =   element.substr( 0, colon );
  attrib = element.substr( colon+1 );
  return true;
}

//
// functions for reading/writing arbitrary data types
//

/** generic version - templated function that returns the specified
    node casted to whatever type is provided. Works with any
    istream-compatible type */
template<class T>
bool 
get( MetaData &m, const char *path, T &val )
{
  // Get the value as a string
  string raw;
  if( !m.get( path, raw ) )
    return false;
  
  // Convert to our datatype
  stringstream ss;
  ss.str( raw );
  ss >> val;
  
  return true;
}


/** generic version - read inlined array of values. Works with any
    istream-compatible type */
template<class T>
bool
get( MetaData &m, const char *path, T *valArray, unsigned count )
{
  // Get the value as a string
  string raw;
  if( !m.get( path, raw ) )
    return false;
  
  // Read values
  stringstream ss;
  ss.str( raw );
  for( unsigned i= 0 ; i< count ; i++ )
    ss >> valArray[i];
  
  return !ss.fail();
}


/** generic version - templated function that sets the specified node
    casted from whatever type is provided.  This function can work on
    any type compatible with the ostream operator */
template<class T>
void
set( MetaData &m, const char *path, const T &val )
{
  // convert value to string
  stringstream oss;
  oss << val;
  
  m.set( path, oss.str() );
}


/** generic version - write inlined array of values. Works with any
    ostream-compatible type */
template<class T>
void 
set( MetaData &m, const char *path, const T *valArray, unsigned count )
{
  // convert value to string
  stringstream oss;
  for( unsigned i= 0 ; i< count ; i++ )
    oss << valArray[i] << ' ';
  
  m.set( path, oss.str() );
}


//
// Explicit instantiations for basic types
//

template bool get( MetaData &, const char *, int & );
template bool get( MetaData &, const char *, unsigned & );
template bool get( MetaData &, const char *, long & );
template bool get( MetaData &, const char *, double & );

template bool get( MetaData &, const char *, IntRange & );
template bool get( MetaData &, const char *, UIntRange & );
template bool get( MetaData &, const char *, LongRange & );
template bool get( MetaData &, const char *, ULongRange & );
template bool get( MetaData &, const char *, FloatRange & );

template bool get( MetaData &, const char *, ChannelList & );
template bool get( MetaData &, const char *, CoordinateVector & );

template bool get( MetaData &, const char *, int *, unsigned );
template bool get( MetaData &, const char *, unsigned *, unsigned );
template bool get( MetaData &, const char *, long *, unsigned );
template bool get( MetaData &, const char *, double *, unsigned );

template void set( MetaData &, const char *, const int & );
template void set( MetaData &, const char *, const unsigned & );
template void set( MetaData &, const char *, const long & );
template void set( MetaData &, const char *, const double & );

template void set( MetaData &, const char *, const IntRange & );
template void set( MetaData &, const char *, const UIntRange & );
template void set( MetaData &, const char *, const LongRange & );
template void set( MetaData &, const char *, const ULongRange & );
template void set( MetaData &, const char *, const FloatRange & );

template void set( MetaData &, const char *, const ChannelList & );
template void set( MetaData &, const char *, const CoordinateVector & );

template void set( MetaData &, const char *, const int *, unsigned );
template void set( MetaData &, const char *, const unsigned *, unsigned );
template void set( MetaData &, const char *, const long *, unsigned );
template void set( MetaData &, const char *, const double *, unsigned );


} /* namespace */

#endif /* METADATA_METADATA_C */
