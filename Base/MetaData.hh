// ==========================================================================
// $Id:$
// XML metadata parser with data accessors
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

#ifndef METADATA_METADATA_H
#define METADATA_METADATA_H

/*! \file  MetadataReader.hh
    \brief XML metadata parser with data accessors
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>
#include <string>

#include "TinyXML/tinyxml.h"

namespace MDA {

    using namespace std;

  /** \class MetaData MetaData.hh
      XML metadata parser with data accessors */
class MetaData {
  
public:
  
  /** default constructor - creates empty XML document */
  MetaData();
  
  /** construct from existing file */
  inline MetaData( const char *fileName )
  {
    read( fileName );
  }
  
  /** default destructor */
  inline ~MetaData()
  {}
  
  /** read metadata from a separate file */
  bool read( const char *fileName );
  
  /** write to file */
  bool write( const char *fileName );
  
  /** returns the specified node as a string */
  bool get( const char *path, string &val );
  
  /** sets the specified node as a string, creating paths as needed */
  void set( const char *path, const string &val );
  
  
protected:
  
  /** ensure that a certain path exists, by creating any portions
      that may be missing */
  TiXmlElement *make( const char *path );
  
  /** returns the specified node as XML element */
  bool get( const char *path, TiXmlElement **node );
  
  /** chop a path into tokens at a separator character */
  void chopString( string const &path, vector<string> &result,
		   char const sep='.' );
  
  /** extract array index from a path string */
  bool chopArrayString( string const &element, string &node,
			unsigned int &index );
  
  /** extract the attribute portion of a path */
  bool chopAttribString( string const &element, string &node, string &attrib );
  
  /** the XML document */
  TiXmlDocument document;
  
private:
  
  /** copy constructors are not permitted... */
  MetaData( const MetaData& ) {};
  
  /** assignment is not permitted */
  MetaData& operator=( MetaData& m ) const {return m;};
};


//
// get and set functions for arbitrary types
//

/** templated function that returns the specified node casted to whatever
    type is provided.  This function can work on any type compatible with
    the istream operator */
template<class T> bool get( MetaData &m, const char *path, T &val );

/** templated function for reading an inlined array of values */
template<class T> bool get( MetaData &m, const char *path,
			    T *valArray, unsigned count );

/** templated function that sets the specified node casted from
    whatever type is provided.  This function can work on any type
    compatible with the ostream operator */
template<class T> void set( MetaData &m, char const *path, const T &val );

/** template function for writing inlined array of values */
template<class T> void set( MetaData &m, const char *path,
			    const T *valArray, unsigned count );
  

} /* namespace */

#endif /* METADATA_METADATA_H */

