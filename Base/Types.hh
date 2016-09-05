// ==========================================================================
// $Id: Types.hh 199 2008-04-08 03:58:26Z heidrich $
// Data types used in MDA files
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

#ifndef BASE_TYPES_H
#define BASE_TYPES_H

/*! \file  Types.hh
    \brief Data types used in MDA files
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <iostream>

#include "CommandlineParser.hh"


namespace MDA {

  using namespace std;

  /** data types for the multi-dimensional array */
  enum DataType 
  {
    UByte,
    Byte,
    UShort,
    Short,
    UInt,
    Int,
    Float,
    Double,
    UndefinedType // always last!
  };
  
  
  /** size (in bytes) of the individual data types */
  extern int dataTypeSizes[8];
  
  /** write data type name to ostream */
  ostream &operator<<( ostream &os, DataType dataType );
  
  /** read data type name from istream */
  istream &operator>>( istream &is, DataType &dataType );
  
  
  /** a union of all supported types */
  typedef union
  {
    char		byteVal;	/**< unsigned byte */
    unsigned char	ubyteVal;	/**< signed byte */
    short		shortVal;	/**< unsigned short (2 bytes) */
    unsigned short	ushortVal;	/**< signed short (2 bytes) */
    int			intVal;		/**< unsigned int (4 bytes) */
    unsigned int	uintVal;	/**< signed int (4 bytes) */
    float		floatVal;	/**< single prec. float (4 bytes) */
    double		doubleVal;	/**< double prec. float (8 bytes) */
  } AnyType;
  
  /** write one data value of the specified type to an ostream */
  ostream &write( ostream &os, const AnyType &value, const DataType &type );
  
  /** read one data value of the specified type from an istream */
  istream &read( istream &is, AnyType &value, const DataType &type );

  /** convert one or more values from one data type to another */
  void typeConvert( void *fromMem, DataType fromType,
		    void *toMem, DataType toType, unsigned long count= 1 );
  
  
  
  
  /** \class TypeOption Types.hh
      parser for data type */
  class TypeOption: public CommandlineOption {
    
  public:
    
    /** constructor from reference to DataType object */
    TypeOption( DataType &d,
		const char *msg= "\tData type\n",
		const char *lTxt= "--type", const char *sTxt= "-t"  )
      : CommandlineOption( msg, lTxt, sTxt ), dataType( d )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
  
    /** reference to the DataType object */
    DataType &dataType;
    
    
  };
  
  
} /* namespace */



#endif /* BASE_TYPES_H */

