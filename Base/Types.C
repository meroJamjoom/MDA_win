// ==========================================================================
// $Id: Types.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef BASE_TYPES_C
#define BASE_TYPES_C

#include <string.h>

#include "Types.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inlcusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** size (in bytes) of the individual data types */
int dataTypeSizes[8]=
{
  1, 1, 2, 2, 4, 4, 4, 8
};

static const int maxDataTypeNameLength= 7;
static const char dataTypeNames[9][maxDataTypeNameLength]=
{
  "ubyte", "byte", "ushort", "short", "uint", "int", "float", "double",
  "(nil)"
};


/** write one data value of the specified type to an ostream */
ostream &
write( ostream &os, const AnyType &value, const DataType &type )
{
  short val;
  switch( type )
  {
  case Byte:
    val= value.byteVal;
    os << val;
    break;
  case UByte:
    val= value.ubyteVal;
    os << val;
    break;
  case Short:
    os << value.shortVal;
    break;
  case UShort:
    os << value.ushortVal;
    break;
  case Int:
    os << value.intVal;
    break;
  case UInt:
    os << value.uintVal;
    break;
  case Float:
    os << value.floatVal;
    break;
  case Double:
    os << value.doubleVal;
    break;
  default:
    // nothing to do...
    break;
  }
  return os;
}
  
/** read one data value of the specified type from an istream */
istream &
read( istream &is, AnyType &value, const DataType &type )
{
  short val;
  switch( type )
  {
  case Byte:
    is >> val;
    value.byteVal= val;
    break;
  case UByte:
    is >> val;
    value.ubyteVal= val;
    break;
  case Short:
    is >> value.shortVal;
    break;
  case UShort:
    is >> value.ushortVal;
    break;
  case Int:
    is >> value.intVal;
    break;
  case UInt:
    is >> value.uintVal;
    break;
  case Float:
    is >> value.floatVal;
    break;
  case Double:
    is >> value.doubleVal;
    break;
  default:
    // nothing to do...
    break;
  }
  return is;
}


/** write data type name to ostream */
ostream &
operator<<( ostream &os, DataType dataType )
{
  return os << ((dataType>= UByte && dataType<= Double) ?
		dataTypeNames[dataType] : dataTypeNames[UndefinedType]);
}


/** read data type name from istream */
istream &
operator>>( istream &is, DataType &dataType )
{
  int i;
  char buffer[maxDataTypeNameLength];

  is >> ws;
  for( i= 0 ; i< maxDataTypeNameLength-1 && isalpha( is.peek() ) ; i++ )
    buffer[i]= is.get();
  buffer[i]= '\0';
  
  for( int i= 0 ; i< UndefinedType ; i++ )
    if( !strncasecmp( buffer, dataTypeNames[i], maxDataTypeNameLength ) )
    {
      dataType= (DataType)i;
      return is;
    }
  
  // if we get here, something went wrong
  warning( "  could not read data type!" );
  is.setstate( ios::failbit );
  return is;
}


/** convert one or more values from one data type to another */
void
typeConvert( void *fromMem, DataType fromType,
	     void *toMem, DataType toType, unsigned long count )
{
  float		f;
  double	d;
  long		i;
  
  // if to and from type are the same, use memcpy
  if( fromType== toType )
  {
    memcpy( toMem, fromMem, count*dataTypeSizes[fromType] );
    return;
  } 
  else switch( toType )
  {
  case UByte:
    // convert to unsigned byte
    {
      unsigned char *dest= (unsigned char *)toMem;
      switch( fromType )
      {
      case Byte:
	// from signed byte to unsigned byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((char *)fromMem)[i] + 128;
	break;
      case UShort:
	// from unsigned short to unsigned byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned short *)fromMem)[i] >> 8;
	break;
      case Short:
	// from signed short to unsigned byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((short *)fromMem)[i] >> 8) + 128;
	break;
      case UInt:
	// from unsigned int to unsigned byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned int *)fromMem)[i] >> 24;
	break;
      case Int:
	// from signed int to unsigned byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]=  (((int *)fromMem)[i] >> 24) + 128;
	break;
      case Float:
	// from float to unsigned byte
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 255.0 + 0.5;
	  dest[i]= f< 0.0 ? 0 : (f> 255.0 ? 255 : (unsigned char)f);
	}
	break;
      case Double:
	// from double to unsigned byte
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 255.0 + 0.5;
	  dest[i]= d< 0.0 ? 0 : (d> 255.0 ? 255 : (unsigned char)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case Byte:
    // convert to signed byte
    {
      char *dest= (char *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to signed byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned char *)fromMem)[i] - 128;
	break;
      case UShort:
	// from unsigned short to signed byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned short *)fromMem)[i] >> 8) - 128;
	break;
      case Short:
	// from signed short to signed byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((short *)fromMem)[i] >> 8;
	break;
      case UInt:
	// from unsigned int to signed byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned int *)fromMem)[i] >> 24) - 128;
	break;
      case Int:
	// from signed int to signed byte
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((int *)fromMem)[i] >> 24;
	break;
      case Float:
	// from float to signed byte
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 128.0 + 0.5;
	  dest[i]= f< -128.0 ? -128 : (f> 127.0 ? 127 : (char)f);
	}
	break;
      case Double:
	// from double to signed byte
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 128.0 + 0.5;
	  dest[i]= d< -128.0 ? -128 : (d> 127.0 ? 127 : (char)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case UShort:
    // convert to unsigned short
    {
      unsigned short *dest= (unsigned short *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to unsigned short
	for( i= 0 ; i < count ; i++ )
	{
	  unsigned char *ucp= (unsigned char *)(dest+i);
	  ucp[0]= ucp[1]= ((unsigned char *)fromMem)[i];
	}
	break;
      case Byte:
	// from signed byte to unsigned short
	for( i= 0 ; i < count ; i++ )
	{
	  char *cp= (char *)(dest+i);
	  cp[0]= cp[1]= ((char *)fromMem)[i] + 128;
	}
	break;
      case Short:
	// from signed short to unsigned short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((short *)fromMem)[i] + 32768;
	break;
      case UInt:
	// from unsigned int to unsigned short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned int *)fromMem)[i] >> 16;
	break;
      case Int:
	// from signed int to unsigned short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((int *)fromMem)[i] >> 16) + 32768;
	break;
      case Float:
	// from float to unsigned short
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 65535.0 + 0.5;
	  dest[i]= f< 0.0 ? 0 : (f> 65535.0 ? 65535 : (unsigned short)f);
	}
	break;
      case Double:
	// from double to unsigned short
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 65535.0 + 0.5;
	  dest[i]= d< 0.0 ? 0 : (d> 65535.0 ? 65535 : (unsigned short)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case Short:
    // convert to signed short
    {
      short *dest= (short *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to signed short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned char *)fromMem)[i] << 8) - 32768;
	break;
      case Byte:
	// from signed byte to signed short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((char *)fromMem)[i] << 8);
	break;
      case UShort:
	// from unsigned short to signed short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned short *)fromMem)[i] - 32768;
	break;
      case UInt:
	// from unsigned int to signed short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned int *)fromMem)[i] >> 16) - 32768;
	break;
      case Int:
	// from signed int to signed short
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((int *)fromMem)[i] >> 16;
	break;
      case Float:
	// from float to signed short
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 32768.0 + 0.5;
	  dest[i]=
	    f< -32768.0 ? -32768 : (f> 32767.0 ? 32767 : (short)f);
	}
	break;
      case Double:
	// from double to signed short
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 32768.0 + 0.5;
	  dest[i]=
	    d< -32768.0 ? -32768 : (d> 32767.0 ? 32767 : (short)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case UInt:
    // convert to unsigned int
    {
      unsigned int *dest= (unsigned int *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  unsigned char *ucp= (unsigned char *)(dest+i);
	  ucp[0]= ucp[1]= ucp[2]= ucp[3]= ((unsigned char *)fromMem)[i];
	}
	break;
      case Byte:
	// from signed byte to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  char *cp= (char *)(dest+i);
	  cp[0]= cp[1]= cp[2]= cp[3]= ((char *)fromMem)[i] + 128;
	}
	break;
      case UShort:
	// from unsigned short to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  unsigned short *usp= (unsigned short *)(dest+i);
	  usp[0]= usp[1]= ((unsigned short *)fromMem)[i];
	}
	break;
      case Short:
	// from signed short to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  short *sp= (short *)(dest+i);
	  sp[0]= sp[1]= ((short *)fromMem)[i]+ 32768;
	}
	break;
      case Int:
	// from signed int to unsigned int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((int *)fromMem)[i] + 2147483648U;
	break;
      case Float:
	// from float to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 4294967295.0 + 0.5;
	  dest[i]=
	    f< 0.0 ? 0 : (f> 4294967295.0 ? 4294967295U : (unsigned int)f);
	}
	break;
      case Double:
	// from double to unsigned int
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 4294967295.0 + 0.5;
	  dest[i]=
	    d< 0.0 ? 0 : (d> 4294967295.0 ? 4294967295U : (unsigned int)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case Int:
    // convert to signed int
    {
      int *dest= (int *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to signed int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned char *)fromMem)[i] << 24) - 2147483648LL;
	break;
      case Byte:
	// from signed byte to signed int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((char *)fromMem)[i] << 24;
	break;
      case UShort:
	// from unsigned short to signed int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= (((unsigned short *)fromMem)[i] << 16) - 2147483648LL;
	break;
      case Short:
	// from signed short to signed int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((short *)fromMem)[i] << 16;
	break;
      case UInt:
	// from unsigned int to signed int
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned int *)fromMem)[i] - 2147483648LL;
	break;
      case Float:
	// from float to signed int
	for( i= 0 ; i < count ; i++ )
	{
	  f= ((float *)fromMem)[i] * 2147483648.0 + 0.5;
	  dest[i]=
	    f< -2147483648.0 ? -2147483648LL :
	    (f> 2147483647.0 ? 2147483647 : (int)f);
	}
	break;
      case Double:
	// from double to signed int
	for( i= 0 ; i < count ; i++ )
	{
	  d= ((double *)fromMem)[i] * 2147483648.0 + 0.5;
	  dest[i]=
	    d< -2147483648.0 ? -2147483648LL :
	    (d> 2147483647.0 ? 2147483647 : (int)d);
	}
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case Float:
    // convert to float
    {
      float *dest= (float *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned char *)fromMem)[i] / 255.0;
	break;
      case Byte:
	// from signed byte to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((char *)fromMem)[i] / 128.0;
	break;
      case UShort:
	// from unsigned short to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned short *)fromMem)[i] / 65535.0;
	break;
      case Short:
	// from signed short to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((short *)fromMem)[i] / 32768.0;
	break;
      case UInt:
	// from unsigned int to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned int *)fromMem)[i] / 4294967295.0;
	break;
      case Int:
	// from signed int to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((int *)fromMem)[i] / 2147483648.0;
	break;
      case Double:
	// from double to float
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((double *)fromMem)[i];
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
    
  case Double:
    // convert to double
    {
      double *dest= (double *)toMem;
      switch( fromType )
      {
      case UByte:
	// from unsigned byte to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned char *)fromMem)[i] / 255.0;
	break;
      case Byte:
	// from signed byte to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((char *)fromMem)[i] / 128.0;
	break;
      case UShort:
	// from unsigned short to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned short *)fromMem)[i] / 65535.0;
	break;
      case Short:
	// from signed short to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((short *)fromMem)[i] / 32768.0;
	break;
      case UInt:
	// from unsigned int to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((unsigned int *)fromMem)[i] / 4294967295.0;
	break;
      case Int:
	// from signed int to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((int *)fromMem)[i] / 2147483648.0;
	break;
      case Float:
	// from float to double
	for( i= 0 ; i < count ; i++ )
	  dest[i]= ((float *)fromMem)[i];
	break;
      default:
	// src = dest, so should never happen (caught above)
	break;
      }
      break;
    }
  default:
    // undefined type; nothing to do
    break;
  }
}



  //
  // TypeOption functions
  // 

/** the actual parsing function */
bool
TypeOption::parse( int &index, int argc, char *argv[] )
{
  if( index>  argc-1 )
    return false;
  string param= argv[index++];
  istringstream optStr( param );
  optStr >> dataType;
  return !optStr.fail();
}

/** output usage string */
void
TypeOption::usage( ostream &os )
{
  int i, length;
  
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <type>\n"
     << helpTxt
     << "\tPossible values:";
  for( i= 0 ; i<= Double ; i++ )
  {
    length= strlen( dataTypeNames[i] );
    os << "\n\t  " << dataTypeNames[i++];
    while( i<= Double )
    {
      if( (length+= strlen( dataTypeNames[i] ))> 50 )
	break;
      os << ", " << dataTypeNames[i++];
    }
  }
  os << "\n\tCurrently: " << dataType << "\n\n";
}


} /* namespace */

#endif /* BASE_TYPES_C */

