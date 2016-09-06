// ==========================================================================
// $Id: mda-fromimage.C 260 2008-12-30 09:21:39Z heidrich $
// swizzle the channels in an MDA file
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

#include <assert.h>
#include <Magick++.h>



#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;
using namespace Magick;


int
main( int argc, char *argv[] )
{
  unsigned long i;
  
  CommandlineParser parser;
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "<options> <input file>" );
    exit( 1 );
  }
  
  
  // setup input image
  Image inImage( argv[argc-1] );
  Pixels view( inImage );   

  unsigned int numChannels;
  const char *formatString;
  switch( inImage.type() )
  {
  case BilevelType:
  case GrayscaleType:
    numChannels= 1;
    formatString= "R";
    break;
#if defined (_WIN32) || defined (_WIN64)
  case GrayscaleAlphaType:
	  numChannels= 2;
	  formatString= "RA";
	  break;
#else 
  case GrayscaleMatteType:
    numChannels= 2;
    formatString= "RA";
    break;
#endif
  case TrueColorType:
  case PaletteType:
    numChannels= 3;
    formatString= "RGB";
    break;
 #if defined (_WIN32) || defined (_WIN64)
  case TrueColorAlphaType:
  case PaletteAlphaType:
	  numChannels= 4;
	  formatString= "RGBA";
	  break;
#else 
  case TrueColorMatteType:
  case PaletteMatteType:
    numChannels= 4;
    formatString= "RGBA";
    break;
#endif
  case ColorSeparationType:
    numChannels= 4;
    formatString= "CMYK";
    break;
  default:
    // whatever- shouldn't happen
    numChannels= 3;
    formatString= "RGB";
    break;
  }
  DataType type;
  StorageType magickType;
  switch( inImage.depth() )
  {
  case 8:
    type= UByte;
    magickType= CharPixel;
    break;
  case 16:
    type= UShort;
    magickType= ShortPixel;
    break;
  case 32:
 #if defined (_WIN32) || defined (_WIN64)
	  type= UInt;
	  magickType= LongPixel;
	  break;
#else 
    type= UInt;
    magickType= IntegerPixel;
    break;
#endif
  default:
    // whatever: shouldn't happen
   type= UByte;
    magickType= CharPixel;
    break;
  }    
  
  // set up MDA writer
  CoordinateVector dim;
  dim.vec.push_back( inImage.baseColumns() );
  dim.vec.push_back( inImage.baseRows() );
  
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
 unsigned long numScanlines= writer.getNumScanlinesLeft();
  //short numScanlines = writer.getNumScanlinesLeft();
  assert( numScanlines== dim.vec[1] );
  unsigned long scanlineSizeOut= writer.getScanlineSize();
  //short scanlineSizeOut = writer.getScanlineSize();
 // unsigned char *scanlineOut= new unsigned char[scanlineSizeOut];
  char *scanlineOut = new char[scanlineSizeOut];
  
  for( i= 0 ; i< numScanlines ; i++ )
  {
    // get a scanline from ImageMagick
    inImage.write( 0, i, dim.vec[0], 1,
		   formatString, magickType, scanlineOut );
    
    // and write it to MDA
    writer.writeScanline( scanlineOut );
  }
  if( !writer.disconnect() )
  {
    cerr << "Error writing data\n";
    exit( 1 );
  }
  
  return 0;
}
