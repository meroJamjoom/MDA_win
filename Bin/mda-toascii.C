// ==========================================================================
// $Id: mda-toascii.C 319 2009-05-27 21:17:01Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: atcheson (Bradley Atcheson)
// Email:   atcheson@cs.ubc.ca
// ==========================================================================

#include <string.h>

#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "\n\n\
Echo the pixel values (in scanline order) to the console in text form\n"


int main( int argc, char** argv )
{
  CommandlineParser parser;

  // cmdline args
  //...
  
  bool coordOutput= false;
  BoolOption coordOpt( coordOutput,
		       "\twhether or not to print the pixel coordinates\n",
		       "--coord", NULL, "--no-coord", NULL );
  parser.registerOption( &coordOpt );
  
  bool collateChannels= false;
  BoolOption collateOpt( collateChannels,
		 "\twhether or not to collate all channels onto one line\n",
			 "--collate", NULL, "--no-collate", NULL );
  parser.registerOption( &collateOpt );
  
  
  
  // parse options
  int index = 1;
  if( !parser.parse(index,argc,argv) ) {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }

  // setup MDA reader and read input header
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() ) {
    cerr << "Cannot read file header" << endl;
    exit( 1 );
  }
  DataType type = reader.getType();
  CoordinateVector dim = reader.getDim();
  int numChannels = reader.getNumChannels();
  unsigned long numScanlines = reader.getNumScanlinesLeft();

  // keep track of the position in the array
  CoordinateVectorIter iter( dim );
  iter.begin();
  
  // read each scanline out to stdout
  unsigned int bytesPerValue = dataTypeSizes[type];
#if defined (_WIN32) || defined (_WIN64)
  char* buffer=new char[bytesPerValue];
#else
  char buffer[bytesPerValue];
#endif
  for( unsigned long i= 0; i< numScanlines; i++ )
  {
    char* scanline = (char*)reader.readScanline();
    // for each pixel
    for( unsigned long p= 0 ; p<dim.vec[0]; p++ ) {
      
      // in collation mode, output one set of coordinates per pixel
      if( collateChannels && coordOutput )
	  for( unsigned h= 0 ; h< dim.vec.size() ; h++ )
	    printf( "%d\t", iter.getPos().vec[h] );
      
      // for each channel
      for( unsigned long c=0; c<numChannels; c++ ) {
	
	// in un-collated mode, output coordinates once per channel
	if( !collateChannels && coordOutput )
	  for( unsigned h= 0 ; h< dim.vec.size() ; h++ )
	    printf( "%d\t", iter.getPos().vec[h] );
	
        // copy just one data value to a temp buffer
        memcpy( buffer,
                scanline + (p*numChannels + c)*bytesPerValue,
                bytesPerValue );
        // output to stdout
        switch( type ) {
          case UByte: {
            unsigned char val = *((unsigned char*)buffer);
            printf( "%d", val );
            break;
          }
          case Byte: {
            char val = *((char*)buffer);
            printf( "%d", val );
            break;
          }
          case UShort: {
            unsigned short val = *((unsigned short*)buffer);
            printf( "%d", val );
            break;
          }
          case Short: {
            short val = *((short*)buffer);
            printf( "%d", val );
            break;
          }
          case UInt: {
            unsigned int val = *((unsigned int*)buffer);
            printf( "%d", val );
            break;
          }
          case Int: {
            int val = *((int*)buffer);
            printf( "%d", val );
            break;
          }
          case Float: {
            float val = *((float*)buffer);
            printf( "%f", val );
            break;
          }
          case Double: {
            double val = *((double*)buffer);
            printf( "%f", val );
            break;
          }
#if defined (_WIN32) || defined (_WIN64)
			delete[]buffer;
#endif

        }
	
	// if there are more channels to print on this line, output a
	// tab, otherwise a newline
	if( !collateChannels || c== numChannels-1 )
	  printf( "\n" );
	else
	  printf( "\t" );
      }
      
      // advance pixel position
      ++iter;
    }
  }

  reader.disconnect();
  return 0;
}

