#include <iostream>
#include <assert.h>

#include "MDA/Base/Errors.hh"
#include "MDA/Base/ChannelList.hh"
#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Array/MDAFileIO.hh"
#include "MDA/Color/ColorSpaceFactory.hh"

using namespace std;
using namespace MDA;

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // setup options

  // channel list
  ChannelList channels;
  ChannelListOption channelOpt( channels );
  parser.registerOption( &channelOpt );
  
  // source color space
  ColorSpaceFactory colorFac;
  colorFac.registerOptions( parser );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], "<options>\n"
		  "Convert from XYZ into a given color space\n" );
    exit( 1 );
  }
  
  // create the color space
  ColorSpace *space= colorFac.makeColorSpace();
  
  // setup MDA reader and read input header
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() )
  {
    cerr << "Cannot read file header\n";
    exit( 1 );
  }
  DataType      type= reader.getType();
  CoordinateVector      dim= reader.getDim();
  unsigned int  numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  char *scanline;
  
  // check channels against dimensionality of space
  if( channels.vec.size()== 0 )
    channels= reader.allChannels();
  errorCond( channels.vec.size()== 3,
	     "  number of channels does not match dimension of XYZ space!" );
  unsigned int numOutChannels= space->getDimension();
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, channels.vec.size(), type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  assert( scanlineSize== writer.getScanlineSize() );
  
  // actually convert the data
  Vector XYZ( 3 );
  Vector dstSpace( numChannels );
  
  unsigned int bytesPerValue= dataTypeSizes[type];
  char *outScanline= new char[3*dim.vec[0]*bytesPerValue];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    // read scanline
    scanline= (char *)reader.readScanline();
    
    // for each pixel, gamma correct each selected channel
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      for( k= 0 ; k< channels.vec.size() ; k++ )
      {
        double value;
        
        // convert current channel value to double
        typeConvert( scanline+ (j*numChannels + channels.vec[k])*bytesPerValue,
                     type, &value, Double, 1 );
	
	XYZ[k]= value;
      }
      
      // convert one pixel
      space->fromXYZ( XYZ, dstSpace );
      
      // write to output scanline
      for( k= 0 ; k< numOutChannels ; k++ )
	typeConvert( &dstSpace[k], Double, 
		     outScanline+ (j*3+k)*bytesPerValue, type, 1 );
    }
    
    // write scanline
    writer.writeScanline( outScanline );
  }
  
  delete [] outScanline;
  return 0;
}
