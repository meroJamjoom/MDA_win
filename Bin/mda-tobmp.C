// ==========================================================================
// $Id: mda-tobmp.C 199 2008-04-08 03:58:26Z heidrich $
// convert MDA stream to image using the ImageMagick library
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

#include <stdio.h>

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


// BMP header datastructure
typedef struct
{
  short		pad;	// need this to avoid byte alignmnet issues
  short		bfType;	// the real header starts here
  unsigned int	bfSize;
  unsigned int	bfReserved;
  unsigned int	bfOffBits;
  
  unsigned int	biSize;
  unsigned int	biWidth;
  unsigned int	biHeight;
  short		biPlanes;
  short		biBitCount;
  unsigned int	biCompression;
  unsigned int	biSizeImage;
  unsigned int	biXPelsPerMeter;
  unsigned int	biYPelsPerMeter;
  unsigned int	biClrUsed;
  unsigned int	biClrImportant;
} BMP_Header;



int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  BMP_Header imageHeader= {
    0, 19778, 0, 0, 54,
    40, 0, 0, 1, 24, 0, 0, 0, 0, 0, 0 
  };
			 
  CommandlineParser parser;
  
  // setup options
  ChannelList	channels;
  ChannelListOption channelOption( channels );
  parser.registerOption( &channelOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "<options> <bmp outfile>" );
    exit( 1 );
  }
  
  // setup MDA reader and read input header
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() )
  {
    cerr << "Cannot read file header\n";
    exit( 1 );
  }
  DataType	type= reader.getType();
  CoordinateVector	dim= reader.getDim();
  unsigned int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSizeIn= reader.getScanlineSize();
  void *scanlineIn;
  
  unsigned int numOutChannels= channels.vec.size();
  if( numOutChannels== 0 )
    // no channels specified? use the first 3 channels!
    for( numOutChannels= 0 ; numOutChannels< 3 ; numOutChannels++ )
      channels.vec.push_back( numOutChannels );
  
  // check that we don't have too many channels
  if( numOutChannels!= 3 )
  {
    cerr << "Can only deal with exactly 3 channels\n";
    exit( 1 );
  }
  // check that all channels in the channel list are vaild
  for( i= 0 ; i< numOutChannels ; i++ )
    if( channels.vec[i]>= numChannels )
    {
      cerr << "Input MDA only has " << numChannels << " channels\n";
      exit( 1 );
    }
  
  // setup output stream
  filebuf bmpFile;
  bmpFile.open( argv[argc-1], ios::out );
  ostream image( &bmpFile );
  if( !image.good() )
  {
    cerr << "Error opening output image file " << argv[argc-1] << endl;
    exit( 1 );
  }
  // byte alignment padding for scanline
  int padding;
  if( (dim.vec[0]*3) % 4 )
    padding= 4- (dim.vec[0]*3)%4;
  else
    padding= 0;
  int paddingData= 0;
  imageHeader.biWidth= dim.vec[0];
  imageHeader.biHeight= numScanlines;
  imageHeader.bfSize=
    (imageHeader.biWidth*3+padding)*imageHeader.biHeight + 54;

  // swap endianness on Macs and other big-endian machines.
  // this is like pulling teeth due to the different data types in the header.
  if( cpuByteOrder()!= LittleEndian )
  {
    swapByteOrder( &(imageHeader.bfType), 2 );
    swapByteOrder( &(imageHeader.bfSize), 4 );
    swapByteOrder( &(imageHeader.bfReserved), 4 );
    swapByteOrder( &(imageHeader.bfOffBits), 4 );
  
    swapByteOrder( &(imageHeader.biSize), 4 );
    swapByteOrder( &(imageHeader.biWidth), 4 );
    swapByteOrder( &(imageHeader.biHeight), 4 );
    swapByteOrder( &(imageHeader.biPlanes), 2 );
    swapByteOrder( &(imageHeader.biBitCount), 2 );
    swapByteOrder( &(imageHeader.biCompression), 4 );
    swapByteOrder( &(imageHeader.biSizeImage), 4 );
    swapByteOrder( &(imageHeader.biXPelsPerMeter), 4 );
    swapByteOrder( &(imageHeader.biYPelsPerMeter), 4 );
    swapByteOrder( &(imageHeader.biClrUsed), 4 );
    swapByteOrder( &(imageHeader.biClrImportant), 4 );
  }
  image.write( ((char *)&imageHeader)+2, 54 );
  
  unsigned char *scanlineOut= new unsigned char[dim.vec[0]*numChannels];
  
  // copy image data
  for( i= 0 ; i< numScanlines ; i++ )
  {
    // read scanline and convert to unsigned byte
    scanlineIn= reader.readScanline();
    typeConvert( scanlineIn, type, scanlineOut, UByte,
		 dim.vec[0]*numChannels );
    
    // now copy all the channels for all pixels in the scanline
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      image.put( scanlineOut[j*numChannels + channels.vec[2]] );
      image.put( scanlineOut[j*numChannels + channels.vec[1]] );
      image.put( scanlineOut[j*numChannels + channels.vec[0]] );
    }
    // output padding
    if( padding )
      image.write( (char *)&paddingData, padding );
  }
  
  reader.disconnect();
  bmpFile.close();
  
  return 0;
}
