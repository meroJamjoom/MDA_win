// ==========================================================================
// $Id: mda-swizzle.C 319 2009-05-27 21:17:01Z heidrich $
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

#include <string.h>
#include <assert.h>

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // setup options
  ChannelList	channels;
  ChannelListOption channelOption( channels );
  parser.registerOption( &channelOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) )
  {
    parser.usage( argv[0] );
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
  char *scanlineIn;
  

  // check if all channels in the channel list are vaild
  unsigned int numOutChannels= channels.vec.size();
  for( i= 0 ; i< numOutChannels ; i++ )
    if( channels.vec[i]>= numChannels )
    {
      cerr << "Input MDA only has " << numChannels << " channels\n";
      exit( 1 );
    }
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numOutChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  unsigned long scanlineSizeOut= writer.getScanlineSize();
  char *scanlineOut= new char[scanlineSizeOut];
  
  
  // actually copy the data
  unsigned int bytesPerValue= dataTypeSizes[type];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    scanlineIn= (char *)reader.readScanline();
    
    // for each pixel in the scanline, copy the required channels
    for( j= 0 ; j< dim.vec[0] ; j++ )
      for( k= 0 ; k< numOutChannels ; k++ )
	memcpy( scanlineOut + (j*numOutChannels + k)*bytesPerValue,
		scanlineIn + (j*numChannels + channels.vec[k])*bytesPerValue,
		bytesPerValue );
    
    writer.writeScanline( scanlineOut );
  }
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
