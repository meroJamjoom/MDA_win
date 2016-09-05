// ==========================================================================
// $Id: mda-lut.C 319 2009-05-27 21:17:01Z heidrich $
// apply a lookup table to an MDA file
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
  channels.vec.push_back( 0 );
  ChannelListOption channelOption( channels );
  parser.registerOption( &channelOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0] );
    exit( 1 );
  }
  
  // read the lookup table
  MDAReader reader;
  reader.connect( argv[index] );
  if( !reader.readHeader() )
  {
    cerr << "Cannot read LUT file\n";
    exit( 1 );
  }
  int lutSize= reader.getDim().vec[0];
  if( reader.getDim().vec.size()!= 1 || lutSize> 65536 )
  {
    cerr << "Lookup table has to be 1D, with at most 65536 entries\n";
    exit( 1 );
  }
  DataType	typeLUT= reader.getType();
  unsigned int	numChannelsLUT= reader.getNumChannels();
  char *lut= new char[reader.getScanlineSize()];
  memcpy( lut, reader.readScanline(), reader.getScanlineSize() );
  reader.disconnect();
  
  // read input header
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
  unsigned short *scanlineIn= new unsigned short[dim.vec[0]*numChannels];
  
  
  // check if all channels in the channel list are vaild
  unsigned int numOutChannels= channels.vec.size();
  for( i= 0 ; i< numOutChannels ; i++ )
    if( channels.vec[i]>= numChannels )
    {
      cerr << "Input MDA only has " << numChannels << " channels\n";
      exit( 1 );
    }
  if( numOutChannels== 1 )
    for( i= 1 ; i< numChannelsLUT ; i++ )
      channels.vec.push_back( channels.vec[0] );
  else
    if( numOutChannels!= numChannelsLUT )
    {
      cerr << "Number of provided channels does not match #LUT channels\n";
      exit( 1 );    
    }
  numOutChannels= numChannelsLUT;
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numOutChannels, typeLUT ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  unsigned long scanlineSizeOut= writer.getScanlineSize();
  char *scanlineOut= new char[scanlineSizeOut];
  
  // actually copy the data
  unsigned int bytesPerValue= dataTypeSizes[typeLUT];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    // LUT indices are in short
    typeConvert( reader.readScanline(), type, scanlineIn, UShort,
		 dim.vec[0]*numChannels );
    
    // for each pixel in the scanline, perform the lookup
    for( j= 0 ; j< dim.vec[0] ; j++ )
      for( k= 0 ; k< numOutChannels ; k++ )
	memcpy( scanlineOut + (j*numOutChannels + k)*bytesPerValue,
		lut+(((int)scanlineIn[j*numChannels+channels.vec[k]]*
		      (lutSize-1))/65535*numOutChannels+k)*bytesPerValue,
		bytesPerValue );
    
    // write output scanline
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
