// ==========================================================================
// $Id: mda-mergechannels.C 199 2008-04-08 03:58:26Z heidrich $
// merge selected channels from multiple MDA files
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


#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

template class vector<MDAReader *>;

#define USAGE_TEXT "<options> <MDA file> [<MDA file>...]"

int
main( int argc, char *argv[] )
{
  unsigned long i, l;
  int j, k;
  vector<MDAReader *> readers;
  vector<int> channels;
  vector<DataType> types;
  vector<unsigned> bytesPerValue;
  int numInputs= 0;
  int numOutChannels= 0;
  int numChannels;
  
  // set up options
  CommandlineParser parser;
  
  // output type
  DataType	outType= UndefinedType;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index>= argc )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  CoordinateVector dim;
  // all remaining command line parameters are MDA files
  for( numInputs= 0 ; index< argc ; index++, numInputs++ )
  {
    // creat an MDA reader for each file
    readers.push_back( new MDAReader );
    if( !readers[numInputs]->connect( argv[index] ) ||
	!readers[numInputs]->readHeader() )
    {
      cerr << "Cannot read file header! ("
	   << argv[index] << ")\n";
      exit( 1 );
    }
    
    if( numInputs== 0 )
    {
      // the first file determines dimension and (possibly) type
      outType= (outType== UndefinedType) ? readers[0]->getType() : outType;
      dim= readers[0]->getDim();
    }
    else
    {
      // for the remaining MDA streams, check if the parameters match
      if( dim.vec!= readers[numInputs]->getDim().vec )
      {
	cerr << "Dimension of MDA files does not match! ("
	     << argv[index] << ")\n";
	exit( 1 );
      }
    }
    
    // store data type and size for that input channel
    types.push_back( readers[numInputs]->getType() );
    bytesPerValue.push_back( dataTypeSizes[types[numInputs]] );
    
    // store number of channels for the latest file
    numChannels= readers[numInputs]->getNumChannels();
    channels.push_back( numChannels );
    numOutChannels+= numChannels;
  }
  
  if( numInputs== 0 )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numOutChannels, outType ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  unsigned long numScanlines= writer.getNumScanlinesLeft();
  unsigned long scanlineSizeOut= writer.getScanlineSize();
  char *scanlineOut= new char[scanlineSizeOut];
  unsigned outBytesPerValue= dataTypeSizes[outType];
  
  // actually copy the data
  for( i= numScanlines ; i> 0 ; i-- )
  {
    for( j= k= 0 ; j< numInputs ; j++ )
    {
      char *scanlineIn= (char *)readers[j]->readScanline();
      
      // for each pixel in the scanline, copy all channels
      for( l= 0 ; l< dim.vec[0] ; l++ )
	typeConvert( scanlineIn + l*channels[j]*bytesPerValue[j],
		     types[j],
		     scanlineOut + (l*numOutChannels + k)*outBytesPerValue,
		     outType,
		     channels[j] );
      
      k+= channels[j];
    }
    
    writer.writeScanline( scanlineOut );
  }
  
  for( i= 0 ; i< numInputs ; i++ )
    readers[i]->disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
