// ==========================================================================
// $Id: mda-topfm.C 79 2007-05-16 05:59:14Z agr $
// convert MDA stream to pfm file
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: agr (Allan Rempel)
// Email:   agr@cs.ubc.ca
// ==========================================================================

#include <stdio.h>
#include <string.h>

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
    unsigned long i;

    CommandlineParser parser;

    // setup options
    ChannelList channels;
    ChannelListOption channelOption( channels );
    parser.registerOption( & channelOption );

    // parse options
    int index= 1;
    if( ! parser.parse( index, argc, argv ) || index != argc - 1 )
    {
        parser.usage( argv[0], "<options> <pfm outfile>" );
        exit( 1 );
    }

    // setup MDA reader and read input header
    MDAReader reader;
    reader.connect( cin );
    if( ! reader.readHeader() )
    {
        cerr << "Cannot read file header\n";
        exit( 1 );
    }
    DataType type = reader.getType();
    CoordinateVector dim = reader.getDim();
    unsigned int numChannels = reader.getNumChannels();
    unsigned long numScanlines = reader.getNumScanlinesLeft();
    void *scanlineIn;

    unsigned int numOutChannels = channels.vec.size();
    if ( numOutChannels == 0 )
    // no channels specified? use the first 3 channels!
    for ( numOutChannels = 0 ; numOutChannels < 3 ; numOutChannels++ )
        channels.vec.push_back( numOutChannels );

    // check that we don't have too many channels
    if ( numOutChannels != 3 )
    {
        cerr << "Can only deal with exactly 3 channels\n";
        exit( 1 );
    }
    // check that all channels in the channel list are valid
    for ( i = 0 ; i < numOutChannels ; i++ )
        if ( channels.vec[i] >= numChannels )
        {
            cerr << "Input MDA only has " << numChannels << " channels\n";
            exit( 1 );
        }

    // setup output stream
    FILE * outfile;
    if ( ! strcmp( argv[argc - 1], "-" ) ) outfile = stdout;
    else outfile = fopen( argv[ argc - 1], "w" );
    if ( ! outfile )
    {
        cerr << "Error opening output image file " << argv[argc - 1] << endl;
        exit( 1 );
    }

    fprintf( outfile, "PF\n%d %d\n%f\n", dim.vec[0], numScanlines,
             cpuByteOrder() == LittleEndian ? -1.0 : 1.0 );

    float * scanlineOut = new float[dim.vec[0] * numChannels];

    // copy image data
    for ( i = 0 ; i < numScanlines ; i++ )
    {
        // read scanline and convert to unsigned byte
        scanlineIn = reader.readScanline();
        typeConvert( scanlineIn, type, scanlineOut, Float,
		 dim.vec[0] * numChannels );
        fwrite( (unsigned char *) scanlineOut, 4, dim.vec[0] * numChannels,
                outfile );
        // output padding - not necessary because floats are 4 bytes
    }

    reader.disconnect();
    if ( ! strcmp( argv[ argc - 1], "-" ) ) fclose( outfile );
    delete [] scanlineOut;

    return 0;
}
