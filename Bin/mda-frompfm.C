// ==========================================================================
// $Id: mda-frompfm.C  $
// convert pfm file to MDA stream
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: nasarouf (Mushfiqur Rouf)
// Email:   nasarouf@cs.ubc.ca
// ==========================================================================

#include <stdio.h>
#include <string.h>

#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
    unsigned long i;

    CommandlineParser parser;

    // setup options
    //no options!!

    // parse options
    int index= 1;
    if( ! parser.parse( index, argc, argv ) || index != argc - 1 )
    {
        parser.usage( argv[0], "<options> <pfm infile>" );
        exit( 1 );
    }

    // setup input stream
    FILE * infile;
    if ( ! strcmp( argv[argc - 1], "-" ) ) infile = stdin;
    else infile = fopen( argv[ argc - 1], "r" );
    if ( ! infile )
    {
        cerr << "Error opening input image file " << argv[argc - 1] << endl;
        exit( 1 );
    }

    CoordinateVector dim(2); //has to be 2D
    int numChannels=3; //has to be 3!
    int numScanlines;
    float thisCpuByteOrder;
    
    if (fscanf( infile, "PF\n%d %d\n%f\n", 
	       &dim.vec[0], &numScanlines, &thisCpuByteOrder)!=3){
        cerr << "Wrong pfm header: " << argv[argc - 1] << endl;
        exit( 1 );
    }
    dim.vec[1]=numScanlines;

    // setup MDA writer and read input header
    MDAWriter writer;

    writer.connect( cout );
    //  cpuByteOrder() == LittleEndian ? -1.0 : 1.0 );
    if( ! writer.writeHeader(dim,numChannels,Float,
			     thisCpuByteOrder<0?LittleEndian:BigEndian) )
    {
        cerr << "Cannot write file header\n";
        exit( 1 );
    }

    float * scanlineIn = new float[dim.vec[0] * numChannels];

    // copy image data
    for ( i = 0 ; i < numScanlines ; i++ )
    {
        // read scanline and convert to unsigned byte
        fread( (unsigned char *) scanlineIn, 4, dim.vec[0] * numChannels,
                infile );
        //scanlineIn = reader.readScanline();

	writer.writeScanline(scanlineIn);
    }

    writer.disconnect();
    if ( ! strcmp( argv[ argc - 1], "-" ) ) fclose( infile );
    delete [] scanlineIn;

    return 0;
}
