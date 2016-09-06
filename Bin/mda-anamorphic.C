// ==========================================================================
// $Id: mda-anamorphic.C 79 2007-08-17 16:51:14Z bradleyd $
// resize mda file in X dimension so that aspect ratio is 1.78:1, performing 
// cubic interpolation in each scanline
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: bradleyd (Derek Bradley)
// Email:   bradleyd@cs.ubc.ca
// ==========================================================================

#include <sstream>
#include <math.h>

#include "MDA/Array/MDAFileIO.hh"
#include "MDA/Resampling/NormalizedDeviceCoordinates.hh"
#include "MDA/Resampling/Resampling.hh"


using namespace MDA;
using namespace std;

#define USAGE_TEXT "\n\n\
Resize the MDA file in the X dimension so that the aspect\n\
ratio becomes 1.777778:1.  Performs cubic interpolation on\n\
each scanline.\n"


int
main( int argc, char *argv[] )
{
  unsigned long i;

  // commandline parsing
  CommandlineParser parser;
  
  // the resampling filter
  ResampleMode resampleMode= CubicResampling;
  ResamplerOption resamplerOpt( resampleMode );
  parser.registerOption( &resamplerOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], "<options>" );
    exit( 1 );
  }
  
  // create resampler (float precision is sufficient for the char and
  // short arrays we primarily expect in this tool)
  ResamplerFactory<float> factory;
  Resampler<float> *resampler= factory.create( resampleMode );
  resampler->setBoundaryMethod( Clamp );
  
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
  unsigned		dimension= dim.vec.size();
  unsigned int	numChannels= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long scanlineSize= reader.getScanlineSize();
  char *scanline;
  char *scanlineFinal;
  float *scanlineIn;
  float *scanlineOut;
  
  // determine new scanline width
  unsigned long newWidth = (int)(4.0/3.0 * dim.vec[0] + 0.5);
  unsigned long newScanlineSize = newWidth*numChannels*dataTypeSizes[type];
  
  // output dimension is same as input, with first component exchanged
  CoordinateVector dimOut= dim;
  dimOut.vec[0]= newWidth;
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dimOut, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  
  scanlineFinal = new char[newScanlineSize];
  scanlineIn = new float[dim.vec[0]*numChannels];
  scanlineOut = new float[dimOut.vec[0]*numChannels];
  
  // set up the transformation (scaling is implicit in pixel to NDC mappings)
  NormalizedDeviceCoordinates in2NDC( dim );
  NormalizedDeviceCoordinates out2NDC( dimOut );
  // full n x n homog. matrix is a bit overkill, but it is the easiest way
  Matrix pixelXform( dimension+1, dimension+1 );
  multMatrixMatrix( in2NDC.getNDCToPixelMatrix(),
		    out2NDC.getPixelToNDCMatrix(), pixelXform );
  
  // resample the data by scanline (along scanlines only)
  CoordinateVectorIter pos( dim );
  for( pos.begin() ; !pos.isAtEnd() ; pos.incrComp( 1 ) )
  {
    // read scanline
    scanline = (char *)reader.readScanline();
    
    // convert to float
    typeConvert(scanline, type, scanlineIn, Float, dim.vec[0]*numChannels);
    
    // resample each channel
    for( i= 0 ; i< numChannels ; i++ )
      resampler->resampleLinear( scanlineIn+i, dim.vec[0], numChannels,
				 scanlineOut+i, dimOut.vec[0], numChannels,
				 pixelXform[0][0], pixelXform[0][dimension] );
    
    // convert back
    typeConvert( scanlineOut, Float, scanlineFinal, type,
		 dimOut.vec[0]*numChannels);
    
    // write new scanline
    writer.writeScanline( scanlineFinal );
  }
  
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  delete [] scanlineFinal;
  delete [] scanlineIn;
  delete [] scanlineOut;

  return 0;
}
