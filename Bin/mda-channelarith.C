// ==========================================================================
// $Id: mda-channelarith.C 260 2008-12-30 09:21:39Z heidrich $
// perform per-pixel arithmetic operatiosn on channels
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
#include <time.h>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Array/MDAFileIO.hh"

#if defined(_WIN32) || defined(_WIN64)
#define srand48(num) srand(num)
#endif

using namespace MDA;
using namespace std;

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  int l, m, n;
  
  CommandlineParser parser;
  
  // setup options
  DataType	outType= Double;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  // variables
  EXPR::ExpressionSequence vars;
  ExpressionOption varOpt( vars,
       "\tdefinitions of temporary variables (a sequence of expressions)\n",
			   "--variables", "-v" );
  parser.registerOption( &varOpt );
  
  // whether or not to initialize the random number generator
  bool random= true;
  BoolOption randomOpt( random,
			"\twhether or not to initialize the random number"
			"generator differently\n\tfor each call\n",
			"--random", NULL, "--deterministic", NULL );
  parser.registerOption( &randomOpt );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "[options] <expression sequence>" );
    exit( 1 );
  }
  
  // initialize random number generator if desired
  if( random )
    srand48( (long)clock() );
  
  // get pixel value
  EXPR::ExpressionSequence color= EXPR::parse( argv[argc-1] );
  if( color.size()== 0 )
  {
    parser.usage( argv[0], "[options] <expression sequence>" );
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
  DataType	inType= reader.getType();
  CoordinateVector	dim= reader.getDim();
  int	numChannelsIn= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long inScanlineSize= reader.getScanlineSize();
  char *inScanline;
  
  int numVars= vars.size();
  int numChannelsOut= color.size();
  
  // make sure all register indices are valid
  // first in the variable definitions
  for( l= 0 ; l< numVars ; l++ )
    if( (m= vars[l]->getMaxVariable( true ))>= numChannelsIn ||
	(n= vars[l]->getMaxVariable( false ))>= l )
    {
      cerr << "Using input channels up to " << m
	   << " and temp variables up to " << n << " in variable " << l << endl
	   << "  (only channels 0.." << numChannelsIn-1
	   << " and variables 0.." << l-1 << " are valid at this point)\n";
      exit( 1 );
    }
  // then in the pixel value expressions
  for( l= 0 ; l< numChannelsOut ; l++ )
    if( (m= color[l]->getMaxVariable( true ))>= numChannelsIn ||
	(n= color[l]->getMaxVariable( false ))>= numVars )
    {
      cerr << "Using input channels up to " << m
	   << " and temp variables up to " << n << " in channel " << l << endl
	   << "  (only channels 0.." << (int)numChannelsIn-1
	   << " and variables 0.." << numVars-1
	   << " are valid at this point)\n";
      exit( 1 );
    }
  
  
  // setup output MDA object
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dim, numChannelsOut, outType ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  assert( numScanlines== writer.getNumScanlinesLeft() );
  unsigned long outScanlineSize= writer.getScanlineSize();
  char *outScanline= new char[outScanlineSize];
  
  double *inPixel= new double[numChannelsIn+numVars];
  double *outPixel= new double[numChannelsOut];
  
  // actually copy the data
  unsigned int bytesPerValue= dataTypeSizes[outType];
  for( i= numScanlines ; i> 0 ; i-- )
  {
    inScanline= (char *)reader.readScanline();
    
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      // convert one pixel to double
      typeConvert( inScanline+j*numChannelsIn*dataTypeSizes[inType], inType,
		   inPixel, Double, numChannelsIn );
      
      // compute the value of all variables
      for( k= 0 ; k< numVars ; k++ )
	inPixel[numChannelsIn+k]=
	  vars[k]->eval( inPixel, inPixel+numChannelsIn );
      
      // evaluate all channels
      for( k= 0 ; k< numChannelsOut ; k++ )
	outPixel[k]= color[k]->eval( inPixel, inPixel+numChannelsIn );
      
      // convert back to target data type
      typeConvert( outPixel, Double,
		   outScanline+j*numChannelsOut*dataTypeSizes[outType],outType,
		   numChannelsOut );
    }
    
    // write scanline out
    writer.writeScanline( outScanline );
  }
  reader.disconnect();
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
