// ==========================================================================
// $Id: mda-newmda.C 640 2010-03-08 23:27:20Z heidrich $
// generate a new MDA stream initialized with user provided pixel data
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

#include <sstream>
#include <math.h>

#include <time.h>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Resampling/NormalizedDeviceCoordinates.hh"
#include "MDA/Array/MDAFileIO.hh"

#if defined (_WIN32) || defined (_WIN64)
#define srand48(num) srand(num)
#endif

using namespace MDA;
using namespace std;

#define USAGE_TEXT "<options> <dimensions>\n\n\
Create a new MDA stream initialized with user provided pixel data\n"

int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  EXPR::ExpressionSequence vars; 
  ExpressionOption varOpt( vars,
       "\tdefinitions of temporary variables (a sequence of expressions)\n",
			   "--variables", "-v" );
  parser.registerOption( &varOpt );
  
  // pixel value
  EXPR::ExpressionParseTree defaultPixel;
  defaultPixel.makeConstant( 0.0 );
  EXPR::ExpressionSequence color;
  color.push_back( &defaultPixel );
  ExpressionOption colorOpt( color,
			     "\tpixel value (sequence of expressions)\n",
			     "--pixel-value", "-p" );
  parser.registerOption( &colorOpt );
  
  // type
  DataType type= Float;
  
  TypeOption typeOpt( type );
  parser.registerOption( &typeOpt );
  
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
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // initialize random number generator if desired
  if( random )
    srand48( (long)clock() );
  
  // the dimensions of the MDA stream
  CoordinateVector      dim;
  istringstream dimString( argv[index] );
  dimString >> dim;
  if( dimString.fail() )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // process user-provided pixel data for channel info
  int dimensions= dim.vec.size();
  int numChannels= color.size();
  int numVars= vars.size();
  // make sure all register indices are valid
  // first in the variable definitions
  for( i= 0 ; i< numVars ; i++ )
    if( vars[i]->getMaxVariable( true )>= (int)dimensions ||
	vars[i]->getMaxVariable( false )>= (int)i )
    {
      cerr << "Using coordinates up to " << vars[i]->getMaxVariable( true )
	   << " and temp variables up to " << vars[i]->getMaxVariable( false )
	   << " in temp variable " << i << endl
	   << "  (only coordinates up to " << (int)2*dimensions-1
	   << " and variables up to " << (int)i-1
	   << " are valid at this point)\n";
      exit( 1 );
    }
  // then in the pixel value expressions
  for( i= 0 ; i< numChannels ; i++ )
  if (color[i]->getMaxVariable(true) >= 2 * dimensions ||
	  color[i]->getMaxVariable(false) >= numVars) {
#if defined (_WIN32) || defined (_WIN64)
		  {
			  cerr << "Using coordinates up to " << i
			  << " and temp variables up to " << color[i]->getMaxVariable(false)
			  << " in temp variable " << color[i]->getMaxVariable(true) << endl
			  << "  (only coordinates up to " << (int)2 * dimensions - 1
			  << " and variables 0.." << numVars - 1
			  << " are valid at this point)\n";
			  exit(1);
		  }
  
	#else
    {
      cerr << "Using coordinates up to " << j
	   << " and temp variables up to " << color[i]->getMaxVariable( false )
	   << " in temp variable " << color[i]->getMaxVariable( true ) << endl
	   << "  (only coordinates up to " << (int)2*dimensions-1
	   << " and variables 0.." << numVars-1
	   << " are valid at this point)\n";
      exit( 1 );
    }
#endif
  }
  // setup output MDA object

	MDAWriter writer;
  writer.connect( cout );

  if( !writer.writeHeader( dim, numChannels, type ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  unsigned long outScanlineSize= writer.getScanlineSize();
  unsigned long numScanlines= writer.getNumScanlinesLeft();
 
  // create a scanline and registers corresponding to the current coordinates
  
   char *scanline= new char[outScanlineSize];
  double *tmpVariables= new double[numVars];
  double *tempPixel= new double[numChannels];
  
  // initial position (in pixel coord), and NDC mapping
  
  NormalizedDeviceCoordinates ndcMap( dim );
 
  // input variables for expressions (first NDC, then pixel)
  Vector inVariables( dimensions*2, 0.0 );
  // pixel coord. subvector
  Vector pixelPos= inVariables.getSubVector( dimensions, 2*dimensions-1 );
  // NDC subvector
  Vector ndcPos= inVariables.getSubVector( 0, dimensions-1 );
  
  // output the data
  for( i= 0 ; i< numScanlines ; i++ )
  {
	
    // map pixel pos to NDC
    ndcMap.pointFromPixelToNDC( pixelPos, ndcPos );
    
    // now create each pixel in the scanline
    for( j= 0 ; j< dim.vec[0] ; j ++ )
    {
      // coordinate within the scanline
      pixelPos[0]= j;
      ndcPos[0]= ndcMap.pointCompFromPixelToNDC( j, 0 );
      
      // compute the value of the variables, and append them to the end
      // of the register array
      for( k= 0 ; k< numVars ; k++ )
	tmpVariables[k]= vars[k]->eval( &ndcPos[0], tmpVariables );
      
      // then evaluate the per-channel expressions
      for( k= 0 ; k< numChannels ; k++ )
	tempPixel[k]= color[k]->eval( &ndcPos[0], tmpVariables );
      
      // and convert the result to the target type
      typeConvert( tempPixel, Double,
		   scanline+j*numChannels*dataTypeSizes[type], type,
		   numChannels );
    }
    
    // write scanline
    writer.writeScanline( scanline );
    
    // update scanline positions
    for( j= 1 ; j< dimensions ; j++ )
      if( (pixelPos[j]+= 1.0)>= dim.vec[j] )
	pixelPos[j]= 0.0;
      else
	break;
  }
 
  if( !writer.disconnect() )
  {
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  // clean up
  delete [] scanline;
  delete [] tmpVariables;
  delete [] tempPixel;
  
  return 0;
}

