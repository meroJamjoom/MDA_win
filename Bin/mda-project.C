// ==========================================================================
// $Id: mda-project.C 319 2009-05-27 21:17:01Z heidrich $
// project an N-dimensional array down to an N-1 dimesnional one
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

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  unsigned long i, j, k;
  
  CommandlineParser parser;
  
  // setup options
  DataType	outType= UndefinedType;
  TypeOption	typeOption( outType );
  parser.registerOption( &typeOption );
  
  // initalization of pixel values
  EXPR::ExpressionSequence pre;
  ExpressionOption preOpt( pre,
	   "\tinitialization of pixel values (a sequence of expressions)\n",
			   "--pre" );
  parser.registerOption( &preOpt );
  
  // initalization of pixel values
  EXPR::ExpressionSequence post;
  ExpressionOption postOpt( post,
	   "\tpost-processing of pixel values (a sequence of expressions)\n",
			   "--post" );
  parser.registerOption( &postOpt );
  
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "[options] <expression sequence>" );
    exit( 1 );
  }
  
  // get reduction operators
  EXPR::ExpressionSequence reduction= EXPR::parse( argv[argc-1] );
  int numChannelsOut= reduction.size();
  if( numChannelsOut== 0 )
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
  outType= (outType== UndefinedType) ? inType : outType;
  CoordinateVector	dim= reader.getDim();
  int numChannelsIn= reader.getNumChannels();
  unsigned long numScanlines= reader.getNumScanlinesLeft();
  unsigned long inScanlineSize= reader.getScanlineSize();
  
  
  // make sure pre, post, and reduction all have the same number of channels
  if( (pre.size()!= 0 && pre.size()!= numChannelsOut) ||
      (post.size()!= 0 && post.size()!= numChannelsOut) )
  {
    cerr << "Number of pre and post expressions needs to match number"
	 << " of reduction operators\n";
    exit( 1 );
  }
  // make sure all variable indices are valid
  // first in the pre sequence
  for( i= 0 ; i< pre.size() ; i++ )
    if( pre[i]->getMaxVariable( true )>= 0 ||
	pre[i]->getMaxVariable( false )>= 0 )
    {
      cerr << "Must not use variables in pre expressions\n";
      exit( 1 );
    }
  // then in the reduction operators
  for( i= 0 ; i< numChannelsOut ; i++ )
    if( reduction[i]->getMaxVariable( true )>= numChannelsIn ||
	reduction[i]->getMaxVariable( false )>= numChannelsOut )
    {
      cerr << "Using input channels up to "
	   << reduction[i]->getMaxVariable( true )
	   << " and temp variables up to "
	   << reduction[i]->getMaxVariable( false ) <<
	" in temp variable " << i
	   << "\n  (only channels 0.." << numChannelsIn-1
	   << " and variables 0.." << numChannelsOut-1 << " are valid)\n";
      exit( 1 );
    }
  // finally, in the post expressions
  for( i= 0 ; i< post.size() ; i++ )
    if( post[i]->getMaxVariable( true )>= 0 ||
	post[i]->getMaxVariable( false )>= numChannelsOut )
    {
      cerr << "Using input channels and/or temp variables up to "
	   << post[i]->getMaxVariable( false )
	   << " in post epxression " << i << endl
	   << "  (no input and only variables 0.." << numChannelsOut-1
	   << " are valid)\n";
      exit( 1 );
    }
  
  
  // setup a buffer for one slice and initialize with pre expressions
  unsigned long numSliceScanlines= numScanlines/dim.vec[dim.vec.size()-1];
  unsigned long numSlicePixels= numSliceScanlines*dim.vec[0];
  
  double *slice= new double[numSlicePixels*numChannelsOut];
  if( pre.size()== 0 )
    for( i= 0 ; i< numSlicePixels*numChannelsOut ; i++ )
      slice[i]= 0.0;
  else
  {
    for( k= 0 ; k< numChannelsOut ; k++ )
      slice[k]= pre[k]->eval(  NULL, NULL );
    for( i= 1 ; i< numSlicePixels ; i++ )
      for( j= 0 ; j< numChannelsOut ; j++, k++ )
	slice[k]= slice[j];
  }
  
  
  
  // read data by scanline, and apply reduction ops to all pixels
  double *inScanline= new double[dim.vec[0]*numChannelsIn];
  double *pixel= new double[numChannelsOut];
  for( i= 0 ; i< numScanlines ; i++ )
  {
    double *sliceScanline=
      slice + ((i%numSliceScanlines) * dim.vec[0]*numChannelsOut);
    
    // read next scanline, and convert to double
    typeConvert( (char *)reader.readScanline(), inType,
		 inScanline, Double, dim.vec[0]*numChannelsIn );
    
    
    // apply reduction operator to all pixels
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      // apply the appropriate reduction operator to each channel
      for( k= 0 ; k< numChannelsOut ; k++ )
      {
	pixel[k]= reduction[k]->eval( inScanline+j*numChannelsIn,
				      sliceScanline+j*numChannelsOut );
      }
      // copy pixel values back to slice buffer
      memcpy( sliceScanline+j*numChannelsOut, pixel,
	      numChannelsOut*sizeof(double) );
    }
  }
  
  
  // apply post operator to all pixels
  if( post.size()> 0 )
    for( i= 0 ; i< numSlicePixels ; i++ )
    {
      for( j= 0 ; j< numChannelsOut ; j++ )
	pixel[j]= post[j]->eval( NULL, slice+i*numChannelsOut );
      memcpy( slice+i*numChannelsOut, pixel, numChannelsOut*sizeof(double) );
    }
  
  // setup output MDA object
  CoordinateVector dimOut;
  for( i= 0 ; i< dim.vec.size()-1 ; i++ )
    dimOut.vec.push_back( dim.vec[i] );
  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dimOut, numChannelsOut, outType ) )
  {
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
  unsigned long numScanlinesOut= writer.getNumScanlinesLeft();
  unsigned long outScanlineSize= writer.getScanlineSize();
  char *outScanline= new char[outScanlineSize];
  
  // write out the data
  for( i= 0 ; i< numScanlinesOut ; i++ )
  {
    // convert a scanline from double to outType
    typeConvert( slice+i*dim.vec[0]*numChannelsOut, Double,
		 outScanline, outType,
		 dim.vec[0]*numChannelsOut );
    
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
