// ==========================================================================
// $Id: mda-reduce.C 319 2009-05-27 21:17:01Z heidrich $
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

namespace MDA {
/** the parsing function, specialized for char */
	template<> bool

		  ScalarOption<char>::parse(int &index, int argc, char *argv[])

	{
			if (index > argc - 1)
				return false;

			value = argv[index++][0];
#if defined (_WIN32) || defined (_WIN64)
			if (value< min_s || value> max_s)
#else
			if (value< min || value> max)
#endif
			{
#ifdef DEBUG
				cerr << "Value " << value << " out of range " << min << "..." << max
					<< " for option " << argv[index - 2] << endl;
#endif
				return false;
			}

			return true;

		}
}

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

  // output to MDA format, or just write values in ASCII
  bool outputMDA= false;
  BoolOption outputMDAOption( outputMDA,
	"\twrite a scalar-valued MDA-file, or just echo results to console?\n",
			      "--output-mda", "-mda",
			      "--output-direct", "-d" );
  parser.registerOption( &outputMDAOption );
  
  // a separator string
  char separator= '\n';
  ScalarOption<char> sepOption( separator,
				"\tcharacter used to separate results\n",
				"--separator", "-sep", -128, 127 );
  parser.registerOption( &sepOption );
  
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
  outType= outType== UndefinedType ? inType : outType;
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
	   <<  reduction[i]->getMaxVariable( true )
	   << " and temp variables up to "
	   <<  reduction[i]->getMaxVariable( false )
	   << " in temp variable " << i
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
  
  
  double *result= new double[numChannelsOut];
  if( pre.size()== 0 )
    for( i= 0 ; i< numChannelsOut ; i++ )
      result[i]= 0.0;
  else
    for( k= 0 ; k< numChannelsOut ; k++ )
      result[k]= pre[k]->eval(  NULL, NULL );
  
  
  // read data py scanline, and apply reduction ops to all pixels
  double *inScanline= new double[dim.vec[0]*numChannelsIn];
  double *pixel= new double[numChannelsOut];
  for( i= 0 ; i< numScanlines ; i++ )
  {
    // read next scanline, and convert to double
    typeConvert( (char *)reader.readScanline(), inType,
	 inScanline, Double, dim.vec[0]*numChannelsIn );
	
    
    // apply reduction operator to all pixels
    for( j= 0 ; j< dim.vec[0] ; j++ )
    {
      // apply the appropriate reduction operator to each channel
      for( k= 0 ; k< numChannelsOut ; k++ )
      {
	pixel[k]= reduction[k]->eval( inScanline+j*numChannelsIn, result );
      }
      // copy pixel values back to result buffer
      memcpy( result, pixel,
	      numChannelsOut*sizeof(double) );
    }
  }
  
  
  // apply post operator to all pixels
  if( post.size()> 0 )
    for( j= 0 ; j< numChannelsOut ; j++ )
    {
	pixel[j]= post[j]->eval( NULL, result );
	memcpy( result, pixel, numChannelsOut*sizeof(double) );
    }


  if( outputMDA )
  {
    // setup output MDA object
    CoordinateVector dimOut;
    dimOut.vec.push_back( 1 );
    MDAWriter writer;
    writer.connect( cout );
    if( !writer.writeHeader( dimOut, numChannelsOut, outType ) )
	
    {
      cerr << "Cannot write file header\n";
      exit( 1 );
    }
    unsigned long outScanlineSize= writer.getScanlineSize();
    char *outScanline= new char[outScanlineSize];

    // convert a scanline from double to outType
   typeConvert( result, Double,
		 outScanline, outType,
		 numChannelsOut );
	
    
    // write scanline out
    writer.writeScanline( outScanline );
    
    reader.disconnect();
    if( !writer.disconnect() )
    {
      cerr << "Error while processing data\n";
      exit( 1 );
    }

    delete[] outScanline;
  }
  else
  {
    cout << result[0];
    for( i= 1; i< numChannelsOut ; i++ )
      cout << separator << result[i];
    cout << endl;
  }

  delete[] pixel;
  delete[] inScanline;
  delete[] result;
  
  return 0;
}
