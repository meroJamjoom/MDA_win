// ==========================================================================
// $Id: mda-resize.C 642 2010-03-09 20:10:31Z heidrich $
// resize an MDA stream
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


#include "MDA/Config.hh"
#include "MDA/Array/Array.hh"
#include "MDA/Resampling/Resampling.hh"
#include "MDA/GeometricTransform/GeometricTransform.hh"

using namespace MDA;
using namespace std;

char usageText[BUFFER_SIZE]= "[<options>...] <new dimensions>";


/** the actual resizing is encapuslated in this template function, so
    we can choose the type for computations based on the desired
    output type */
template<class T>
void
doScale( char *progName, CoordinateVector outDim,
	 ResampleMode resampleMode, DataType outType )
{
  // create resampler and configure it
  ResamplerFactory<T> factory;
  Resampler<T> *resampler= factory.create( resampleMode );
  resampler->setBoundaryMethod( Clamp );
  
  // Scaling object
  Scaling<T> scaler( resampler );
  scaler.setScaling( Vector( outDim.vec.size(), 1.0 ) );
  
  // read MDA stream
  Array<T> inArray;
  if( !inArray.read() )
  {
    cerr << progName << ": Cannot read input MDA stream!\n\n";
    exit( 1 );
  }
  
  // create output array
  Array<T> outArray( outDim );
  outArray.setNativeType( inArray.getNativeType() );
  
  // apply scaling
  scaler.GeometricTransformation<T>::apply( inArray, outArray, true );
  
  // write MDA stream
#if defined (_WIN32) || defined (_WIN64) 
  cout.seekp(ios_base::beg);
  outArray.write( cout, outType );
#else
  outArray.write(cout, outType);
#endif 
}


int
main( int argc, char *argv[] )
{
  // setup options
  CommandlineParser parser;
  
  // resampler mode
  ResampleMode resampleMode= CubicResampling;
  ResamplerOption resampleOpt( resampleMode );
  parser.registerOption( &resampleOpt );
  
  // output type
  DataType	outType= UndefinedType;
  TypeOption	typeOption( outType,
			    "\tOutput type (default: same as input)\n" );
  parser.registerOption( &typeOption );
  
  // process commandline options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], usageText );
    exit( 1 );
  }
  
  CoordinateVector outDim;
  istringstream dimString( argv[index] );
  dimString >> outDim;
  if( dimString.fail() )
  {
    parser.usage( argv[0], usageText );
    exit( 1 );
  }
  
  if( outType== Double )
    doScale<double>( argv[0], outDim, resampleMode, outType );
  else
    doScale<float>( argv[0], outDim, resampleMode, outType );
  
  return 0;
}
