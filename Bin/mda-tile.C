// ==========================================================================
// $Id: mda-tile.C 144 2007-07-20 22:46:40Z heidrich $
// tile an image
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

#include <iostream>
#include <string.h>

#include "MDA/Array/MDAFileIO.hh"
#include "MDA/Array/Array.hh"
#include "MDA/Array/Boundary.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "<PaddingVector>\n"

typedef float TYPE;


int
main( int argc, char *argv[] )
{
  int	i, j, k;
  
  CommandlineParser parser;
  
  CoordinateVector padding;
  padding.vec.push_back(50);
  padding.vec.push_back(25);
  padding.vec.push_back(100);
  padding.vec.push_back(75);
  
  CoordinateOption paddingOpt( padding,
                       "\t padding"
                       "(if size<=0, default)\n",
                       "--radius", "-r" );
  
  parser.registerOption( &paddingOpt );
  
  
  BoundaryMethod boundary= Cyclic;
  BoundaryOption boundaryOpt( boundary );
  parser.registerOption( &boundaryOpt );
  
  int index= 1;
  if( !parser.parse( index, argc, argv ) ) {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // read MDA stream

  Array<TYPE> array;

  if( !array.read() )  {
    cerr << argv[0] << ": Cannot read input MDA stream!\n\n";
    exit( 1 );
  }
  CoordinateVector srcDim= array.getDimension();
  CoordinateVector dstDim= srcDim;
  
  if (srcDim.vec.size()*2 != padding.vec.size()){
    cerr << argv[0] << ": wrong number of paddings!\n\n";
    exit( 1 );
  }
  int dstSize= 1;
  for( i= 0 ; i< dstDim.vec.size() ; i++ ) {
    dstDim.vec[i]+= (padding.vec[i*2]+padding.vec[i*2+1]);
    dstSize*= dstDim.vec[i];
  }
  
  
  int n=array.getNumChannels();
  int c=dstDim.vec[0];
  //n=1;

  MDAWriter writer;
  writer.connect( cout );
  if( !writer.writeHeader( dstDim, n, Float)){//array.getNativeType() ) ){
    cerr << "Cannot write file header\n";
    exit( 1 );
  }
#if defined (_WIN32) || defined (_WIN64)
  TYPE *buf=new TYPE[n]; 
#else
  TYPE *buf[n]; 
#endif
  for (i=0; i<n; i++){

#if defined (_WIN32) || defined (_WIN64)
	  TYPE* buf = new TYPE[dstSize];
	 
	fetchChannel<TYPE>(&buf[i], &(*array[i])[0],
		srcDim, padding, boundary, 0.5f );
#else
    buf[i]=new TYPE[dstSize];

    fetchChannel<TYPE>( buf[i], &(*array[i])[0], 
			srcDim, padding, boundary, 0.5 );
#endif
  }
#if defined (_WIN32) || defined(_WIN64)
  TYPE* tmpbuf=new TYPE[n*c] , *p;
#else
  TYPE tmpbuf[n*c], *p;
#endif
  while( writer.getNumScanlinesLeft() ){
    p=tmpbuf;
    for (j=0; j<c; j++)
      for (i=0; i<n; i++)
#if defined (_WIN32) || defined (_WIN64)
		  *p++=buf[i]++;
#else
	*p++=*buf[i]++;
#endif
    writer.writeScanline( tmpbuf );
  }
  
  if( !writer.disconnect() ){
    cerr << "Error while processing data\n";
    exit( 1 );
  }
  
  return 0;
}
