#include <iostream>

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Color/Spectrum.hh"
#include "MDA/Color/ColorSpaceFactory.hh"

using namespace std;
using namespace MDA;

int
main( int argc, char *argv[] )
{
  CommandlineParser parser;
  
  // setup options
  ColorSpaceFactory colorFac;
  colorFac.registerOptions( parser );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0],
		  "<options> [range]\n"
		  "\nConvert a spectrum to a color in a specified space.\n"
		  "The range indicates the frequency range provided in the\n"
		  "input MDA file representing the spectral samples\n" );
    exit( 1 );
  }

  Array<float> spec;
  IntRange sRange;
  istringstream rangeString( argv[index] );
  rangeString >> sRange;
  
  Vector XYZ( 3 );
  spec.read( cin );
  Spectrum s( spec, 0, sRange, sRange );
  s.toXYZ( XYZ, false );
  
  // normalize XYZ intensity, i.e. convert to xyz
  XYZ/= XYZ.norm( 1 );
  
  Vector tri( 3 );
  ColorSpace *colorSpace= colorFac.makeColorSpace();
  colorSpace->fromXYZ( XYZ, tri );
  
  cout << tri[0] << ' '  << tri[1] << ' '  << tri[2] << endl;
  
  return 0;
}
