// ==========================================================================
// $Id: mda-newcaltag.C 199 2008-04-08 03:58:26Z heidrich $
// create a new CALTag pattern (postscript and metadata file)
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010, UBC
//
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#include <string>

#include <MDA/Base/CommandlineParser.hh>

#include "MDA/CameraCalibration/CALTagPattern.hh"
#include "MDA/CameraCalibration/PointCorrespondence.hh"

#define USAGE_TEXT "<options> <file base>\n  Output new CALTag pattern as <file base>.ps and <file base>.xml\n  (at least one of \"-p\" or \"-r\" must be specified)"

using namespace std;
using namespace MDA;

int
main( int argc, char *argv[] )
{
  // set up options
  CommandlineParser parser;
  
  // possible MetaData file
  char *pattern= NULL;
  FileOption patternOpt( pattern,
			 "\tExisting pattern file (metadata file format)\n",
			 "--pattern", "-p" );
  parser.registerOption( &patternOpt );
  
  //!! this prefix needs to be exposed on the commandline...
  const char *mdPrefix= "MDA.CALTag";
  
  // pattern style/layout
  int layout= 0;
  IntOption layoutOpt( layout,
		       "\tPattern layout (0: default, 1: checker, 2: bowtie)\n",
		       "--layout", "-l" );
  parser.registerOption( &layoutOpt );
  
  // start ID
  int startID= 0;
  IntOption startOpt( startID, "\tOffset for used marker IDs\n",
		      "--start-id", NULL );
  parser.registerOption( &startOpt );
  
  // scale
  double scale= -1.0;
  DoubleOption sizeOpt( scale,
			"\tMarker size (in inches), -1.0 means default\n",
			"--scale", "-s" );
  parser.registerOption( &sizeOpt );
  
  // marker re-use
  bool newMarkers= false;
  BoolOption newOpt( newMarkers,
		     "\tWhether or not to re-use existing marker IDs\n",
		     "--fresh-markers", NULL, "--reuse-markers", NULL );
  parser.registerOption( &newOpt );
  
  // are manual markers raw or with CRC?
  bool makeCRC= true;
  BoolOption rawOpt( makeCRC,
	     "\tWhether to use custom markers as is (raw) or add a CRC\n",
		     "--crc-markers", NULL, "--raw-markers", NULL );
  parser.registerOption( &rawOpt );
  
  // custom marker IDs
  CoordinateVector customMarkers;
  CoordinateOption markerOpt( customMarkers,
			      "\tCustom marker numbers (comma separated)\n",
			      "--markers", "-m" );
  parser.registerOption( &markerOpt );
  
  // marker columns, rows
  CoordinateVector grid;
  grid.vec.push_back( 0 );
  grid.vec.push_back( 0 );
  CoordinateOption gridOpt( grid, "\tGrid resolution: columns,rows\n",
			    "--res", "-r" );
  parser.registerOption( &gridOpt );
  
  // do the parsing
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index>= argc )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // some consistency checks
  errorCond( grid.vec.size()== 2, "Wrong grid dimension (not 2D)" );
  errorCond( (grid.vec[0]> 0 && grid.vec[1]> 0) || pattern!= NULL,
	     "Must have at least one of \"-p\" or \"-r\" options" );
  
  // Create an empty CALTag pattern and metadata structure
  CALTagPattern p;
  MetaData md;
  // read existing pattern if one has been specified; and initialize pattern
  if( pattern!= NULL )
  {
    md.read( pattern );
    p.configure( md, mdPrefix );
  }
  
  // then make configuration changes according to commandline options
  p.changeLayout( grid.vec[1], grid.vec[0], scale, layout );
  p.setStart( startID );
  
  // load custom IDs if any are specified, otherwise generate any
  // missing IDs
  if( customMarkers.vec.size()> 0 )
    p.setIDs( customMarkers, makeCRC );
  else
    p.makeIDs( newMarkers );
  
  // create the postscript pattern
  p.makePattern( ((string)(argv[argc-1])+".ps").c_str(), PaperSize( Letter ) );
  
  // create new metadata file
  p.writeConfig( md, mdPrefix );
  md.write( ((string)(argv[argc-1])+".xml").c_str() );
  
  return 0;
}
