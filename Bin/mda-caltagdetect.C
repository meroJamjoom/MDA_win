
#include <string>

#include <MDA/Base/CommandlineParser.hh>

#include "MDA/CameraCalibration/CALTagPattern.hh"
#include "MDA/CameraCalibration/PointCorrespondence.hh"

#define USAGE_TEXT "<options> <imagefile.mda\n  Output detected points as <dest> given a pattern specification file <pattern>\n  Both \"-p\" and \"-d\" must be specified."

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

  char* dest = NULL;
  FileOption destOpt( dest,
                      "\tDestination XML file for correspondences (include .xml suffix)\n",
                      "--dest", "-d" );
  parser.registerOption( &destOpt );
  
  //!! this prefix needs to be exposed on the commandline...
  //no stringoption in commandlineparser
  //this has to be "MDA.CALTag" not "MDA.CALTag.corr" because otherwise the
  //configure() call will fail when looking for "rows" and "columns" attribs
  //the ".corr" path extension is automatically added
  const char *mdPrefix= "MDA.CALTag";

  // do the parsing
  int index= 1;
  if( !parser.parse( index, argc, argv ) )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // some consistency checks
  errorCond( pattern!=NULL, "Must specify a pattern file" );
  errorCond( dest!=NULL, "Must specify a destination file" );
  
  // Create an empty CALTag pattern and metadata structure
  CALTagPattern p;
  MetaData md;
  // read existing pattern if one has been specified; and initialize pattern
  md.read( pattern );
  p.configure( md, mdPrefix );
  
  Array<double> inputImage;
  if( !inputImage.read() )
  {
    cerr << argv[0] << ": Cannot read input MDA stream!\n\n";
    exit( 1 );
  }
  p.setImage( &inputImage );
  if( p.detect() )
  {
    p.constructCorrespondences();
    p.writeCorrespondences( md, mdPrefix, dest );
    cout << "Detected " << p.getCorrespondences().size() << " points" << endl;
  }
  else
  {
    cout << "Failed to detect any points" << endl;
  }
  
  return 0;
}
