// ==========================================================================
// $Id: mda-saveseq.C 638 2010-03-08 19:34:24Z heidrich $
// save a sequence of MDA objects into separate files
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
#include <math.h>
#include <stdlib.h>

#include "MDA/Config.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "[<options>] <file name> [<file name>...]\n"

int
main( int argc, char *argv[] )
{
  char	fileName[BUFFER_SIZE];
  char  cmdName[BUFFER_SIZE];
  unsigned long	i;
  
  CommandlineParser parser;
  
  // integer offset for label of first MDA object
  int num= 0;
  IntOption firstOption( num,
			 "\toffset for numbering of output files"
			 " (label for first extracted MDA object)\n",
			 "--first", NULL );
  parser.registerOption( &firstOption );
  
  // increment for the number of output objects
  int incr= 1;
  IntOption incrOption( incr,
			"\tincrement for output numbers\n",
			"--incr", NULL );
  parser.registerOption( &incrOption );
  
  // whether or not to compress the MDA files (default true)
  bool compression= true;
  BoolOption comprOption( compression,
                          "\ttoggle gzip compression of results\n",
                          "--compression", "-c", "--no-compression", "-nc" );
  parser.registerOption( &comprOption );
  
#ifdef DEBUG
  bool verbose= true;
#else
  bool verbose= false;
#endif
  BoolOption verboseOption( verbose,
			    "\tverbose output\n",
			    "--verbose", "-v", "--quiet", "-q" );
  parser.registerOption( &verboseOption );
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index> argc-1 )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  
  // extract as many MDA streams as possible from stdin
  MDAReader reader;
  while( reader.connect( cin ) )
  {
    if( !reader.readHeader() )
      break;
    
    if( index< argc-1 )
      sprintf( fileName, "%s.mda", argv[index++] );
    else
    {
      sprintf( fileName, "%s.%06d.mda", argv[argc-1], num );
      num+= incr;
    }
    sprintf( cmdName, "gzip -9 %s", fileName );
    
    // setup a writer with the exact properties of the input stream
    MDAWriter writer;
    writer.connect( fileName );
    if( !writer.writeHeader( reader.getDim(),
			     reader.getNumChannels(),
			     reader.getType() ) )
    {
      cerr << argv[0] << ": Cannot write file header\n";
      exit( 1 );
    }
    
    // copy all scanlines
    unsigned long numScanlines= reader.getNumScanlinesLeft();
    for( i= 0 ; i< numScanlines ; i++ )
      writer.writeScanline( reader.readScanline() );
    
    // close both reader and writer
    if( !writer.disconnect() )
    {
      cerr << argv[0] << ": Error writing...\n";
      exit( 1 );
    }
    reader.disconnect();
    
    // run file compression if desired
    if( compression )
      system( cmdName );
    
    if( verbose )
      cerr << argv[0] << ": Wrote " << fileName << endl;
  }
  
  return 0;
}
