// ==========================================================================
// $Id: mda-fromyuv.C 638 2010-03-08 19:34:24Z heidrich $
// split an n-dimensional MDA stream into multiple slices of dimension n-1
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
#include <string.h>

#include "MDA/Config.hh"
#include "MDA/Base/Range.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;

#define USAGE_TEXT "[<options>] [<file name base>]\n"


void
upsampleLine( unsigned char *src, unsigned char *dst, unsigned long numElem )
{
  for( unsigned long i= 0 ; i< numElem-1 ; i++, src++ )
  {
    *dst++= *src;
    *dst++= (unsigned char)(((unsigned)(*src) + (unsigned)(*(src+1))) >> 1);
  }
  *dst= *(dst+1)= *src;
}


void
interpolateLines( unsigned char *src1, unsigned char *src2, unsigned char *dst,
                  int weight1, int weight2, unsigned long numElem )
{
  unsigned norm= weight1+weight2;
  for( unsigned long i= 0 ; i< numElem ; i++ )
    *dst++= (unsigned char)(((*src1++)*weight1 + (*src2++)*weight2) / norm);
}


int
main( int argc, char *argv[] )
{
  char	fileName[BUFFER_SIZE];
  char	cmdName[BUFFER_SIZE];
  char	dotY[3];
  unsigned long	i, j, k;
  
  CommandlineParser parser;
  
  // Luma resolution; defaults are for SONY SR1 AVCHD
  CoordinateVector yRes;
  yRes.vec.push_back( 1440 );
  yRes.vec.push_back( 1080 );
  CoordinateOption yResOption( yRes, "\tresolution of Luma plane\n",
                               "--y-res", NULL );
  parser.registerOption( &yResOption );
  
  // save fields or frames? (default fields)
  bool fields= true;
  BoolOption fieldOption( fields,
                          "\tchoose between saving fields and frames\n",
                          "--fields", "-fi", "--frames", "-fr" );
  parser.registerOption( &fieldOption );

  // whether to write only luminosity channel or colour too (default all)
  bool lumaOnly= false;
  BoolOption lumaOption( lumaOnly,
                          "\tchoose between extracting L or LUV channels\n",
                          "--luma-only", "-l", "--luma-chroma", "-luv" );
  parser.registerOption( &lumaOption );
  
  // whether or not to compress the MDA files (default true)
  bool compression= true;
  BoolOption comprOption( compression,
                          "\ttoggle gzip compression of results\n",
                          "--compression", "-c", "--no-compression", "-nc" );
  parser.registerOption( &comprOption );
  
  // integer offset for frame numbers
  int frame= 0;
  IntOption frameOption( frame,
                         "\toffset for frame numbers"
                         " (label for first extracted frame)\n",
                         "--frame-number-offset", NULL );
  parser.registerOption( &frameOption );
  
  // first frame to extract
  int startFrame= 0;
  IntOption startOption( startFrame,
                         "\tfirst frame to extract "
                         "(relative to beginning of file)\n",
                         "--first-frame", NULL );
  parser.registerOption( &startOption );
  
  // last frame to extract
  int endFrame= -1;
  IntOption endOption( endFrame,
                       "\tlast frame to extract "
                       "(relative to beginning of file)\n",
                       "--last-frame", NULL );
  parser.registerOption( &endOption );

  // interval between extracted frames
  int stepFrame= 1;
  IntOption stepOption( stepFrame,
                          "\textract every nth frame\n",
                          "--step-frame", NULL );
  parser.registerOption( &stepOption );
  
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
  if( !parser.parse( index, argc, argv ) || index< argc-1 || index> argc  )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }
  bool toFile= true;
  if( index== argc )
    toFile= false;
  
  if( yRes.vec.size()!= 2 )
  {
    cerr << argv[0] << ": luma must be a 2D array!\n";
    exit( 1 );
  }
  
  // output resolution depends on field vs. frame mode
  CoordinateVector uvRes;
  uvRes.vec.push_back( yRes.vec[0]/2 );
  uvRes.vec.push_back( yRes.vec[1]/2 );
  CoordinateVector yOutRes;
  yOutRes.vec.push_back( yRes.vec[0] );
  yOutRes.vec.push_back( fields ? yRes.vec[1]/2 : yRes.vec[1] );
  CoordinateVector uvOutRes;
  uvOutRes.vec.push_back( uvRes.vec[0] );
  uvOutRes.vec.push_back( fields ? uvRes.vec[1]/2 : uvRes.vec[1] );
  
  // I/O buffers
  unsigned char *yBuffer= new unsigned char[yRes.vec[0]*yRes.vec[1]];
  unsigned char *uvBuffer= new unsigned char[uvRes.vec[0]*uvRes.vec[1]*2];
  unsigned char *upsampledU= new unsigned char[yRes.vec[0]*yRes.vec[1]];
  unsigned char *upsampledV= new unsigned char[yRes.vec[0]*yRes.vec[1]];
  unsigned char *outScanline= new unsigned char[yRes.vec[0]*3];
  
  // skip initial frames until user specified startFrame
  for( i= 0 ; i< startFrame ; i++ )
  {
    if( cin.good() ) {
      //cin.read( (char *)yBuffer, yRes.vec[0]*yRes.vec[1] );
      //cin.read( (char *)uvBuffer, 2*uvRes.vec[0]*uvRes.vec[1] );
      // calling "ignore" is supposed to be faster than reading to dummy buffer,
      // but it's C++ implementation-dependent
      cin.ignore( yRes.vec[0]*yRes.vec[1] + 2*uvRes.vec[0]*uvRes.vec[1] );
    }
  }
  
  // extract as many frames as we can
  while( cin.good() )
  {
    // read luminance and chrominace planes for a full frame
    cin.read( (char *)yBuffer, yRes.vec[0]*yRes.vec[1] );
    cin.read( (char *)uvBuffer, uvRes.vec[0]*uvRes.vec[1]*2 );
    
    // write two fields or one frame
    for( i= 0 ; i<= (int)fields ; i++ )
    {
      // assemble the file name if we output to a file rather than
      // standard out
      if( toFile )
      {
        if( fields )
          if( lumaOnly )
            sprintf( fileName, "%s.%06d.%d.y.mda", argv[argc-1], frame, i );
          else
            sprintf( fileName, "%s.%06d.%d.yuv.mda", argv[argc-1], frame, i );
        else
          if( lumaOnly )
            sprintf( fileName, "%s.%06d.y.mda", argv[argc-1], frame, i );
          else
            sprintf( fileName, "%s.%06d.yuv.mda", argv[argc-1], frame, i );
        sprintf( cmdName, "gzip -9 %s", fileName );
      }
      
      
      // set up writer
      MDAWriter writer;
      if( toFile )
        writer.connect( fileName );
      else
        writer.connect( cout );
      if( !writer.writeHeader( yOutRes, lumaOnly ? 1 : 3, UByte ) )
      {
        cerr << argv[0] << ": Cannot write file header\n";
        exit( 1 );
      }
      
      // write field/frame
      if( lumaOnly )
        // one channel luminance
        for( j= 0 ; j< yOutRes.vec[1] ; j++ )
          writer.writeScanline( yBuffer + ((fields+1)*j+i)*yOutRes.vec[0] );
      else
      {
        // upsample u, v
        
        // first, upsample the lines where we have original data
        // these are lines 4*k and 4*k+1 of the final resolution
        for( j= 0 ; j< uvRes.vec[1]/2 ; j++ )
        {
          // even line of u
          upsampleLine( uvBuffer+(2*j)*uvRes.vec[0],
                        upsampledU+(4*j)*yRes.vec[0],
                        uvRes.vec[0] );
          // odd line of u
          upsampleLine( uvBuffer+(2*j+1)*uvRes.vec[0],
                        upsampledU+(4*j+1)*yRes.vec[0],
                        uvRes.vec[0] );
          // even line of v
          upsampleLine( uvBuffer+(2*j+uvRes.vec[1])*uvRes.vec[0],
                        upsampledV+(4*j)*yRes.vec[0],
                        uvRes.vec[0] );
          // odd line of v
          upsampleLine( uvBuffer+(2*j+1+uvRes.vec[1])*uvRes.vec[0],
                        upsampledV+(4*j+1)*yRes.vec[0],
                        uvRes.vec[0] );
        }
        
        // now, interpolate the lines where we do NOT have original data
        if( fields )
        {
          // in field mode, line 4k+2 is 1:1 blend of line 4k and 4k+4
          // and line 4k+3 is a 1:1 blend of 4k+1 and 4k+5
          for( j= 0 ; j< uvRes.vec[1]/2-1 ; j++ )
          {
            // line 4k+2 in u
            interpolateLines( upsampledU+(4*j)*yRes.vec[0],
                              upsampledU+(4*j+4)*yRes.vec[0],
                              upsampledU+(4*j+2)*yRes.vec[0],
                              1, 1, yRes.vec[0] );
            // line 4k+2 in v
            interpolateLines( upsampledV+(4*j)*yRes.vec[0],
                              upsampledV+(4*j+4)*yRes.vec[0],
                              upsampledV+(4*j+2)*yRes.vec[0],
                              1, 1, yRes.vec[0] );
            // line 4k+3 in u
            interpolateLines( upsampledU+(4*j+1)*yRes.vec[0],
                              upsampledU+(4*j+5)*yRes.vec[0],
                              upsampledU+(4*j+3)*yRes.vec[0],
                              1, 1, yRes.vec[0] );
            // line 4k+3 in v
            interpolateLines( upsampledV+(4*j+1)*yRes.vec[0],
                              upsampledV+(4*j+5)*yRes.vec[0],
                              upsampledV+(4*j+3)*yRes.vec[0],
                              1, 1, yRes.vec[0] );
          }
          // the last two lines are identical copies of the previous two
          memcpy( upsampledU+(yRes.vec[1]-2)*yRes.vec[0],
                  upsampledU+(yRes.vec[1]-4)*yRes.vec[0],
                  2*yRes.vec[0] );
          memcpy( upsampledV+(yRes.vec[1]-2)*yRes.vec[0],
                  upsampledV+(yRes.vec[1]-4)*yRes.vec[0],
                  2*yRes.vec[0] );
        }
        else
        {
          // in frame mode, line 4k+2 is 2:1 blend of line 4k+1 and 4k+4
          // and line 4k+3 is a 1:2 blend of 4k+1 and 4k+4
          for( j= 0 ; j< uvRes.vec[1]/2-1 ; j++ )
          {
            // line 4k+2 in u
            interpolateLines( upsampledU+(4*j+1)*yRes.vec[0],
                              upsampledU+(4*j+4)*yRes.vec[0],
                              upsampledU+(4*j+2)*yRes.vec[0],
                              2, 1, yRes.vec[0] );
            // line 4k+2 in v
            interpolateLines( upsampledV+(4*j+1)*yRes.vec[0],
                              upsampledV+(4*j+4)*yRes.vec[0],
                              upsampledV+(4*j+2)*yRes.vec[0],
                              2, 1, yRes.vec[0] );
            // line 4k+3 in u
            interpolateLines( upsampledU+(4*j+1)*yRes.vec[0],
                              upsampledU+(4*j+4)*yRes.vec[0],
                              upsampledU+(4*j+3)*yRes.vec[0],
                              1, 2, yRes.vec[0] );
            // line 4k+3 in v
            interpolateLines( upsampledV+(4*j+1)*yRes.vec[0],
                              upsampledV+(4*j+4)*yRes.vec[0],
                              upsampledV+(4*j+3)*yRes.vec[0],
                              1, 2, yRes.vec[0] );
          }
          // the last two lines are identical copies of the third-to-last
          memcpy( upsampledU+(yRes.vec[1]-2)*yRes.vec[0],
                  upsampledU+(yRes.vec[1]-3)*yRes.vec[0],
                  yRes.vec[0] );
          memcpy( upsampledU+(yRes.vec[1]-1)*yRes.vec[0],
                  upsampledU+(yRes.vec[1]-3)*yRes.vec[0],
                  yRes.vec[0] );
          memcpy( upsampledV+(yRes.vec[1]-2)*yRes.vec[0],
                  upsampledV+(yRes.vec[1]-3)*yRes.vec[0],
                  yRes.vec[0] );
          memcpy( upsampledV+(yRes.vec[1]-1)*yRes.vec[0],
                  upsampledV+(yRes.vec[1]-3)*yRes.vec[0],
                  yRes.vec[0] );
        }
        
        // assemble and write each scanline
        for( j= 0 ; j< yOutRes.vec[1] ; j++ )
        {
          // assemble the scanline
          unsigned long base= ((fields+1)*j+i)*yOutRes.vec[0];
          for( k= 0 ; k< yOutRes.vec[0] ; k++ )
          {
            outScanline[k*3]= yBuffer[base+k];
            outScanline[k*3+1]= upsampledU[base+k];
            outScanline[k*3+2]= upsampledV[base+k];
          }
          // and write it
          writer.writeScanline( outScanline );
        }
      }
      
      // disconnect writer, and compress file if desired by user
      writer.disconnect();
      if( compression && toFile )
        system( cmdName );
      if( verbose )
        if( toFile )
          cerr << argv[0] << ": Wrote " << fileName << endl;
        else
          if( fields )
            cerr << argv[0] << ": Wrote frame " << frame << ", field "
                 << i << endl;
          else
            cerr << argv[0] << ": Wrote frame " << frame << endl;
    }
    
    // keep track of #frames extracted;
    // abort if we have reached user-specified last frame
    if( startFrame++== endFrame )
      break;
    frame++;
    
    // skip over some frames
    if( stepFrame > 1 ) {
      for( int i=0; i<stepFrame-1; i++ ) {
        if( cin.good() ) {
          //cin.read( (char*)yBuffer, yRes.vec[0]*yRes.vec[1] );
          //cin.read( (char*)uvBuffer, 2*uvRes.vec[0]*uvRes.vec[1] );
          cin.ignore( yRes.vec[0]*yRes.vec[1] + 2*uvRes.vec[0]*uvRes.vec[1] );
        }
        if( startFrame++== endFrame )
          break;
      }
    }
    
    // check if there's more data to come (will set cin.good flag)
    cin.peek();
  }
  
  return 0;
}
