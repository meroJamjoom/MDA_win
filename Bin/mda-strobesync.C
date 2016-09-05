// ==========================================================================
// $Id: mda-strobesync.C $
// read a sequence of MDA fields captured with the new strobes
// and output a sequence of MDA full frames that start with the first flash
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: bradleyd (Derek Bradley)
// Email:   bradleyd@cs.ubc.ca
// ==========================================================================

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <memory.h>
#include <math.h>

#include "MDA/Array/MDAFileIO.hh"


using namespace MDA;
using namespace std;

#define USAGE_TEXT "[<options>]\n"

// I could not make it work robustly without some thresholds
#define DARK_VALUE_OFFSET 5      // anything above (darkValue+DARK_VALUE_OFFSET) is scene pixels.
#define CONSISTENT_SCANLINES 20  // the number of consecutive scanlines that need to agree about
                                 // being dark or light.
#define SAFETY_SCANLINES 30      // how many scanlines to backup in order to safely be in a dark region

// parameters
int startField = -1, startLine = -1; 

// other temporary variables
int curField=-1;
int darkValue=21;
float lowthresh;
int avgIntensity; //average intensity of non-zero scanlines
bool even=true;   // are we currently processing an even or odd field
double *F, *B;     // current frame and a buffer (for output)
double *Fu, *Fv, *Bu, *Bv; // for U and V components
double *Wu, *Wv; // weights
unsigned long fieldW, fieldH;
int numChannels;
unsigned char* curY=NULL; // input fields
unsigned char* prevY=NULL;
unsigned char* curU=NULL;
unsigned char* prevU=NULL;
unsigned char* curV=NULL;
unsigned char* prevV=NULL;
unsigned char* tmpBuf=NULL;
bool usePrev = false;
float GAMMA, GTHRESH, SLOPE, GAIN, BIAS;
bool CLAMP;

enum State {LOW1=0, HIGH, LOW2};

double gammaFunc( double value, double gamma, double threshold, double slope,
	   double gain, double bias, bool clamp )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( value< 0.0 )
    if( clamp )
      value= 0.0;
    else
    {
      inversion= -1.0;
      value= -value;
    }
  // clamp against 1
  if( value> 1.0 && clamp )
    value= 1.0;
  
  if( value<= threshold )
    return inversion * value * slope;
  else
    return inversion * (gain * pow( value, 1.0/gamma ) + bias);
}

double gammaFuncInv( double value, double gamma, double threshold, double slope,
	      double gain, double bias, bool clamp )
{
  // clamp or reflect against zero
  double inversion= 1.0;
  if( value< 0.0 )
    if( clamp )
      value= 0.0;
    else
    {
      inversion= -1.0;
      value= -value;
    }
  // clamp against 1
  if( value> 1.0 && clamp )
    value= 1.0;
  
  if( value<= slope*threshold )
    return inversion * value / slope;
  else
    return inversion * pow( (value - bias) / gain, gamma );
}

// gamma to linear
void g2l(double* f, unsigned int _w, unsigned int _h)
{
  for (int i=0; i<_w; i++)
    {
      for (int j=0; j<_h; j++)
	{
	  double val = f[j*_w+i]/255.0;
	  val = gammaFuncInv(val, GAMMA, GTHRESH, SLOPE, GAIN, BIAS, CLAMP);
	  f[j*_w+i] = val*255.0;
	}
    }
}

// linear to gamma
void l2g(double* f, unsigned int _w, unsigned int _h)
{
  for (int i=0; i<_w; i++)
    {
      for (int j=0; j<_h; j++)
	{
	  double val = f[j*_w+i]/255.0;
	  val = gammaFunc(val, GAMMA, GTHRESH, SLOPE, GAIN, BIAS, CLAMP);
	  f[j*_w+i] = val*255.0;
	}
    }
}


bool loadArray( unsigned long &width, unsigned long &height, int &numChannels,
	        bool firstTime= false )
{
  MDAReader reader;
  reader.connect( cin );
  if( !reader.readHeader() )
    return false;
  
  // first time setup
  CoordinateVector dim;
  if( firstTime )
  {
    dim= reader.getDim();
    if( dim.vec.size()== 2 )
    {
      width= dim.vec[0];
      height= dim.vec[1];
    }
    numChannels= reader.getNumChannels();
  }
  
  // run a consistency check on the MDA
  dim= reader.getDim();
  if( dim.vec.size()!= 2 )
  {
    cerr << "Expecting 2D array!\n";
    return false;
  }
  if( dim.vec[0]!= width || dim.vec[1]!= height )
  {
    cerr << "Dimensions do not match previous frames\n";
    return false;
  }
  if( reader.getNumChannels()!= numChannels )
  {
    cerr << "Number of channels does not match previous frames\n";
    return false;
  }
  
  // allocate buffer if not yet present
  if( curY== NULL )
    {
      curY= new unsigned char[width*height];
      prevY = new unsigned char[width*height];
    }
  if (numChannels > 1)
    {
      if (curU==NULL)
	{
	  curU= new unsigned char[width*height];
	  prevU = new unsigned char[width*height];
	}
      if (curV==NULL)
	{
	  curV= new unsigned char[width*height];
	  prevV = new unsigned char[width*height];
	}
    }
  if (tmpBuf==NULL)
    {
      tmpBuf = new unsigned char[width*height*numChannels];
    }

  unsigned long numScanlines= reader.getNumScanlinesLeft();

  // move current image to previous
  for( unsigned long i= 0 ; i< numScanlines ; i++ )
    {
      memcpy( prevY+i*width, curY+i*width, width );
      if (curU)
	memcpy( prevU+i*width, curU+i*width, width );
      if (curV)
	memcpy( prevV+i*width, curV+i*width, width );
    }

  // read the MDA contents into the temp buffer

  for( unsigned long i= 0 ; i< numScanlines ; i++ )
    memcpy( tmpBuf+i*numChannels*width,
	    reader.readScanline(), numChannels*width );
  
  reader.disconnect();

  if (numChannels==1)
    memcpy(curY, tmpBuf, width*height);
  else
    {
      // split up Y, U and V
      for (unsigned long i=0; i<width*height; i++)
	{
	  curY[i] = tmpBuf[i*3];
	  curU[i] = tmpBuf[i*3+1];
	  curV[i] = tmpBuf[i*3+2];
	}
    }
  
  curField++;

  return true;
}

void determineDarkValue()
{
  if( !loadArray( fieldW, fieldH, numChannels, true ) )
  {
    cerr << "mda-strobeSync: Cannot load first frame\n";
    exit( 1 );
  }

  long sum=0;
  for (int y=0; y<fieldH; y++)
    {
      for (int x=0; x<fieldW; x++)
	{
	  sum += curY[y*fieldW+x];
	}
    }

  darkValue = (int)round(sum/(float)(fieldW*fieldH));
  fprintf(stderr, "dark value is %d\n", darkValue);
}

void detectStartScanline()
{
  if( !loadArray( fieldW, fieldH, numChannels) )
  {
    cerr << "mda-strobeSync: Cannot load frame\n";
    exit( 1 );
  }
    
  State state = LOW1;
  int consistent;

  while (state == LOW1)
    {
      consistent = 0;
      // go through the scanlines of the current image
      for (int y=0; y<fieldH; y++)
	{
	  if (state==HIGH) continue;

	  // find the max value
	  int max=0;
	  for (int x=0; x<fieldW; x++)
	    {
	      if (curY[y*fieldW+x] > max) max = curY[y*fieldW+x] ;
	    }
	  
	  // see if we have exited the low state
	  if (max > lowthresh)
	    {
	      consistent++;
	      if (consistent >= CONSISTENT_SCANLINES)
		{
		  state = HIGH;
		  if (y>=CONSISTENT_SCANLINES)
		    fprintf(stderr, "entered high state in field %d line %d\n", curField, y-CONSISTENT_SCANLINES);
		  else
		    fprintf(stderr, "entered high state in field %d line %d\n", curField-1, fieldH + (y-CONSISTENT_SCANLINES));
		}
	    }
	  else
	    consistent = 0;
	}
      
      // load the next image
      loadArray(fieldW, fieldH, numChannels);
    }

  while (state == HIGH)
    {
      consistent = 0;
      // go through the scanlines of the current image
      for (int y=0; y<fieldH; y++)
	{
	  if (state==LOW2) continue;

	  // find the max value
	  int max=0;
	  for (int x=0; x<fieldW; x++)
	    {
	      if (curY[y*fieldW+x] > max) max = curY[y*fieldW+x] ;
	    }
	  
	  // see if we have exited the high state
	  if (max < lowthresh)
	    {      
	      consistent++;
	      if (consistent >= CONSISTENT_SCANLINES)
		{
		  state = LOW2;
		  if (y>=CONSISTENT_SCANLINES)
		    fprintf(stderr, "entered low2 state in field %d line %d\n", curField, y-CONSISTENT_SCANLINES);
		  else
		    fprintf(stderr, "entered low2 state in field %d line %d\n", curField-1, fieldH + (y-CONSISTENT_SCANLINES));
		}
	    }
	  else
	    consistent = 0;
	}
      
      // load the next image
      loadArray(fieldW, fieldH, numChannels);
    }
  
  while (state == LOW2)
    {
      consistent = 0;
      // go through the scanlines of the current image
      for (int y=0; y<fieldH; y++)
	{
	  // find the max value
	  int max=0;
	  for (int x=0; x<fieldW; x++)
	    {
	      if (curY[y*fieldW+x] > max) max = curY[y*fieldW+x] ;
	    }
	  
	  // see if we have exited the low state
	  if (max > lowthresh)
	    {      
	      consistent++;
	      
	      if (consistent >= CONSISTENT_SCANLINES)
		{
		  // we found it.  Subtract pixels to make sure we're in the black region
		  startField = curField;
		  startLine = y - CONSISTENT_SCANLINES - SAFETY_SCANLINES;
		  
		  if (startLine < 0)
		    {
		      usePrev=true;
		      startField--;
		      startLine = fieldH + startLine;
		    }
		  
		  fprintf(stderr, "found start field: %d and start line: %d\n", startField, startLine);
		  return;
		}
	    }
	  else
	    consistent = 0;
	}
      // load the next image
      loadArray(fieldW, fieldH, numChannels);
    }

  fprintf(stderr, "Could not detect start line\n");
  delete [] curY;
  delete [] prevY;
  exit(1);
}

// copy the field from scanline S to scanline T into the buffer B
// (take into account if this is an even or odd field)
void copyToBuf(double* Buf, double* field, int S, int T, unsigned int w)
{
  int fs = S*2;
  if (!even)
      fs++;

  for (int i=S, j=fs; i<=T; i++, j+=2)
    {
      for (int x=0; x<w; x++)
	{
	  Buf[j*w+x] = field[i*w+x];
	}
    }
}

double CubicInterpolate(
   double y0,double y1,
   double y2,double y3,
   double mu)
{
   double a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

void interpolate(double* img, double thresh, unsigned int w, unsigned int h2)
{
  for (int y=0; y<h2; y++)
    {
      for (int x=0; x<w; x++)
	{
	  double val = img[y*w+x];
	  
	  if(val <= thresh)
	    {
	      // interpolate
	      if (y<=2) 
		{
		  img[y*w+x] = img[(y+1)*w+x];
		}
	      else if (y>=h2-3) 
		{
		  img[y*w+x] = img[(y-1)*w+x];
		}
	      else 
		{
		  // linear
		  //img[y*w+x] = (img[(y-1)*w+x] + img[(y+1)*w+x])/ 2.0;

		  // cubic
		  img[y*w+x] = CubicInterpolate(img[(y-3)*w+x], img[(y-1)*w+x], img[(y+1)*w+x], img[(y+3)*w+x], 0.5);
		}
	    }
	}
    }

//   fprintf(stderr, "leaving interp\n");
}

void addBtoF(double* _B, double* _F, unsigned int w, unsigned int h2)
{
  for (int i=0; i<w; i++)
    {
      for (int j=0; j<h2; j++)
	{
	  _F[j*w+i] += _B[j*w+i];
	}
    }
}

void multiply(double* _Bu, double *_B, double* _W, unsigned int w, unsigned int h2)
{
    for (int i=0; i<w; i++)
    {
      for (int j=0; j<h2; j++)
	{
	  _Bu[j*w+i] *= _B[j*w+i];
	  _W[j*w+i] += _B[j*w+i];
	}
    }
}

void outputFrames()
{
  // my code uses w, h and h2
  int w = fieldW;
  int h = fieldH;
  int h2 = fieldH*2;

  // is the first field even or odd?
  even = (startField %2 == 0);
  
  // need to have the field in double form for gamma correction
  double * dfield=NULL;
  dfield = new double[w*h];

  double * dfieldu=NULL;
  double * dfieldv=NULL;
  if (numChannels == 3)
    {
      dfieldu = new double[w*h];
      dfieldv = new double[w*h];
    }
  unsigned char* final;
  final = new unsigned char[w*h2*numChannels];

  bool moreToCome=true;
  unsigned char* img, *imgu, *imgv;

  // which field do I start with?
  if (usePrev)
    {
      img = prevY;
      imgu = prevU;
      imgv = prevV;
    }
  else
    {
      img = curY;
      imgu = curU;
      imgv = curV;
    }

  // prepare writer
  MDAWriter writer;
  CoordinateVector dimOut;
  dimOut.vec.push_back( w );
  dimOut.vec.push_back( h2 );

  while (moreToCome)
    {
      // first field ////////////////////////////////
      // make double copy
      for (int i=0; i<w; i++)
	{
	  for (int j=0; j<h; j++)
	    {
	      dfield[j*w+i]=(double)img[j*w+i];
	      if (numChannels == 3)
		{
		  dfieldu[j*w+i]=(double)imgu[j*w+i];
		  dfieldv[j*w+i]=(double)imgv[j*w+i];
		}
	    }
	}

      // gamma to linear
      g2l(dfield, w, h);

      // zero out F's and B's and weights
      memset(F, 0, w*h2*sizeof(double));
      memset(B, 0, w*h2*sizeof(double));
      if (numChannels == 3)
	{
	  memset(Fu, 0, w*h2*sizeof(double));
	  memset(Bu, 0, w*h2*sizeof(double));
	  memset(Fv, 0, w*h2*sizeof(double));
	  memset(Bv, 0, w*h2*sizeof(double));
	  memset(Wu, 0, w*h2*sizeof(double));
	  memset(Wv, 0, w*h2*sizeof(double));
	}

      // Y //
      // copy the current field to B from the startline down
      copyToBuf(B, dfield, startLine, h-1, w);

      // interpolate B
      interpolate(B, 0.0, w, h2);

      // add B to F
      addBtoF(B, F, w, h2);

      if (numChannels == 3)
	{
	  // U //
	  // copy the current field to Bu from the startline down
	  copyToBuf(Bu, dfieldu, startLine, h-1, w);
	  
	  // interpolate Bu
	  interpolate(Bu, 0.0, w, h2);
	  
	  // multiply Bu by B (Y value as weights), and add B to Wu
	  multiply(Bu, B, Wu, w, h2);
	  
	  // add Bu to Fu
	  addBtoF(Bu, Fu, w, h2);
	  
	  // V //
	  // copy the current field to Bv from the startline down
	  copyToBuf(Bv, dfieldv, startLine, h-1, w);
	  
	  // interpolate Bv
	  interpolate(Bv, 0.0, w, h2);
	  
	  // multiply Bv by B (Y value as weights), and add B to Wv
	  multiply(Bv, B, Wv, w, h2);
	  
	  // add Bv to Fv
	  addBtoF(Bv, Fv, w, h2);
	}

      // next field (2nd) /////////////////////////////
      even = !even;

      // get next field
      if (!usePrev)
	{
	  moreToCome = loadArray(fieldW, fieldH, numChannels);
	  if (!moreToCome) break;
	}
      img = curY;
      imgu = curU;
      imgv = curV;
      usePrev = false;

      // make double copy
      for (int i=0; i<w; i++)
	{
	  for (int j=0; j<h; j++)
	    {
	      dfield[j*w+i]=(double)img[j*w+i];

	      if (numChannels == 3)
		{
		  dfieldu[j*w+i]=(double)imgu[j*w+i];
		  dfieldv[j*w+i]=(double)imgv[j*w+i];
		}
	    }
	}

      // subtract dark value from this one (clamping at 0)
      for (int i=0; i<w; i++)
	{
	  for (int j=0; j<h; j++)
	    {
	      double val=dfield[j*w+i];
	      val -= (double)darkValue;
	      if (val<0) dfield[j*w+i] = 0;
	      else dfield[j*w+i] = val;
	    }
	}

      // gamma to linear
      g2l(dfield, w, h);

      // zero out B's
      memset(B, 0, w*h2*sizeof(double));
      if (numChannels == 3)
	{
	  memset(Bu, 0, w*h2*sizeof(double));
	  memset(Bv, 0, w*h2*sizeof(double));
	}

      // Y //
      // copy this whole field to B
      copyToBuf(B, dfield, 0, h-1, w);
      
      // interpolate B
      interpolate(B, 0.0, w, h2);

      // add B to F
      addBtoF(B, F, w, h2);

      if (numChannels == 3)
	{
	  // U //
	  // copy this whole field to Bu
	  copyToBuf(Bu, dfieldu, 0, h-1, w);
	  
	  // interpolate Bu
	  interpolate(Bu, 0.0, w, h2);
	  
	  // multiply Bu by B (Y value as weights), and add B to Wu
	  multiply(Bu, B, Wu, w, h2);
	  
	  // add Bu to Fu
	  addBtoF(Bu, Fu, w, h2);
	  
	  // V //
	  // copy this whole field to Bv
	  copyToBuf(Bv, dfieldv, 0, h-1, w);
	  
	  // interpolate Bv
	  interpolate(Bv, 0.0, w, h2);
	  
	  // multiply Bv by B (Y value as weights), and add B to Wv
	  multiply(Bv, B, Wv, w, h2);
	  
	  // add Bv to Fv
	  addBtoF(Bv, Fv, w, h2);
	}

      // next field (3rd) /////////////////////////////
      even = !even;

      // get next field
      moreToCome = loadArray(fieldW, fieldH, numChannels);
      if (!moreToCome) break;
      img = curY;
      imgu = curU;
      imgv = curV;

      // make double copy
      for (int i=0; i<w; i++)
	{
	  for (int j=0; j<h; j++)
	    {
	      dfield[j*w+i]=(double)img[j*w+i];

	      if (numChannels == 3)
		{
		  dfieldu[j*w+i]=(double)imgu[j*w+i];
		  dfieldv[j*w+i]=(double)imgv[j*w+i];
		}
	    }
	}

      // gamma to linear
      g2l(dfield, w, h);

      // zero out B's
      memset(B, 0, w*h2*sizeof(double));
      if (numChannels == 3)
	{
	  memset(Bu, 0, w*h2*sizeof(double));
	  memset(Bv, 0, w*h2*sizeof(double));
	}


      // Y //
      // copy the current field to B from the startLine up
      copyToBuf(B, dfield, 0, startLine-1, w);

      // interpolate B
      interpolate(B, 0.0, w, h2);

      // add B to F
      addBtoF(B, F, w, h2);

      if (numChannels == 3)
	{  
	  // U //
	  // copy the current field to Bu from the startLine up
	  copyToBuf(Bu, dfieldu, 0, startLine-1, w);
	  
	  // interpolate Bu
	  interpolate(Bu, 0.0, w, h2);
	  
	  // multiply Bu by B (Y value as weights), and add B to Wu
	  multiply(Bu, B, Wu, w, h2);
	  
	  // add Bu to Fu
	  addBtoF(Bu, Fu, w, h2);
	  
	  // V //
	  // copy the current field to Bv from the startLine up
	  copyToBuf(Bv, dfieldv, 0, startLine-1, w);
	  
	  // interpolate Bv
	  interpolate(Bv, 0.0, w, h2);
	  
	  // multiply Bv by B (Y value as weights), and add B to Wv
	  multiply(Bv, B, Wv, w, h2);
	  
	  // add Bv to Fv
	  addBtoF(Bv, Fv, w, h2);
	}
      
      // FINALIZE 
       // linear to gamma for Y
      l2g(F, w, h2);

     // switch to unsigned char and clamp
      double val;
      for (int i=0; i<w; i++)
	{
	  for (int j=0; j<h2; j++)
	    {
	      int ind = j*w+i;
	      
	      // do Y
	      val=F[ind];
	      if (val<0) final[ind*numChannels] = 0;
	      else if (val>255) final[ind*numChannels] = 255;
	      else final[ind*numChannels] = (unsigned char)round(val);

	      if (numChannels == 3)
		{
		  // do U
		  if (Wu[ind] == 0 ) val = 127;
		  else val = Fu[ind] / Wu[ind];
		  final[ind*3+1] = (unsigned char)round(val);

		  // do V
		  if (Wv[ind] == 0 ) val = 127;
		  else val = Fv[ind] / Wv[ind];
		  final[ind*3+2] = (unsigned char)round(val);
		}
	    }
	}

      // write out this frame
      writer.connect( cout );
      if( !writer.writeHeader( dimOut, numChannels, UByte ) )
	{
	  cerr << "mda-strobesync: cannot write header\n";
	  exit( 1 );
	}
      for( int i= 0 ; i< h2 ; i++ )
	{
	  writer.writeScanline( final+i*w*numChannels );
	}
      if( !writer.disconnect() )
	{
	  cerr << "mda-strobesync: error during writing\n";
	  exit( 1 );
	}
      
//       // setup for next one
//       moreToCome = loadArray(fieldW, fieldH, numChannels);
//       img = curY;

    }

  delete [] final;
  delete [] dfield;
  if (numChannels == 3)
    {
      delete [] dfieldu;
      delete [] dfieldv;
    }
}

int main( int argc, char *argv[] )
{
  CommandlineParser parser;

  // setup options
  int standard= 1;
  list<const char *>standards;
  standards.push_back( "--srgb" );
  standards.push_back( "--bt709" );
  standards.push_back( "--xyYCC" );
  SelectionOption standardOpt( standard,
			       "\tselect parameters according to image/video"
			       " standard:\n"
			       "\tsRGB, ITU-R BT.709 (HDTV), or xyYCC (wide "
			       "gamut HDTV)\n",
			       standards );
  parser.registerOption( &standardOpt );

  IntOption startLineOption(startLine, "\tthe start scanline of the sequence (just before the first flash)\n", "--startline", "-sl");
  parser.registerOption(&startLineOption);

  IntOption startFieldOption(startField, "\twhich field contains the start scanline (first field is 0)\n", "--startfield", "-sf");
  parser.registerOption(&startFieldOption);

  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc )
  {
    parser.usage( argv[0], USAGE_TEXT );
    exit( 1 );
  }

  // apply standard curve if specified
  CLAMP=true;
  switch( standard )
  {
  case 0: // sRGB
    GTHRESH= 0.00304;
    SLOPE= 12.92;
    BIAS= -0.055;
    GAIN= 1.055;
    GAMMA= 2.4;
    break;
  case 2: // xvYCC: like BT.709 but without clamping
    CLAMP= false;
  case 1: // BT.709
    GTHRESH= 0.018;
    SLOPE= 4.5;
    BIAS= -0.099;
    GAIN= 1.099;
    GAMMA= 1.0 / 0.45;
    break;
  }

  // determine dark value - the average intensity when strobes are off
  // (just take average of the first field)
  determineDarkValue();

  // lowthresh = above this we have something interesting, below this is only dark value
  lowthresh = darkValue + DARK_VALUE_OFFSET; // need to account for noise

  // did the user specify the starting point?
  if (startLine < 0 || startField < 0)
    {
      // nope
      detectStartScanline();

      // a little debugging, output the start line and field as a filename
      char tmps[256];
      sprintf(tmps,"touch line_%d_field_%d\n", startLine, startField);
      system(tmps);
    }
  else
    {
      // yes, move to start field
      while (curField < startField)
	{
	  if( !loadArray( fieldW, fieldH, numChannels) )
	    {
	      cerr << "mda-strobeSync: Cannot load frame\n";
	      exit( 1 );
	    }
	}
    }

  // allocate my buffers
  F = new double[fieldW*fieldH*2];
  B = new double[fieldW*fieldH*2];
  Fu = new double[fieldW*fieldH*2];
  Bu = new double[fieldW*fieldH*2];
  Fv = new double[fieldW*fieldH*2];
  Bv = new double[fieldW*fieldH*2];
  Wu = new double[fieldW*fieldH*2];
  Wv = new double[fieldW*fieldH*2];


  // create and output the frames
  outputFrames();

  delete [] curY;
  delete [] prevY;
  delete [] curU;
  delete [] prevU;
  delete [] curV;
  delete [] prevV;
  delete [] tmpBuf;
  delete [] F;
  delete [] B;
  delete [] Fu;
  delete [] Bu;
  delete [] Fv;
  delete [] Bv;
  delete [] Wu;
  delete [] Wv;
  return 0;
}
