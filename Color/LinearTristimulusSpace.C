// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef COLOR_LINEARTRISTIMULUSSPACE_C
#define COLOR_LINEARTRISTIMULUSSPACE_C

#include "LinearTristimulusSpace.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

//
// definitions and data local to this file
//
  
class NamedTristimulusSpace {
public:
  const char *name;
  double toXYZMatrix[9];
  double fromXYZMatrix[9];
};
    
NamedTristimulusSpace colorSpaces[]= {
  { "XYZ", // identity matrix - this is the default
    {1.0, 0.0, 0.0,
     0.0, 1.0, 0.0,
     0.0, 0.0, 1.0},
    {1.0, 0.0, 0.0,
     0.0, 1.0, 0.0,
     0.0, 0.0, 1.0}},
  { "sRGB", // linear sRGB space
    {0.4124, 0.3576, 0.1805,
     0.2126, 0.7152, 0.0722,
     0.0193, 0.1192, 0.9505},
    {3.2405, -1.5371, -0.4985,
     -0.9693, 1.8706, 0.0416,
     0.0556, -0.2040, 1.0572}},
  { "Rec709", // linear HDTV - same as sRGB
    {0.4124, 0.3576, 0.1805,
     0.2126, 0.7152, 0.0722,
     0.0193, 0.1192, 0.9505},
    {3.2405, -1.5371, -0.4985,
     -0.9693, 1.8706, 0.0416,
     0.0556, -0.2040, 1.0572}},
  { "xvYCC", // wide gamut HDTV - same as sRGB, but pussibly with negative RGB
    {0.4124, 0.3576, 0.1805,
     0.2126, 0.7152, 0.0722,
     0.0193, 0.1192, 0.9505},
    {3.2405, -1.5371, -0.4985,
     -0.9693, 1.8706, 0.0416,
     0.0556, -0.2040, 1.0572}},
  { "Adobe", // linear Adobe RGB
    {0.5767, 0.1856, 0.1882,
     0.2974, 0.6273, 0.0753,
     0.0270, 0.0707, 0.9911},
    {2.0414, -0.5649, -0.3447,
     -0.9693, 1.8760, 0.0416,
     0.0134, -0.1184, 1.0154}},
  { "Rec601", // linear standard definition TV
    {0.6069, 0.1735, 0.2003,
     0.2989, 0.5866, 0.1145,
     0.0000, 0.0661, 1.1162},
    {1.9100, -0.5325, -0.2882,
     -0.9847, 1.9992, -0.0283,
     0.0583, -0.1184, 0.8976}},
  { "HPE", // Hunt-Pointer-Estevez cone space
    {1.91019, -1.11214, 0.20195,
     0.37095, 0.62905, 0.0000,
     0.00000, 0.00000, 1.0000},
    {0.38971, 0.68898, -0.07868,
     -0.22981, 1.18340, 0.04641,
     0.00000, 0.00000, 1.00000}},
  { "BFD", // linear part of Bradford space
    {0.9870, -0.1471, 0.1600,
     0.4323, 0.5184, 0.0493,
     -0.0085, 0.0400, 0.9685},
    {0.8951, 0.2664, -0.1614,
     -0.7502, 1.7135, 0.0367,
     0.0389, -0.0685, 1.0296}},
  { "CAT02", // CIE CAT02 chromatic adaptation space
    {1.0961, 0.2789, 0.1827,
     0.4544, 0.4735, 0.0721,
     -0.0096, -0.0057, 1.0153},
    {0.7328, 0.4296, -0.1624
     -0.7036, 1.6975, 0.0061,
     0.0030, 0.0136, 0.9834}},
  { "", // end of list
  }
};


//
// LinearTristimulusSpace members
//

/** constructor for named space - doubles as default constructor */
LinearTristimulusSpace::LinearTristimulusSpace( const char *name )
  : toXYZMatrix( 3, 3 ), fromXYZMatrix( 3, 3 )
{
  unsigned i= 0;
  // find standard with matching name, or last default entry
  while( colorSpaces[i].name[0]!= '\0' )
    if( !strcmp( colorSpaces[i].name, name ) )
      break;
    else
      i++;
  if( colorSpaces[i].name[0]== '\0' )
  {
    i= 0;
    cerr << "Color space \"" << name << "\" not found - using \""
	 << colorSpaces[i].name << '\"' << endl;
  }
  
  // copy color matrices from table
  double *m= colorSpaces[i].toXYZMatrix;
  toXYZMatrix[0][0]= m[0]; toXYZMatrix[0][1]= m[1]; toXYZMatrix[0][2]= m[2];
  toXYZMatrix[1][0]= m[3]; toXYZMatrix[1][1]= m[4]; toXYZMatrix[1][2]= m[5];
  toXYZMatrix[2][0]= m[6]; toXYZMatrix[2][1]= m[7]; toXYZMatrix[2][2]= m[8];
  m= colorSpaces[i].fromXYZMatrix;
  fromXYZMatrix[0][0]= m[0];fromXYZMatrix[0][1]= m[1];fromXYZMatrix[0][2]= m[2];
  fromXYZMatrix[1][0]= m[3];fromXYZMatrix[1][1]= m[4];fromXYZMatrix[1][2]= m[5];
  fromXYZMatrix[2][0]= m[6];fromXYZMatrix[2][1]= m[7];fromXYZMatrix[2][2]= m[8];
}


/** convert from this space to XYZ */
void
LinearTristimulusSpace::toXYZ( const Vector &tristimulus, Vector &XYZ ) const
{
  toXYZMatrix.rightMultiply( tristimulus, XYZ );
}   


/** convert from XYZ to this space */
void
LinearTristimulusSpace::fromXYZ( const Vector &XYZ, Vector &tristimulus ) const
{
  fromXYZMatrix.rightMultiply( XYZ, tristimulus );
}


/** return a list of the names of known standard gammas */
void
LinearTristimulusSpace::getStandardNames( list<const char *> &names )
{
  unsigned i= 0;
  while( colorSpaces[i].name[0]!= '\0' )
    names.push_back( colorSpaces[i++].name );
}



} /* namespace */

#endif /* COLOR_LINEARTRISTIMULUSSPACE_C */

