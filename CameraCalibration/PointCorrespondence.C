// ==========================================================================
// $Id: PointCorrespondence.C 755 2010-09-10 22:48:47Z bradleyd $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIB_POINTCORRESPONDENCE_C
#define CAMERACALIB_POINTCORRESPONDENCE_C

#include <sstream>

#include "PointCorrespondence.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


/** get a single PointCorrespondence tag */
bool
get( MetaData &m, const char *path, PointCorrespondence &c )
{
  string p( path );
  
  if( !warnCond( c.iPt.getSize()== 2 && c.wPt.getSize()== 3,
		"dimensional mismatch in PointCorrespondence" ) )
    return false;
  
  bool success= true;
  success&= get( m, (p+":u").c_str(), c.iPt[0] );
  success&= get( m, (p+":v").c_str(), c.iPt[1] );
  success&= get( m, (p+":x").c_str(), c.wPt[0] );
  success&= get( m, (p+":y").c_str(), c.wPt[1] );
  success&= get( m, (p+":z").c_str(), c.wPt[2] );
  
  return success;
}
  

/** set a single PointCorrespondence tag */
void
set( MetaData &m, const char *path, const PointCorrespondence &c )
{
  string p( path );
  
  bool realNums=true;
  for (int i=0; i<2; i++)
    {
      if ( isnan(c.iPt[i]) || isinf(c.iPt[i]) ||
	   isnan(c.wPt[i]) || isinf(c.wPt[i]) )
	realNums=false;
    }
  if (isnan(c.wPt[2]) || isinf(c.wPt[2]))
    realNums=false;

  if (!realNums)
    {
      fprintf(stderr, "Non-real number found in point correspondence when writing\n");
    }
  else
    {
      set( m, (p+":u").c_str(), c.iPt[0] );
      set( m, (p+":v").c_str(), c.iPt[1] );
      set( m, (p+":x").c_str(), c.wPt[0] );
      set( m, (p+":y").c_str(), c.wPt[1] );
      set( m, (p+":z").c_str(), c.wPt[2] );
    }
}
  
  
/** get a vector of PointCorrespondence tags that are all siblings
    under a common <corr> tag */
bool
get( MetaData &m, const char *path, vector<PointCorrespondence> &c )
{
  // read points as long as there are <i> and <w> tags
  bool found;
  unsigned i= 0;
  do
  {
    // check if there is another PointCorrespondence at the same level
    // as the previous ones
    ostringstream p;
    p << path << ".corr.c[" << i << ']';
    PointCorrespondence pc;
    found= get( m, p.str().c_str(), pc );
    // if so, add it to the vector
    if( found )
      c.push_back( pc );
    i++;
  }
  while( found );
  
  return true;
}
  
  
/** set a vector of PointCorrespondence tags that are all siblings
    under a common <corr> tag */
void
set( MetaData &m, const char *path, const vector<PointCorrespondence> &c )
{
  unsigned numPoints= c.size();
  // write all points as siblings
  for( unsigned i= 0 ; i< numPoints ; i++ )
  {
    ostringstream p;
    p << path << ".corr.c[" << i << ']';
    set( m, p.str().c_str(), c[i] );
  }
}


} /* namespace */

#endif /* CAMERACALIB_POINTCORRESPONDENCE_C */

