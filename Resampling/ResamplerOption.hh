// ==========================================================================
// $Id: ResamplerOption.hh 200 2008-05-12 23:46:28Z heidrich $
// meta include file for resampling interface
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

#ifndef RESAMPLING_RESAMPLEOPTION_H
#define RESAMPLING_RESAMPLEOPTION_H

/*! \file  Resampling.hh
    \brief meta include file for resampling interface
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "Resampler.hh"
#include "MDA/Base/CommandlineParser.hh"


namespace MDA {

using namespace std;

/** supported resamplers */
enum ResampleMode
{
  NearestNeighborResampling,
  BoxResampling,
  LinearResampling,
  CubicResampling,
  GaussResampling,
  MinResampling,
  MaxResampling,
  UndefinedResampling // this one should always be last
};

/** names of resamplers in the order of the enum above */
extern const char *resamplerNames[];

/** write resampler name to ostream */
ostream &operator<<( ostream &os, ResampleMode resampler );

/** read resampler name from istream */
istream &operator>>( istream &is, ResampleMode &resampler );


  /** \class ResamplerOption Resampling.hh
      parser for resampling filter option */
  class ResamplerOption: public CommandlineOption {

  public:
    
    /** constructor from reference to ResampleMode object */
    ResamplerOption( ResampleMode &mode,
		     const char *msg= "\tResampling filter\n",
		     const char *lTxt= "--resampling-filter",
		     const char *sTxt= "-rf"  )
      : CommandlineOption( msg, lTxt, sTxt ), resampleMode( mode )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );

    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the resampling mode variable */
    ResampleMode &resampleMode;
    
  };


} /* namespace */



#endif /* RESAMPLING_RESAMPLEOPTION_H */

