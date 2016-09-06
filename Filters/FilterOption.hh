// ==========================================================================
// $Id: FilterOption.hh 954 2012-06-05 20:37:36Z krim $
// commandline option for filter type
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

#ifndef FILTERS_FILTEROPTION_H
#define FILTERS_FILTEROPTION_H

/*! \file  FilterOption.hh
    \brief commandline option for filter type
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CommandlineParser.hh"


namespace MDA {

using namespace std;

/** filter type */
enum FilterType
{
  BilateralFiltering,
  BilateralFilteringMasked,
  BilateralGridFiltering,
  BilateralGridWeightedFiltering,
  BoxFiltering,
  ConnectedComponentFiltering,
  CornerFiltering,
  DilateFiltering,
  DistanceFiltering,
  EdgeFiltering,
  ErodeFiltering,
  EulerNumberFiltering,
  ExtremaFiltering,
  FastGaussFiltering,
  FirstDerivativeFiltering,
  FloodfillFiltering,
  GaussFiltering,
  HatFiltering,
  HitOrMissFiltering,
  LoGFiltering,
  MedianFiltering,
  MedianMaskFiltering,
  PruneFiltering,
  SecondDerivativeFiltering,
  SeparableFiltering,
  SobelFiltering,
  ThinningFiltering,
  UnsharpMaskFiltering,
  Thinning3DFiltering,
  UndefinedFiltering // this one should always be last
};


/** names of filters in the order of the enum above */
extern const char *filterNames[];

/** write filter name to ostream */
ostream &operator<<( ostream &os, FilterType filter );
  
/** read filter name from istream */
istream &operator>>( istream &is, FilterType &filter );
  


  /** \class FilterOption FilterOption.hh
      commandline option for filter type */
  class FilterOption: public CommandlineOption {
    
  public:

    /** default constructor */
    FilterOption( FilterType &type,
		  const char *msg= "\tFilter\n",
		  const char *lTxt= "--filter", const char *sTxt= "-f"  )
      : CommandlineOption( msg, lTxt, sTxt ), filterType( type )
    {}

    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );

    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:

    /** reference to the filter type variable */
    FilterType &filterType;
  
  };


} /* namespace */



#endif /* FILTERS_FILTEROPTION_H */

