// ==========================================================================
// $Id: PaperSize.C 661 2010-03-22 01:19:14Z heidrich $
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

#ifndef CAMERACALIB_PAPERSIZE_C
#define CAMERACALIB_PAPERSIZE_C

#include "PaperSize.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  //  using namespace std;

  /** repository of paper dimensions */
  const double PaperSize::dimensions[][2]= {
    { 11.69, 16.54 },	// A3
    {  8.27, 11.69 },	// A4
    {  5.83,  8.27 },	// A5
    {  8.50, 11.00 },	// Letter
    {  8/50, 14.00 }	// Legal
  };


} /* namespace */

#endif /* CAMERACALIB_PAPERSIZE_C */

