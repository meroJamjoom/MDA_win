// ==========================================================================
// $Id:$
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: Felix Heide
// Email:   heide@informatik.uni-siegen.de
// ==========================================================================

#ifndef CONNCOMPPROPS_SADDLEPOINTFINDER_H
#define CONNCOMPPROPS_SADDLEPOINTFINDER_H

/*! \file  SaddlePointFinder.hh
    \brief A 2D saddlepoint finder. Based on the matlab cameracalibration toolbox file cornerfinder_saddle_point.m 
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Array/Array.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "MDA/LinearAlgebra/Matrix.hh"
#include "MDA/LinearAlgebra/LAPACKBridge.hh"
#include "MDA/DataAnalysis/SampleVector.hh"
#include "MDA/Filters/Filter.hh"
#include "MDA/Filters/Linear1DFilter.hh"

namespace MDA {

  using namespace std;

  /** \class SaddlePointFinder SaddlePointFinder.hh
      A 2D saddlepoint finder. Based on the matlab cameracalibration toolbox file cornerfinder_saddle_point.m */
  template<class T> 
  class SaddlePointFinder {

  public:

    /** default constructor, mask size is wintx * winty, a  central zero zone of the size  wx2, wy2 can be added */
    SaddlePointFinder( Array<T> *image, int wintx = 5, int winty = 5, int wx2 = -1 , int wy2 = -1);

    /** default destructor */
    ~SaddlePointFinder();

      /** finds the saddlepoints in the image defined in the constructor from initial guess positions
      the results are returned in the guesses SampleVector. The boolean vector of the same size as SampleVector
      holds information about which points diverged */
    void findSaddles( SampleVector *guesses,  vector<bool>* divergedGuesses );

  protected:

    /// image and its dimensions the saddlepointfinder works on
    Array<T> *image;
    CoordinateVector imgDim;

    /// mask dimensions and actual mask
    int wintx; 
    int winty;
    double* mask;

    /// resolution of the subpixel saddlefinder
    double resolution;

    /// maximum number of iterations
    int MaxIter;

  };


} /* namespace */

#endif /* CONNCOMPPROPS_SADDLEPOINTFINDER_H */

