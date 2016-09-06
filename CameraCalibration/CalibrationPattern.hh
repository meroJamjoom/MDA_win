// ==========================================================================
// $Id: CalibrationPattern.hh 702 2010-04-15 23:39:58Z atcheson $
// abstract baseclass for (2D or 3D) calibration patterns
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

#ifndef CAMERACALIBRATION_CALIBRATIONPATTERN_H
#define CAMERACALIBRATION_CALIBRATIONPATTERN_H

/*! \file  CalibrationPattern.hh
    \brief abstract baseclass for (2D or 3D) calibration patterns
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif


#include <MDA/LinearAlgebra/LinAlg.hh>

#include "PointCorrespondence.hh"


namespace MDA {

  using namespace std;
  
  /** \class CalibrationPattern CalibrationPattern.hh
      abstract baseclass for (2D or 3D) calibration patterns */
  class CalibrationPattern {

  public:
    
    /** default constructor */
    CalibrationPattern()
    {
      correspondences.clear();
    }
    
    /** read configuration from MetaData (pure virtual) */
    virtual bool configure( MetaData &md, const char *path )= 0;
    
    /** write configure to MetaData (pure virtual ) */
    virtual void writeConfig( MetaData &md, const char *path )= 0;
    
    /** run the detection algorithm for the grid (pure virtual) */
    virtual bool detect()= 0;
    
    /** return whatever corresponences have been computed */
    inline const vector<PointCorrespondence> &getCorrespondences() const
    {
      return correspondences;
    }

    /** write point correspondences to XML file 'dest' (pure virtual) */
    virtual void writeCorrespondences( MetaData& md, const char* path, const char* dest )= 0;

  protected:
    
    /** list of point correspondences */
    vector<PointCorrespondence> correspondences;
    
  };


} /* namespace */

#endif /* CAMERACALIBRATION_CALIBRATIONPATTERN_H */

