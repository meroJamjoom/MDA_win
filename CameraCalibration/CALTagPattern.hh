// ==========================================================================
// $Id: CALTagPattern.hh 754 2010-09-10 00:33:50Z bradleyd $
// a single, planar, tightly-packed CALTag grid
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

#ifndef CAMERACALIBRATION_CALTAGPATTERN_H
#define CAMERACALIBRATION_CALTAGPATTERN_H

/*! \file  CALTagPattern.hh
    \brief a single, planar, tightly-packed CALTag grid
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "CalibrationPattern.hh"
#include "PaperSize.hh"


///////////// From Felix's code //////////////////////////////
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <vector>
#include <string>

#include "MDA/GeometricTransform/GeometricTransform.hh"
#include "MDA/Filters/ConnectedComponentProperties.hh"
#include "MDA/DataAnalysis/SampleVector.hh"
#include "MDA/DataAnalysis/lloydQuadFinder.hh"
#include "MDA/DataAnalysis/ClusterModifier.hh"
#include "SaddlePointFinder.hh"
#include "MDA/LinearAlgebra/LinAlg.hh"
#include "MDA/LinearAlgebra/Matrix.hh"
#include "MDA/LinearAlgebra/JacobiRotation.hh"
#include "MDA/Resampling/PointSampler.hh"

//Debugging includes
#include "MDA/GeometricPrimitives/XFigWriter.hh"
#include "MDA/GeometricPrimitives/Line.hh"
#include "MDA/GeometricPrimitives/Point.hh"
//////////////////////////////////////////////////////////////


namespace MDA {

  /** \class CALTagPattern CALTagPattern.hh
      a single, planar, tightly-packed CALTag grid */
  
  class CALTagPattern: public CalibrationPattern {
    
  public:
    
    /** default constructor */
    inline CALTagPattern()
      : layout( 2 ), columns( 0 ), rows( 0 ), scale( 1.0 ), firstIndex( 0 ),
        bits( 16 ), cBits( 6 ), genPoly( 3 ), codes( NULL ), valid( NULL ),
        //[brad]
        unitCalibPoints( NULL ), bowtieIntersections(NULL), imageCalibPoints( NULL ),
        resCode( 4 ), resSquare( 8 ), dotThresh( -16 ), sampleDist( 0 ),
        angleThreshold( 90.0 ), successfulTracked( false ),
	unitSquare(NULL), codeBlobSamplePos(NULL), code(NULL)
    {}
    
    /** destructor */
    virtual ~CALTagPattern()
    {
      if( codes )
        delete [] codes;
      if( valid )
        delete [] valid;

      //[brad]
      if( unitCalibPoints )
        delete unitCalibPoints;
      if( imageCalibPoints )
        delete imageCalibPoints;

      if (unitSquare)
	delete unitSquare;
      if (bowtieIntersections)
	delete bowtieIntersections;
      if (codeBlobSamplePos)
	delete codeBlobSamplePos;
      if (code)
	delete code;
    }
    
    /** read configuration from MetaData */
    virtual bool configure( MetaData &md, const char *path );
    
    /** write configuration to MetaData */
    virtual void writeConfig( MetaData &md, const char *path );
    
    /** change pattern layout
        (0 or less for any parameter means to keep the current value)
        \param r: new #rows
        \param c: new #columns
        \param s: new scale factor
        \param l: new layout/pattern type (1: checker, 2: bowtie)
    */
    void changeLayout( unsigned r= 0, unsigned c= 0,
                       double s= -1.0, unsigned l= 0 );
    
    /** generate the ID list for the current geometric layout, either
        from scratch, or using existing IDs for the first markers */
    bool makeIDs( bool fromScratch= false );
    
    /** set IDs to specified list (either with or without CRC portion)  */
    bool setIDs( CoordinateVector &indices, bool makeChecksum );
    
    /** set start index for newly generated IDs */
    inline void setStart( unsigned index )
    {
      firstIndex= index;
    }
    
    /** run the detection algorithm for the grid */
    virtual bool detect();
    
    /** create a PostScript file with the calibration pattern */
    virtual bool makePattern( const char *fileName, const PaperSize &paper );

    ////////////////// Added by Brad ///////////////////////////////////////////////////////////////////////////////////////
    
    // after successfully running detect() the unitCalibPoints and imageCalibPoints properties will be filled
    void constructCorrespondences();

    // output the correspondences to xml file
    // this should probably be moved into the CalibrationPattern class
    virtual void writeCorrespondences( MetaData& md, const char* path, const char* dest );

    ////////////////// Taken from Felix's code /////////////////////////////////////////////////////////////////////////////
  protected:
    // comparison function for the indexed sorting
    static bool compareIndexPairs ( Vector i, Vector j);

    /** \class CodeBlobProps CALTag.hh
        Properties of the marker codes*/
    class CodeBlobProps {
    public:
      long ID;
      bool valid;
      long firstCorner;
      Vector orientation;
      SampleVector corners;
      
      /** default constructor */
      CodeBlobProps();

      /** default destructor */
      ~CodeBlobProps();
    };

  public:

    #define T double
    /** changes the input image to the parameterized image*/
    void setImage( Array<T> *array );

  protected:
    /** wellners gaussian adaptive 2D thresholding function, default mode is relative(mode = 0), mode = 1 is fixed.  */
    void adaptiveThreshold( Array<T> *inArray, double t = 15.0, int fsize = -1, int mode = 0 );

    /** rotates the code pattern by 90 degrees counterclockwise */
    void rotCode90();

    /** Checkerboard corner estimator using method described in:
        "Robust Recognition of Checkerboard Pattern for Deformable Surface Matching in Multiple Views"
	Parameter delta added: Delta is the maximum allowed percentage of differring opposite pixels in a valid circle.
    */
    bool validateSaddle( Array<T> *in, Vector saddle, int minR = 2, int maxR = 10, int window = 1, double alpha = 1.0, double kappa = 0.2, double delta = 0.3 );

    // Homography 2D functions:
    /** computes the 2d homography h from the points in to the points out. H needs to be a 3x3 matrix */
    void homography2D( const SampleVector* in, const SampleVector *out, Matrix* h  );

    /** applys a 2D homography h to a bunch of points */
    void homotrans( const SampleVector *points, Matrix* h, SampleVector *result  );

    /** normalizes a set of points ˆXi = T*Xi such that the centroid of the new points 
        has the coordinate (0, 0) and their average distance from the origin is sqrt(2). */
    void normalizePoints( SampleVector* points, Matrix* t );

    // Ransac functions:
    /** Function to determine if a set of 4 2D points give rise
      to a degeneracy in the calculation of a homography as needed by RANSAC.
      This involves testing whether any 3 of the 4 points the set is collinear. */ 
    bool isDegenerate( SampleVector* points, double threshold = 0.01 );

    /** Function to evaluate the symmetric transfer error of a homography with
        respect to a set of matched points as needed by RANSAC. */
    void homogdist2d(const SampleVector *pointsIn, const SampleVector *pointsOut, Matrix* H, double t, vector<unsigned long>* inliers );

    /** fits 2D homography using the RANSAC algorithm. Returns true if successful*/
    bool ransacFitHomography(const SampleVector *pointsIn, const SampleVector *pointsOut, Matrix* H, 
                             vector<unsigned long>* inliers, double t, unsigned long s = 4, unsigned long feedback = 0,
                             unsigned long maxDataTrials = 100, unsigned long maxTrials = 1000 );

    /** corrects distortion for a image coordinate vector */
    void correctDistortion( Vector *point, Vector *imgCenter, double halfDiagonal, double r);

    /** corrects distortion for a couple of points */
    void correctDistortion( SampleVector *points, Vector *imgCenter, double halfDiagonal, double r);

    /** find an r which makes the points collinear */
    double findUndistortR( const SampleVector *points, Vector *imgCenter, double halfDiagonal, double precision );

    /** find an r which maps the pointsIn to pointsOut by minimizing the distance between mapped points and pointsOut */
    double findInversionR( const SampleVector *pointsIn, const SampleVector *pointsOut, Vector *imgCenter, double halfDiagonal, double precision );

    /** rate collinearity */
    double computeDistanceToBestFittedLine( const SampleVector *points );

    /** compute summed distance between points */
    double computePointSetDifference( const SampleVector *pointsIn, const SampleVector *pointsOut  );

  protected:

    /** array with input image */
    Array<T> *array;

    /** dimensions of the input array */
    CoordinateVector arrayDim;

    /** number of entries in the array */
    unsigned long numEntries;

    /** resolution of interior binary code pattern, in arbitrary units. 
        For example if the code identifying each marker uses a 4x4 pattern 
        of dots then set resCode to 4. Currently only NxN patterns are supported. */
    //[brad] redundant: should always equal sqrt(bits)
    unsigned int resCode;

    /** resolution of the marker squares, in units. This must be greater than 
        resCode. For example, if the code identifing each marker uses a 4x4 pattern of 
        dots embedded in an 8x8 pixel square, then set resSquare to 8. 
        Currently only MxM patterns are supported. */
    //[brad] postscript file uses a hardcoded 8x8 marker size
    const unsigned int resSquare;
    
    /** width of border in units around the code dots in each square */
    double borderWidth;

    /** current code pattern */
    vector<bool>* code;

    /** approx max num holes (code dots) per square (note that this is higher than the
        actual number of holes possible in the marker patterns, because as part of the
        image processing we can sometimes create additional small holes */
    int dotThresh;

    /** Ouline sampling distance for speeding up things */
    unsigned int sampleDist;

    /** maximum allowed difference in degrees from median orientation of all markers */
    double angleThreshold;

    /** Derived parameters: */
    unsigned int numMarkers;

    /** number of bits in the resCode X resCode pattern that are dedicated to 
        encoding the numerical ID of the marker. The remaining bits are used for the checksum. */
    unsigned int idBits;

    /** precomputed unit square for code-blob sampling */
    SampleVector* unitSquare;

    /** A vector of the the bowtie intersection points in unit coordinates */
    SampleVector* bowtieIntersections;

    /** code-blob sampling positions */
    SampleVector* codeBlobSamplePos;

    /** A vector of the the calibration points in unit coordinates */
    SampleVector* unitCalibPoints;

    /** A vector of the corresponding points found in the image */
    SampleVector* imageCalibPoints;

    //[brad] adding member variables to replace parameters to track() function since detect() has none
    //[brad] need to add methods + interface for actually setting these
    int debug;
    string debugOutputFilename;
    double debugScaleFactor;
    bool successfulTracked;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  protected:
    
    /** start index for generating new marker IDs */
    unsigned firstIndex;
    
    /** layout type (1: checker, 2: bowties) */
    unsigned layout;
    
    /** # grid columns */
    unsigned columns;
    
    /** # grid rows */
    unsigned rows;
    
    /** grid scale */
    double scale;
    
    /** numeric marker IDs (full payload, including CRC where applicable) */
    CoordinateVector ids;
    
    //
    // parameters for the CRC codes
    //
    
    /** precompute the full code array and which codes are valid */
    void computeCodes();
    
    /** code rotation - hard coded to 4x4 configurations right now */
    unsigned rotate( unsigned code );
    
    /** array of codes (with CRC portion) */
    unsigned *codes;
    
    /** array of bits indicating which codes are valid */
    bool *valid;
    
    /** total number of bits */
    unsigned bits;
    
    /** number of CRC bits included in "bits" */
    unsigned cBits;
    
    /** generating polynomial for CRC */
    unsigned long genPoly;
    
  };


} /* namespace */

#endif /* CAMERACALIBRATION_CALTAGPATTERN_H */

