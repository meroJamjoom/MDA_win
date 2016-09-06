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

#ifndef CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_H
#define CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_H

/*! \file  ConnectedComponentProperties.hh
    \brief A class capsuling properties of the labels in a connected components image. For example covered pixel area, eulernumber, boundingbox etc.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Array/Array.hh"
#include "MDA/Filters/Filters.hh"
#include "MDA/Base/Range.hh"
#include "MDA/LinearAlgebra/Vector.hh"
#include "MDA/DataAnalysis/SampleVector.hh"

#include <math.h>
#include <vector>
#include <limits.h>


namespace MDA {

  enum Properties {
    NoProp= 0,
    NumLabels= 1,
    Area= 2,
    BoundingBox= 4,
    EulerNumber= 8,
    Segments= 16,
    FilledRegion= 32,
    Outline=64
  };

  /** \class ConnectedComponentProperties ConnectedComponentProperties.hh
      A class capsuling properties of the labels in a connected components image. Properties are:
      "NumLabels" ( The number of labels in the connected component image )
      "Area" ( The pixel area of all non-zero pixels of each component )
      "BoundingBox" ( The boundingbox of all each component )
      "EulerNumber" ( The eulernumber of each component )
      "Segments" ( The extracted components ) 
      "FilledRegion" ( The the filled region of each component )
      "Outline" ( The the outline ( 2D silhouette edge ) of each component )*/
  template<class T> 
  class ConnectedComponentProperties {

  public:

    /** default constructor. Takes the connected components in the first channel of the connected components image. */
    ConnectedComponentProperties(Array<T> *connComp);

    /** default destructor */
    ~ConnectedComponentProperties();

    /** computes the number of labels  */
    void computeNumLabels();

    /** computes the boundingboxes and / or area */
    void computeIndependentProps( bool computeBBox, bool computeArea);

    /** puts every boundingbox region in a segment array with a symetrical padding if in regions vector.
        if regions are not set, every boundingbox region is put into a segment array.*/
    void computeSegments( int pad, std::vector<unsigned int>* regions = NULL );

    /** compute eulernumbers of all valid segments  */
    void computeEulernumbers();

    /** computes the filled region of all valid segments if padded at least 1 pixel  */
    void computeFilledRegion();

    /** computes the outline (silhouette edge) of all valid segments in 2D if padded at least 1 pixel  */
    void computeOutline();

    /** computes samples outline (silhouette edge) of all 2D segments  in regions if padded at least 1 pixel 
	sampledistance is the distance on the outline between the distinct samples. */
    vector<SampleVector*>* extractOutlineSamples( unsigned int sampleDist, std::vector<unsigned int>* regions = NULL );

    /** reassembles the image from a the segments with indices in a vector*/
    Array<T>* reassemble( std::vector<unsigned int>* regions );

    /** returns if the specified property is computed. */
    bool getPropertyComputed( unsigned int property);

    /** returns the number of labels if computed, else return -1 */
    int getNumLabels();
  
    /** returns the boundingbox of the defined region*/
    LongRangeList getBoundingbox( unsigned int region );

    /** returns the area of the defined region*/
    unsigned long getArea( unsigned int region );

    /** returns the extracted segment of the defined region*/
    Array<T>* getSegment( unsigned int region );

    /** returns the symmetrical padding of segment of the defined region*/
    unsigned int getPadding( unsigned int region );

    /** returns the eulernumber of the defined region*/
    double getEulerNumber( unsigned int region );

    /** returns the number of extracted segments*/
    unsigned int getNumSegments();

    /** returns if the segment is valid or not*/
    bool isValid( unsigned int region );

  protected:

    /** array with the connected component labels */
    Array<T> *connComp;

    /** dimension vector of the array */
    CoordinateVector srcDim;

    /** total number of elements in first channel of array */
    unsigned long numEntries;

    /** segmented array */
    std::vector< Array<T>* > *segments;

    /** Entries in the segments */
    vector< unsigned long >* regionEntries;

    /** a map of all properties, if they have been computed or not */
    unsigned int computedProps;

    /** number of the connected components in the array */
    int numLabels;
   
    /** eulernumbers for connected components in the array */
    std::vector<double> *eulerNum;

    /** Boundingboxes for each connected component in the array */
    std::vector<LongRangeList*> *bboxes;

    /** Area for each connected component in the array */
    std::vector<unsigned long> *area;

    /** Symetrical padding for each segment */
    std::vector<unsigned int> *padding;

  };


} /* namespace */

#endif /* CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_H */

