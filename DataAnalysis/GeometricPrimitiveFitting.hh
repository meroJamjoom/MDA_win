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

#ifndef GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_H
#define GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_H

/*! \file  GeometricPrimitiveFitting.hh
    \brief Fitting methods for geometric primitives.
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/GeometricPrimitives/Point.hh"
#include "MDA/GeometricPrimitives/Line.hh"

#include "SampleVector.hh"

namespace MDA {

    using namespace std;

    /** fit a point to a SampleVector of data-points */
    void fitToDataPoints( Point &p, const SampleVector &dataPoints );

    /** fit a line to a SampleVector of data-points */
    void fitToDataPoints( Line &l, const SampleVector &dataPoints );


} /* namespace */

#endif /* GEOMETRICPRIMITIVES_GEOMETRICPRIMITIVEFITTING_H */

