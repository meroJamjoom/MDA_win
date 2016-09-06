// ==========================================================================
// $Id: ExpressionHelperFunctions.hh 319 2009-05-27 21:17:01Z heidrich $
// simple helper functions for the expresison parser
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

#ifndef EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_H
#define EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_H

/*! \file  ExpressionHelperFunctions.hh
    \brief simple helper functions for the expresison parser
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif


#if defined (_WIN32) || defined (_WIN64)
typedef double(*UNIVARIATE1) (double);
typedef double(*BIVARIATE1) (double, double);
#endif
/* callbacks for computing minumum and maximum of two numbers */
extern double minOp( double a, double b );
extern double maxOp( double a, double b );

/* callbacks for "missing" arithmetic functions */
extern double fracOp( double a );
extern double sgnOp( double a );

/* helper functions that are used as callbacks for built-in operators */
extern double plusOp( double a, double b );
extern double minusOp( double a, double b );
extern double multOp( double a, double b );
extern double divOp( double a, double b );
extern double uminusOp( double a );
extern double lessOp( double a, double b );
extern double greaterOp( double a, double b );
extern double eqOp( double a, double b );
extern double neqOp( double a, double b );
extern double leqOp( double a, double b );
extern double geqOp( double a, double b );
extern double andOp( double a, double b );
extern double orOp( double a, double b );
extern double notOp( double a );

#endif /* EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_H */

