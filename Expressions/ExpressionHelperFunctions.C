// ==========================================================================
// $Id: ExpressionHelperFunctions.C 319 2009-05-27 21:17:01Z heidrich $
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

#ifndef EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_C
#define EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_C

#include <math.h>

#include "ExpressionHelperFunctions.hh"

/* callbacks for computing minumum and maximum of two numbers */
double minOp( double a, double b ) {return a< b ? a : b;}
double maxOp( double a, double b ) {return a> b ? a : b;}

/* callbacks for "missing" arithmetic functions */
double fracOp( double a ) {return a-floor( a );}
double sgnOp( double a ) {return a<0.0 ? -1.0 : (a>0.0 ? 1.0 : 0.0 );}
  
/* helper functions that are used as callbacks for built-in operators */
double plusOp( double a, double b ) {return a+b;}
double minusOp( double a, double b ) {return a-b;}
double multOp( double a, double b ) {return a*b;}
double divOp( double a, double b ) {return a/b;}
double uminusOp( double a ) {return -a;}
double lessOp( double a, double b ) {return a<b ? 1.0 : 0.0;}
double greaterOp( double a, double b ) {return a>b ? 1.0 : 0.0;}
double eqOp( double a, double b ) {return a==b ? 1.0 : 0.0;}
double neqOp( double a, double b ) {return a!=b ? 1.0 : 0.0;}
double leqOp( double a, double b ) {return a<=b ? 1.0 : 0.0;}
double geqOp( double a, double b ) {return a>=b ? 1.0 : 0.0;}
double andOp( double a, double b ) {return (a> 0.0 && b> 0.0) ? 1.0 : 0.0;}
double orOp( double a, double b ) {return (a> 0.0 || b> 0.0) ? 1.0 : 0.0;}
double notOp( double a ) {return (a<= 0.0) ? 1.0 : 0.0;}

#endif /* EXPRESSIONS_EXPRESSIONHELPERFUNCTIONS_C */

