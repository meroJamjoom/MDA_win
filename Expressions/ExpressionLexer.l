/* ==========================================================================
 * $Id: ExpressionLexer.l 243 2008-11-16 01:20:52Z heidrich $
 * a lexical analyzer for arithmetic expressions
 * ==========================================================================
 * License: Internal use at UBC only! External use is a copyright violation!
 * ==========================================================================
 * (C)opyright:
 *
 * 1995-2007, Wolfgang Heidrich, UBC
 *
 * Creator: heidrich (Wolfgang Heidrich)
 * Email:   heidrich@cs.ubc.ca
 * ==========================================================================
 */

%option prefix="EXPR"

%{
/* this macro implements parsing from a string, character by character */
#undef YY_INPUT
#define YY_INPUT(buf,result,max_size) \
{\
buf[0]= EXPRparseString[EXPRparseOffset++];\
result= (buf[0]== '\0') ? YY_NULL : 1;\
}

#include <math.h>
#include <iostream>

#if defined (_WIN32) || defined (_WIN64)
#define _USE_MATH_DEFINES
#include <math.h>
#endif

/* the string to be parsed */
char *EXPRparseString;
/* the current position inside that string */
int EXPRparseOffset;

#if defined (_WIN32) || defined (_WIN64)
#include "ExpressionParseTree.hh"
#else
#include "../ExpressionParseTree.hh"
#endif
using namespace EXPR;
#include "ExpressionParser.h"
using namespace std;
#if defined (_WIN32) || defined (_WIN64)
#include "ExpressionHelperFunctions.hh"
#else
#include "../ExpressionHelperFunctions.hh"
#endif

void EXPRerror (const char *s);
%}

%%

(([0-9]+)|([0-9]*\.[0-9]+)|([0-9]+\.[0-9]*))(e\-?[0-9]+)?		{EXPRlval.constant= atof(yytext); return CONSTANT;}
#[0-9]+		{EXPRlval.variable= atoi(yytext+1); return IN_VARIABLE;}
\%[0-9]+	{EXPRlval.variable= atoi(yytext+1); return OUT_VARIABLE;}
[ \t\n]		;
\=\=		{return EQ;}
\!\=		{return NEQ;}
\<\=		{return LEQ;}
\>\=		{return GEQ;}
\!		{return NOT;}
\|\|		{return OR;}
\&\&		{return AND;}
[-+*/(),?:<>]	{return yytext[0];}
(PI)|(pi)|(Pi)	{EXPRlval.constant= M_PI; return CONSTANT;}
[Ee]		{EXPRlval.constant= M_E; return CONSTANT;}
sqrt		{EXPRlval.univariate= &sqrt; return UNIVARIATE;}
exp		{EXPRlval.univariate= &exp; return UNIVARIATE;}
log		{EXPRlval.univariate= &log; return UNIVARIATE;}
sin		{EXPRlval.univariate= &sin; return UNIVARIATE;}
cos		{EXPRlval.univariate= &cos; return UNIVARIATE;}
tan		{EXPRlval.univariate= &tan; return UNIVARIATE;}
asin		{EXPRlval.univariate= &asin; return UNIVARIATE;}
acos		{EXPRlval.univariate= &acos; return UNIVARIATE;}
atan		{EXPRlval.univariate= &atan; return UNIVARIATE;}
atan2		{EXPRlval.bivariate= &atan2; return BIVARIATE;}
floor		{EXPRlval.univariate= &floor; return UNIVARIATE;}
ceil		{EXPRlval.univariate= &ceil; return UNIVARIATE;}
round		{EXPRlval.univariate= &round; return UNIVARIATE;}
frac		{EXPRlval.univariate= &fracOp; return UNIVARIATE;}
abs		{EXPRlval.univariate= &fabs; return UNIVARIATE;}
sgn		{EXPRlval.univariate= &sgnOp; return UNIVARIATE;}
pow		{EXPRlval.bivariate= &pow; return BIVARIATE;}
min		{EXPRlval.bivariate= &minOp; return BIVARIATE;}
max		{EXPRlval.bivariate= &maxOp; return BIVARIATE;}
rand		{return RANDOM;}
.		{ECHO; EXPRerror( "Unexpected character" );}

%%
int EXPRwrap (void) {return 1;}
