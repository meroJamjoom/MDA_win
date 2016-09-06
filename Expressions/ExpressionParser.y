/* ==========================================================================
 * $Id: ExpressionParser.y 243 2008-11-16 01:20:52Z heidrich $
 * a parser for arithmetic expressions
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
%{
#include <iostream>
#include <vector>
#if defined (_WIN32) || defined (_WIN64)
#include "ExpressionParseTree.hh"
#include "ExpressionHelperFunctions.hh"
#else
#include "../ExpressionParseTree.hh"
#include "../ExpressionHelperFunctions.hh"
#endif

using namespace EXPR;
using namespace std;

vector<ExpressionParseTree *> EXPRresult;

int yylex();
void yyerror (const char *s);

%}

%union {
  double constant;
  int variable;
  double (*univariate)(double);
  double (*bivariate)(double,double);
  ExpressionParseTree *node;
}
%token <constant> CONSTANT
%token <variable> IN_VARIABLE
%token <variable> OUT_VARIABLE
%token <univariate> UNIVARIATE
%token <bivariate> BIVARIATE
%token <node> RANDOM
%left ','
%left '?' ':'
%left OR
%left AND
%left '<' '>' GEQ LEQ EQ NEQ
%left '+' '-'
%left '*' '/'
%left NOT
%type <node> sequence expr

%%

sequence: expr		    {$1->optimize(); EXPRresult.push_back($1);}
	| sequence ',' expr {$3->optimize(); EXPRresult.push_back($3);
			     $$= $1;}
	;

expr	: CONSTANT	    {$$= new ExpressionParseTree; 
			     $$->makeConstant( $1 );}
	| RANDOM	    {$$= new ExpressionParseTree;
			     $$->makeRandom();}
	| IN_VARIABLE	    {$$= new ExpressionParseTree;
			     $$->makeInVariable( $1 );}
	| OUT_VARIABLE	    {$$= new ExpressionParseTree;
			     $$->makeOutVariable( $1 );}
        | UNIVARIATE '(' expr ')'
			    {$$= new ExpressionParseTree;
			     $$->makeUnivariate( $1, $3 );}
        | BIVARIATE '(' expr ',' expr ')'
			    {$$= new ExpressionParseTree;
			     $$->makeBivariate( $1, $3, $5 );}
        | '-' expr	    {$$= new ExpressionParseTree;
			     $$->makeUnivariate( &uminusOp, $2 );}
	| NOT expr          {$$= new ExpressionParseTree;
			     $$->makeUnivariate( &notOp, $2 );}
	| expr AND expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &andOp, $1, $3 );}
	| expr OR expr      {$$= new ExpressionParseTree;
			     $$->makeBivariate( &orOp, $1, $3 );}
	| expr EQ expr      {$$= new ExpressionParseTree;
			     $$->makeBivariate( &eqOp, $1, $3 );}
	| expr NEQ expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &neqOp, $1, $3 );}
	| expr LEQ expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &leqOp, $1, $3 );}
	| expr GEQ expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &geqOp, $1, $3 );}
	| expr '<' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &lessOp, $1, $3 );}
	| expr '>' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &greaterOp, $1, $3 );}
	| expr '+' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &plusOp, $1, $3 );}
        | expr '-' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &minusOp, $1, $3 );}
	| expr '*' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &multOp, $1, $3 );}
        | expr '/' expr     {$$= new ExpressionParseTree;
			     $$->makeBivariate( &divOp, $1, $3 );}
	| expr '?' expr ':' expr
			    {$$= new ExpressionParseTree;
			     $$->makeIfThenElse( $1, $3, $5 );}
	| '(' expr ')'      {$$ = $2;}
        ;

%%

void yyerror (const char *s) {cerr << s << endl;}
