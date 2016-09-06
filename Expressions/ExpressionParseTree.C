// ==========================================================================
// $Id: ExpressionParseTree.C 999 2014-05-28 15:07:31Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 1995-2008, Wolfgang Heidrich, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef EXPRESSIONS_EXPRESSIONPARSETREE_C
#define EXPRESSIONS_EXPRESSIONPARSETREE_C

#include <string.h>
#include <math.h>
#include <vector>
#include "ExpressionParseTree.hh"
#include "ExpressionHelperFunctions.hh"

#if defined (_WIN32) || defined (_WIN64)
#define drand48() rand()
#endif


// Unfortunately, the following declarations cannot be in the namespace.
// We use the EXPR prefix as a workaround

/** running the parser */
extern int EXPRparse(void);

/** the result of the parsing operation */
extern std::vector<EXPR::ExpressionParseTree *> EXPRresult;

/** the string to be parsed */
extern char *EXPRparseString;

/** current location in parse string */
#if defined (_WIN32) || defined (_WIN64)
extern int EXPRparseOffset;
#else
extern int *EXPRparseOffset;
#endif

namespace EXPR {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

/** destructor */
ExpressionParseTree::~ExpressionParseTree()
{
  cleanup();
}

/** remove the children and set type to unknown */
void
ExpressionParseTree::cleanup()
{
  if( children[0] )
    delete children[0];
  if( children[1] )
    delete children[1];
  if( children[2] )
    delete children[2];
  children[0]= children[1]= children[2]= NULL;
  type= Undefined;
}


/** optimize the tree by executing operations only involving constants */
ExpressionParseTree &
ExpressionParseTree::optimize()
{
  switch( type )
  {
  case UnivariateFunc:
    children[0]->optimize();
    if( children[0]->type == Constant )
    {
      double val= (*value.univariate)( children[0]->value.constant );
      type= Constant;
      value.constant= val;
      delete children[0];
      children[0]= NULL;
    }
    break;
  case BivariateFunc:
    children[0]->optimize();
    children[1]->optimize();
    if( children[0]->type == Constant && children[1]->type == Constant )
    {
      double val= (*value.bivariate)( children[0]->value.constant,
				      children[1]->value.constant );
      type= Constant;
      value.constant= val;
      delete children[0];
      children[0]= NULL;
      delete children[1];
      children[1]= NULL;
    }
    break;
  case IfThenElse:
    children[0]->optimize();
    children[1]->optimize();
    children[2]->optimize();
    if( children[0]->type == Constant )
    {
      ExpressionParseTree *trueBranch=
	children[0]->value.constant> 0.0 ? children[1] : children[2];
      ExpressionParseTree *falseBranch=
	children[0]->value.constant> 0.0 ? children[2] : children[1];
      
      // copy info from the true branch
      type= trueBranch->type;
      children[0]= trueBranch->children[0];
      children[1]= trueBranch->children[1];
      children[2]= trueBranch->children[2];
      memcpy( &value, &(trueBranch->value), sizeof( value ) );
      
      // delete the false branch
      falseBranch->children[0]=		// so we don't recursively
	falseBranch->children[1]=	// delete nodes that are
	falseBranch->children[2]= NULL; // still required
      delete falseBranch; // this is now just the one node by itself
    }
    break;
  default:
    break;
  }
  
  return *this;
}


/** return the largest variable index used by the parse tree */
int
ExpressionParseTree::getMaxVariable( bool in )
{
  int a, b, c;
  
  switch( type )
  {
  case InVariable:
    if( in )
      return value.variable;
    else
      return -1;
  case OutVariable:
    if( in )
      return -1;
    else
      return value.variable;
  case UnivariateFunc:
    return children[0]->getMaxVariable( in );
  case BivariateFunc:
    a= children[0]->getMaxVariable( in );
    b= children[1]->getMaxVariable( in );
    return a> b ? a : b;
  case IfThenElse:
    a= children[0]->getMaxVariable( in );
    b= children[1]->getMaxVariable( in );
    c= children[2]->getMaxVariable( in );
    return a>b ? (a>c ? a:c) : (b>c ? b:c);
  default:
    return -1;
  }
}


/** evaluate the tree */
double
ExpressionParseTree::eval( double *inVariables, double *outVariables )
{
  switch( type )
  {
  case Constant:
    return value.constant;
  case Random:
    return drand48();
  case InVariable:
    return inVariables[value.variable];
  case OutVariable:
    return outVariables[value.variable];
  case UnivariateFunc:
    return (*value.univariate)( children[0]->eval(inVariables, outVariables) );
  case BivariateFunc:
    return (*value.bivariate)( children[0]->eval( inVariables, outVariables ),
			       children[1]->eval( inVariables, outVariables ));
  case IfThenElse:
    if( children[0]->eval( inVariables, outVariables )> 0.0 )
      return children[1]->eval( inVariables, outVariables );
    else
      return children[2]->eval( inVariables, outVariables );
  default:
    cerr << "Unknown node in parse tree -- cannot evaluate\n";
#if defined (_WIN32) || defined (_WIN64)
	return 0.0/1.0;
#else
    return 0.0/0.0;
#endif
  }
}


/** output the tree */
void
ExpressionParseTree::printTree( ostream &os )
{
  switch( type )
  {
  case Constant:
    os << value.constant;
    break;
  case Random:
    os << "rand";
    break;
  case InVariable:
    os << "#" << value.variable;
    break;
  case OutVariable:
    os << "%" << value.variable;
    break;
  case UnivariateFunc:
#if defined (_WIN32) || defined (_WIN64)
	  if (value.univariate == (UNIVARIATE1)&sqrt)
		  os << "sqrt(";
	  else if (value.univariate == (UNIVARIATE1)&exp)
		  os << "exp(";
	  else if (value.univariate == (UNIVARIATE1)&log)
		  os << "log(";
	  else if (value.univariate == (UNIVARIATE1)&sin)
		  os << "sin(";
	  else if (value.univariate == (UNIVARIATE1)&cos)
		  os << "cos(";
	  else if (value.univariate == (UNIVARIATE1)&tan)
		  os << "tan(";
	  else if (value.univariate == (UNIVARIATE1)&asin)
		  cerr << "asin(";
	  else if (value.univariate == (UNIVARIATE1)&acos)
		  os << "acos(";
	  else if (value.univariate == (UNIVARIATE1)&atan)
		  os << "atan(";
	  else if (value.univariate == (UNIVARIATE1)&floor)
		  os << "floor(";
	  else if (value.univariate == (UNIVARIATE1)&ceil)
		  os << "ceil(";
	  else if (value.univariate == (UNIVARIATE1)&round)
		  os << "round(";
	  else if (value.univariate == &fracOp)
		  os << "frac(";
	  else if (value.univariate == (UNIVARIATE1)&fabs)
		  os << "fabs(";
	  else if (value.univariate == &sgnOp)
		  os << "sgn(";
	  else if (value.univariate == &uminusOp)
		  os << "-(";
	  else if (value.univariate == &notOp)
		  os << "!(";
	  else
		  os << "unknown(";
#else

    if( value.univariate== &sqrt )
      os << "sqrt(";
    else if( value.univariate== &exp )
      os << "exp(";
    else if( value.univariate== &log )
      os << "log(";
    else if( value.univariate== &sin )
      os << "sin(";
    else if( value.univariate== &cos )
      os << "cos(";
    else if( value.univariate== &tan )
      os << "tan(";
    else if( value.univariate== &asin )
      cerr << "asin(";
    else if( value.univariate== &acos )
      os << "acos(";
    else if( value.univariate== &atan )
      os << "atan(";
    else if( value.univariate== &floor )
      os << "floor(";
    else if( value.univariate== &ceil )
      os << "ceil(";
    else if( value.univariate== &round )
      os << "round(";
    else if( value.univariate== &fracOp )
      os << "frac(";
    else if( value.univariate== &fabs )
      os << "fabs(";
    else if( value.univariate== &sgnOp )
      os << "sgn(";
    else if( value.univariate== &uminusOp )
      os << "-(";
    else if( value.univariate== &notOp )
      os << "!(";
    else
      os << "unknown(";
#endif
    children[0]->printTree( os );
    os << ')';
    break;
  case BivariateFunc:
#if defined (_WIN32) || defined (_WIN64)
	  // two types of bivariates: methods and infix operators
	  if (value.bivariate == (BIVARIATE1)&atan2 || value.bivariate == (BIVARIATE1)&pow ||
		  value.bivariate== &minOp || value.bivariate== &maxOp )
	  {
		  // methods
		  if (value.bivariate == (BIVARIATE1)&atan2)
			  os << "atan2(";
		  else if (value.bivariate == (BIVARIATE1)&pow)
			  os << "pow(";
		  else if( value.bivariate== &minOp )
			  os << "min(";
		  else
			  os << "max(";

		  children[0]->printTree( os );
		  os << ',';
		  children[1]->printTree( os );
		  os << ')';
	  }
#else
    // two types of bivariates: methods and infix operators
    if( value.bivariate== &atan2 || value.bivariate== &pow ||
	value.bivariate== &minOp || value.bivariate== &maxOp )
    {
      // methods
      if( value.bivariate== &atan2 )
	os << "atan2(";
      else if( value.bivariate== &pow )
	os << "pow(";
      else if( value.bivariate== &minOp )
	os << "min(";
      else
	os << "max(";
      
      children[0]->printTree( os );
      os << ',';
      children[1]->printTree( os );
      os << ')';
    }
#endif
    else
    {
      // infix operators
      if( children[0]->type== BivariateFunc )
      {
	// use brackets for left child if it is a bivariate function
	// (it could be an infix operator...)
	os << '(';
	children[0]->printTree( os );
	os << ')';
      }
      else
	children[0]->printTree( os );
	
      if( value.bivariate== &plusOp )
	os << '+';
      else if( value.bivariate== &minusOp )
	os << '-';
      else if( value.bivariate== &multOp )
	os << '*';
      else if( value.bivariate== &divOp )
	os << '/';
      else if( value.bivariate== &lessOp )
	os << '<';
      else if( value.bivariate== &greaterOp )
	os << '>';
      else if( value.bivariate== &eqOp )
	os << "==";
      else if( value.bivariate== &neqOp )
	os << "!=";
      else if( value.bivariate== &leqOp )
	os << "<=";
      else if( value.bivariate== &geqOp )
	os << ">=";
      else if( value.bivariate== &andOp )
	os << "&&";
      else if( value.bivariate== &orOp )
	os << "||";
      
      if( children[1]->type== BivariateFunc )
      {
	// use brackets for right child if it is a bivariate function
	// (it could be an infix operator...)
	os << '(';
	children[1]->printTree( os );
	os << ')';
      }
      else
	children[1]->printTree( os );
    }
    break;
  case IfThenElse:
    os << "(";
    children[0]->printTree( os );
    os << ") ? (";
    children[1]->printTree( os );
    os << ") : (";
    children[2]->printTree( os );
    os << ")";
    break;
  default:
    os << "Unknown node type\n";
  }
}

  
ExpressionSequence &
parse( char *string )
{
  EXPRparseString= string;
  EXPRparseOffset= 0;
  EXPRresult.clear();
  EXPRparse();
  return EXPRresult;
}

} /* namespace */


namespace MDA {

/** the actual parsing function */
bool
ExpressionOption::parse( int &index, int argc, char *argv[] )
{
  if( index>  argc-1 )
    return false;
  value= EXPR::parse( argv[index++] );
  return value.size()> 0;
}


/** output usage string */
void
ExpressionOption::usage( ostream &os )
{
  os << "  ";
  if( shortTxt!= NULL )
    os << shortTxt << " | ";
  os << longTxt << " <expr>,<expr>,...,<expr>\n"
     << helpTxt
     << "\tCurrently: ";
  for( unsigned i= 0 ; i< value.size() ; i++ )
  {
    if( i> 0 )
      os << ", ";
    value[i]->printTree( os );
  }
  os << "\n\n";
}

} /* namespace */

#endif /* EXPRESSIONS_EXPRESSIONPARSETREE_C */

