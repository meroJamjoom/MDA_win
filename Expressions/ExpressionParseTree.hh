// ==========================================================================
// $Id: ExpressionParseTree.hh 243 2008-11-16 01:20:52Z heidrich $
// nodes for a parse tree
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 1995-2007, Wolfgang Heidrich, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef EXPRESSIONS_EXPRESSIONPARSETREE_H
#define EXPRESSIONS_EXPRESSIONPARSETREE_H

/*! \file  ExpressionParseTree.hh
    \brief nodes for a parse tree
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <vector>
#include <iostream>

#include "MDA/Base/CommandlineParser.hh"

namespace EXPR {

  /** \class ExpressionParseTree ExpressionParseTree.hh
      nodes for a parse tree */
  
  class ExpressionParseTree {
    
  public:
    /** types of nodes in the parse tree */
    typedef enum
    {
      Undefined,
      Constant,
      Random,
      InVariable,
      OutVariable,
      UnivariateFunc,
      BivariateFunc,
      IfThenElse
    } Type;
    
    /** default constructor */
    ExpressionParseTree()
    {
      type= Undefined;
      children[0]= children[1]= children[2]= NULL;
    }
    
    /** destructor */
    ~ExpressionParseTree();
    
    /** remove the children and set type to unknown */
    void cleanup();

    /** make the node a constant */
    inline void makeConstant( double val )
    {
      cleanup();
      value.constant= val;
      type= Constant;
    }
    
    /** make the node a random number */
    inline void makeRandom()
    {
      cleanup();
      type= Random;
    }
    
    /** make the node an input variable reference */
    inline void makeInVariable( int val )
    {
      cleanup();
      value.variable= val;
      type= InVariable;
    }
    
    /** make the node an output variable reference */
    inline void makeOutVariable( int val )
    {
      cleanup();
      value.variable= val;
      type= OutVariable;
    }
    
    /** make the node a univariate function */
    inline void makeUnivariate( double (*func)(double),
				ExpressionParseTree *arg )
    {
      cleanup();
      value.univariate= func;
      children[0]= arg;
      type= UnivariateFunc;
    }
    
    /** make the node a biivariate function */
    inline void makeBivariate( double (*func)(double,double),
			       ExpressionParseTree *arg1,
			       ExpressionParseTree *arg2 )
    {
      cleanup();
      value.bivariate= func;
      children[0]= arg1;
      children[1]= arg2;
      type= BivariateFunc;
    }
    
    /** make the node a biivariate function */
    inline void makeIfThenElse( ExpressionParseTree *cond,
				ExpressionParseTree *arg1,
				ExpressionParseTree *arg2 )
    {
      cleanup();
      children[0]= cond;
      children[1]= arg1;
      children[2]= arg2;
      type= IfThenElse;
    }

    /** output the tree */
    void printTree( std::ostream &os= std::cout );

    /** optimize the tree by executing operations only involving constants */
    ExpressionParseTree &optimize();
    
    /** evaluate the tree */
    double eval( double *inVariables, double *outVariables );
    
    /** return the largest input variable index used by the parse tree */
    int getMaxVariable( bool in );
    
  protected:
    
    /** the type of this node */
    Type type;
    
    /** the value of this node */
    union
    {
      double	constant;
      unsigned	variable;
      double	(*univariate)(double);
      double	(*bivariate)(double,double);
    } value;
    
    /** the children of this node */
    ExpressionParseTree *children[3];

  };
  
  /** sequence of expressions */
  typedef std::vector<ExpressionParseTree *> ExpressionSequence;
  
  /** parse the string; return a ointer to the parse tree */
  ExpressionSequence &parse( char *string );
}

namespace MDA {

  using namespace EXPR;
  
  /** \class ExpressionOption ExpressionParseTree.hh
      commandline option providing an arithmetic expression */
  
  class ExpressionOption: public CommandlineOption {
    
  public:
    /** constructor from reference to DataType object */
    ExpressionOption( ExpressionSequence &val, const char *msg,
		      const char *longOpt, const char *shortOpt= NULL )
      : CommandlineOption( msg, longOpt, shortOpt ), value( val )
    {}
    
    /** the actual parsing function */
    virtual bool parse( int &index, int argc, char *argv[] );
    
    /** output usage string */
    virtual void usage( ostream &os= cerr );
    
  protected:
    
    /** reference to the DataType object */
    ExpressionSequence &value;
  };
    
  
} /* namespace */


#endif /* EXPRESSIONS_EXPRESSIONPARSETREE_H */

