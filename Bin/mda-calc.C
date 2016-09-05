// ==========================================================================
// $Id: mda-calc.C 91 2007-05-24 07:16:47Z heidrich $
// a commandline calculator
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================


#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;


int
main( int argc, char *argv[] )
{
  int i, j;
  
  CommandlineParser parser;
  
  // parse options
  int index= 1;
  if( !parser.parse( index, argc, argv ) || index!= argc-1 )
  {
    parser.usage( argv[0], "[options] <expression sequence>" );
    exit( 1 );
  }
  
  // get pixel value
  EXPR::ExpressionSequence expr= EXPR::parse( argv[argc-1] );
  if( expr.size()== 0 )
  {
    parser.usage( argv[0], "[options] <expression sequence>" );
    exit( 1 );
  }
  
  // make sure all register indices are valid
  for( i= 0 ; i< expr.size() ; i++ )
  {
    if( expr[i]->getMaxVariable( true )>= 0 )
    {
      cerr << "No input variables allowed in mda-calc\n";
      exit( 1 );
    }
    if( (j= expr[i]->getMaxVariable( false ))>= i )
    {
      cerr << "Using register " << j << " in expression " << i << endl
	   << "  (only registers 0.." << i-1 << " are valid at this point)\n";
      exit( 1 );
    }
  }
  
  double *regs= new double[expr.size()];
  
  // compute the value of all expessions
  for( i= 0 ; i< expr.size() ; i++ )
    regs[i]= expr[i]->eval( NULL, regs );
  
  // output results
  cout.precision( 14 );
  cout << regs[0];
  for( i= 1 ; i< expr.size() ; i++ )
    cout << ',' << regs[i];
  cout << endl;
  
  return 0;
}
