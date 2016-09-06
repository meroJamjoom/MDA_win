// ==========================================================================
// $Id: mda-poisson.C 638 2011-08-29 19:34:24Z heidrich $
// apply poisson solver to an MDA stream
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

#include "MDA/Array/MDAFileIO.hh"
#include "MDA/LinearAlgebra/PoissonSystem.hh"
#include "MDA/LinearAlgebra/MultigridPreconditioner.hh"
#include "MDA/LinearAlgebra/ConjugateGradient.hh"
#if defined (_WIN32) || defined (_WIN64)
#include "MDA/LinearAlgebra/MultigridPreconditioner.c"
#endif
#include <limits>
#include <list>
#include <algorithm>

using namespace MDA;
using namespace std;

#define MAX_TEXTLEN 1024
char usageText[MAX_TEXTLEN] =
 "<params> [<params> ...] \n"
 "Find the closest scalar field (in the least squares sense) to the provided vector field of gradients. \n"
 "Each parameter set specifies one particular channel to integrate";

#define TYPE double

template<typename T>
inline void setBit(T & flag, unsigned short index, bool value)
{
  flag &= ~(1 << index);
  flag |= value << index;
}

template<typename T>
inline bool getBit(T flag, unsigned index)
{
  return (flag >> index & 1);
}

template<typename T>
inline bool hasChannel(Array<T> & array, unsigned toFind)
{
  if(array[toFind]!=NULL && &(array[toFind][0])!=NULL)
    return true;
  return false;
}


template<typename T>
unsigned addMaskChannel(Array<T> & array, bool valueConstrained)
{
  long i;
  long totalPts;
  unsigned dimension = array.getDimension().vec.size();
  unsigned currDim;
  unsigned boundaryFlag = 0;
  long * nextBoundaryChange = new long[dimension];
  long * offsetTable = new long[dimension];
  unsigned maskCh = array.addChannel();
  offsetTable[0] = 1;
  for(i=1; i<dimension; i++)
    offsetTable[i] = array.getDimension().vec[i-1]*offsetTable[i-1];
  totalPts = offsetTable[dimension-1]*array.getDimension().vec[dimension-1];
  long currOffset;
  int valueConstrainedMask = PoissonSystem<T>::valueMask;
  int gradientConstrainedMask = PoissonSystem<T>::gradientMask;
  int unconstrainedMask = PoissonSystem<T>::unconstrainedMask;
  for(currDim =0; currDim < dimension; currDim++)
  {
    nextBoundaryChange[currDim] = offsetTable[currDim] - 1;
    setBit(boundaryFlag, 2*currDim+1, true);
  }
  for(i=0; i<totalPts; i++)
  {
    for( currDim= 0; currDim < dimension; currDim++)
    {
      currOffset = offsetTable[currDim];
      if( getBit(boundaryFlag, 2*currDim) == true || 
	  getBit(boundaryFlag, 2*currDim+1) == true)
      {
	array[maskCh][0][i] = valueConstrained ? 
	  valueConstrainedMask : gradientConstrainedMask ;
	break;
      }
    }
    array[maskCh][0][0] = valueConstrainedMask;
    if(currDim == dimension) 
      array[maskCh][0][i] = unconstrainedMask;
    for( currDim= 0; currDim < dimension; currDim++)
    {
      if( i==nextBoundaryChange[currDim] )
      {
	currOffset = offsetTable[currDim];
        switch( (boundaryFlag >> 2*currDim) & 3)
	{
	case 0: 
	  setBit(boundaryFlag, 2*currDim, 1);
	  nextBoundaryChange[currDim] += currOffset;
	  break;
	case 1: 
	  setBit(boundaryFlag,2*currDim,0);
	  setBit(boundaryFlag,2*currDim+1,1);
	  nextBoundaryChange[currDim] += currOffset;
	  break;
	case 2: 
	  setBit(boundaryFlag,2*currDim+1,0);
	  nextBoundaryChange[currDim] += (array.getDimension().vec[currDim]-2)*currOffset;
	  break;
	}
      }
    }
    
  }

  delete [] nextBoundaryChange;
  delete [] offsetTable;
  return maskCh;
}

int 
main(int argc, char* argv[])
{ 
  unsigned long i;
  int numChannels, numOldChannels;
  CommandlineParser parser; 

  ChannelList channels;
  ChannelListOption channelOpt(channels,
			       "\tGradient channels (number must match array dimension. Necessarily only if gradient constraints used or generated potential uses gradient values)");
  parser.registerOption( &channelOpt );

  DataType outType= UndefinedType;
  TypeOption typeOption( outType );
  parser.registerOption( &typeOption );

  unsigned constraintCh;
  ScalarOption<unsigned> constraintChOpt( constraintCh, "\tConstraint channel (values at value constrained pixels need only be set)" , "--constraint" , "-c", numeric_limits<unsigned>::min(), numeric_limits<unsigned>::max());
  
  parser.registerOption( &constraintChOpt );

  unsigned targetCh;
  ScalarOption<unsigned> targetChOpt( targetCh, "\tTarget channel (default: create a new channel)","--target",NULL, numeric_limits<unsigned>::min(), numeric_limits<unsigned>::max());  

  parser.registerOption( &targetChOpt );

  unsigned sourceCh;
  ScalarOption<unsigned> sourceChOpt( sourceCh, "\tSource channel (rhs in Poisson Equation)","--source","-s",numeric_limits<unsigned>::min(), numeric_limits<unsigned>::max());
  parser.registerOption( &sourceChOpt );

  double constSource;
  DoubleOption constSourceOpt( constSource, "\tConstant source value","--const-source","-cs",numeric_limits<double>::min(),numeric_limits<double>::max());
  parser.registerOption( &constSourceOpt );

  bool computeDiv;
  BoolOption computeDivOpt( computeDiv, "\tUse the divergence of the vector field as the potential function","--use-div","-div","","");
  parser.registerOption( &computeDivOpt );

  bool multigrid=false;
  BoolOption multigridOpt( multigrid, 
			   "\tUse multigrid preconditioner",
			   "--multigrid","-mg","","");
  parser.registerOption( &multigridOpt );

  unsigned maxNumIter=-1;
  ScalarOption<unsigned> maxNumIterOpt( maxNumIter, "\tMaximum number of CG iterations" , "--numiter" , "-n", numeric_limits<unsigned>::min(), numeric_limits<unsigned>::max());
  
  parser.registerOption( &maxNumIterOpt );


  int gamma = 2; 
  IntOption gammaOpt( gamma, "\tGamma option - this variable controls how many iteration of CG are performed before downsampling the system to a lower resolution for solving","--gamma","-g",numeric_limits<int>::min(),numeric_limits<int>::max());
  parser.registerOption( &gammaOpt );

  double thresholdError = -1;
  DoubleOption thresholdErrorOpt( thresholdError, "\tAmount of error we are willing to tolerate before terminating conjugate gradient solving","--thresh","-t",numeric_limits<double>::min(),numeric_limits<double>::max());
  parser.registerOption( &thresholdErrorOpt );

  
  unsigned maskCh; 
  ScalarOption<unsigned> maskChOpt( maskCh, "\tMask channel. Pixels constrained according to the following criteria):" 
				    "\n\t\t 0: value constrained"
				    "\n\t\t > 0: unconstrained"
				    "\n\t\t < 0: gradient constrained ",
				    "--mask" , NULL, numeric_limits<unsigned>::min(),    numeric_limits<unsigned>::max());

  parser.registerOption( &maskChOpt );

  bool valConstrained;
  BoolOption valConstrainedOpt( valConstrained, "\tUse value constraints on boundaries",  "--value-constrained","-vc","-gradient-constrained","-gc");
  parser.registerOption( &valConstrainedOpt );

  int index = 1;
  Array<TYPE> array;
  if( !array.read() ) { 
    errorCond(1,": Cannot read input MDA stream!\n\n");
    exit( 1 );
  }

  unsigned long totalPts = 1;
  for(i=0; i<array.getDimension().vec.size(); i++)
    totalPts *= array.getDimension().vec[i];
  Vector rhs(totalPts);
  Array<double> solutionArray(array.getDimension());
  solutionArray.addChannel();
  Vector solution(totalPts,&(solutionArray[0][0][0])); 
  PoissonSystem<TYPE> imageSys;
  ConjugateGradient solver;      
  bool addedMask = false;
   
  constraintCh = numeric_limits<unsigned>::max();
  targetCh = numeric_limits<unsigned>::max();
  sourceCh = numeric_limits<unsigned>::max(); 
  computeDiv = false;
  constSource = numeric_limits<double>::min();
  maskCh = numeric_limits<unsigned>::max();
  
  valConstrained = true;
  multigrid=false;
  maxNumIter=-1;
   
  if(!parser.parse( index, argc, argv )){
    errorCond(index > argc-1,"Index overshoot ");
    parser.usage( argv[0], usageText );
    exit( 1 );
  }

  if(thresholdError == -1){
    int maskChannelSize = 1;
    for(int i = 0; i < array.getDimension().vec.size();++i){
      maskChannelSize *= array.getDimension().vec[i];
    }
    int numUnknowns = 0;
    for(int i = 0; i < maskChannelSize; ++i){
      if((*array[maskCh])[i] == 1){
        numUnknowns += 1;  
      }
    }
    thresholdError = numUnknowns*pow((1.0/256.0),2);
  }

  int numberOfDimensions = array.getDimension().vec.size();
  int numberOfChannels = array.getNumChannels();
  if(computeDiv == false){
    errorCond(numberOfChannels == numberOfDimensions+2,"Poisson solver input has missing channel information - check your input to make sure there are no missing constraint, gradient, or mask channels.");    
  }else{
    errorCond(numberOfChannels != 3,"Poisson solver input has missing channel information - check your input to make sure there is no missing constraint, divergence potential field, or mask channel.");
  }

  if( !computeDiv && constSource==numeric_limits<double>::min() 
      && !hasChannel(array,sourceCh) )
    computeDiv = true;

  if( !hasChannel(array,constraintCh) ) {
    errorCond(index > argc-1,"Error, constraint channelcould not be read ");
    exit( 1 );
  }
   
  if( !hasChannel(array,maskCh) ) {
    maskCh = addMaskChannel(array, valConstrained ); 
    addedMask = true;
  }

  imageSys.setup( array, rhs, channels, constraintCh, maskCh, 
		  computeDiv, sourceCh, multigrid);
 
  if (multigrid){
    MultigridPreconditioner* mg = new MultigridPreconditioner(imageSys,gamma);
    solver.setPreconditioner(mg);
    solver.setThreshold(thresholdError);
  }
  
  solver.setMaxNumIter(maxNumIter); 
  solver.solve( imageSys, rhs, solution);
  if( addedMask )
    array.deleteChannel(maskCh);
  if( !hasChannel(array,targetCh) )
    targetCh = array.addChannel();
  double error = 0;
  double maxerror = 0;
  double currErr;
  for( i=0 ; i<totalPts ; i++ ) {
    array[targetCh][0][i] = solution[i];
    currErr = (solution[i]-array[constraintCh][0][i]);
    error += currErr*currErr;
    maxerror = max(fabs(currErr),maxerror);
  }

  if(outType==UndefinedType)
    outType = array.getNativeType();

  array.write( cout, outType );
  return 0;

}
