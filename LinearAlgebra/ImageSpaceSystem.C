// ==========================================================================
// $Id: ImageSpaceSystem.C 999 2014-05-28 15:07:31Z heidrich $
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

#ifndef IMAGESPACESYSTEMS_IMAGESPACESYSTEM_C
#define IMAGESPACESYSTEMS_IMAGESPACESYSTEM_C

#include "ImageSpaceSystem.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;
 

  /** destructor */
  ImageSpaceSystem::~ImageSpaceSystem()
  {
    if( nonZeroes!= NULL )
      delete [] nonZeroes;
    if( supportSize!= NULL )
      delete [] supportSize;
    if (offsetTable!=NULL) 
      delete [] offsetTable;
    if (neighbourConfig!=NULL)
      delete [] neighbourConfig;
    if( elems!= NULL )
      {
	for( unsigned i= 0 ; i< numCases ; i++ )
	  delete [] elems[i];
	delete [] elems;
      }
  }
    

  /** multiply a vector with the image-space system 
   * this method deals with interior pixels and constant pixels
   * (diagonal entry only). Boundary cases need to be handled by
   * subclasses.
   */
  Vector &
  ImageSpaceSystem::rightMultiply( const Vector &v, Vector &res ) const
  {
    unsigned dim= v.getSize();
    unsigned config;
    unsigned num;
    double h;
  
    for( long i= 0 ; i< dim ; i++ )
      if( nonZeroes[i]== 0 )
	// constant pixel (diagonal element =1, all others =0 )
	res[i]= v[i];
      else
	{
	  h= 0.0;
	  config= nonZeroes[i];	// bit vector of pixel configuration
	  num= supportSize[config];	// number of pixels in support
	  for( unsigned j= 0 ; j <= num ; j++ )
	    h+= elems[config][j].second * v[i+elems[config][j].first];
	  res[i]= h;
	}
  
    return res;
  }

  void ImageSpaceSystem::setupCaseTable()
  {
    unsigned i, j;
  
    // compute number of different configurations of local
    // neighborhoods
    unsigned numNeighbors = 2*dimension; 
    //unsigned numNeighbors = pow(dimension,2); 
    ImageSpaceSystem::numCases = 1 << numNeighbors;
  
    unsigned nonZero, currDim;
  
    // initialize the number of non-zero elements for each case
    if( ImageSpaceSystem::supportSize!= NULL )
      delete [] ImageSpaceSystem::supportSize;
    ImageSpaceSystem::supportSize= 
      new unsigned char[ImageSpaceSystem::numCases];
    // count # one-bits for each case
    for(i = 0 ; i< ImageSpaceSystem::numCases ; i++ ){
      j = i;
      ImageSpaceSystem::supportSize[j]= 0;
      for(; i!=0; i>>=1) {
	if(i & 1)
	  ImageSpaceSystem::supportSize[j]++;
      }	
      i = j;
    }
  
    // initialize elems (depends on array dimensions)
    if( ImageSpaceSystem::elems!= NULL ){
      for( i= 0 ; i< ImageSpaceSystem::numCases ; i++ )
	delete [] ImageSpaceSystem::elems[i];
      delete [] ImageSpaceSystem::elems;
    }
  
    ImageSpaceSystem::elems= 
      new pair<long,double>*[ImageSpaceSystem::numCases];
    for( i= 0; i < ImageSpaceSystem::numCases; i++)
      ImageSpaceSystem::elems[i]= 
	new pair<long,double>[ImageSpaceSystem::supportSize[i] + 1];
  
    for( i= 0; i < ImageSpaceSystem::numCases; i++){
      //diagonal entry
      ImageSpaceSystem::elems[i][0].first= 0;
      ImageSpaceSystem::elems[i][0].second= ImageSpaceSystem::supportSize[i];
      
      //go through each dimension, determining the adjacent pixels 
      for(nonZero= 1, currDim = 0; currDim < dimension; currDim++){
	if(getBit(i,2*currDim)==1){
	  ImageSpaceSystem::elems[i][nonZero].first = offsetTable[currDim];
	  ImageSpaceSystem::elems[i][nonZero].second = -1;
	  nonZero++;
	}
	if(getBit(i,2*currDim+1)==1){
	  ImageSpaceSystem::elems[i][nonZero].first = -offsetTable[currDim];
	  ImageSpaceSystem::elems[i][nonZero].second = -1;
	  nonZero++;
	}
      }
    } 
  }



  /** determine the neighbourhood configuration for each pixel */
  void ImageSpaceSystem::computeNeighborCounts ( )
  {
    // initialize nonZeroes (depends on specific array)
    if(ImageSpaceSystem::nonZeroes != NULL)
      delete [] ImageSpaceSystem::nonZeroes;
    ImageSpaceSystem::nonZeroes= new unsigned [totalPts];
    if(neighbourConfig != NULL)
      delete [] neighbourConfig;
    neighbourConfig = new unsigned[totalPts];

    unsigned currDim, neighbourFlag;
    unsigned boundaryFlag = 0;
    long * nextBoundaryChange = new long[dimension];
    long currOffset;
    //iteration starts at the negative boundary of all dimensions
    //(i.e. 0) so initialize helper variables to reflect this
    for( currDim= 0; currDim < dimension; currDim++){
      nextBoundaryChange[currDim]= offsetTable[currDim] - 1;
      //initial point is on the negative boundary of each dimension
      setBit(boundaryFlag, 2*currDim+1, 1);
    }
    for( long i= 0; i < totalPts; i++ ){
      neighbourFlag= 0;
      for( currDim= 0; currDim < dimension; currDim++ ){
	currOffset = offsetTable[currDim];
	//Check for a neighbour in the positive direction, ensuring
	//point is away away from the positive boundary in that
	//dimension
	if( getBit(boundaryFlag, 2*currDim) == 0 )
	  setBit(neighbourFlag, 2*currDim, true);
	//Check for a neighbour in the negative direction, in a
	//similar fashion
	if( getBit(boundaryFlag, 2*currDim+1) == 0 )
	  setBit(neighbourFlag, 2*currDim+1,true);
      }
    
      ImageSpaceSystem::nonZeroes[i] = neighbourFlag;
      neighbourConfig[i] = neighbourFlag;
      //Maintain the array which determines the index at which the
      //boundary will change in a particular dimension
      for( currDim= 0; currDim < dimension; currDim++)
	if(i==nextBoundaryChange[currDim]){
	  currOffset = offsetTable[currDim];
	  switch( (boundaryFlag >> 2*currDim ) & 3){
	  case 0: //not on any boundary currently, moving to the positive one
	    setBit(boundaryFlag ,2*currDim,1);
	    nextBoundaryChange[currDim] += currOffset;
	    break;
	  case 1: //on the positive boundary, moving to the negative one
	    setBit(boundaryFlag,2*currDim,0);
	    setBit(boundaryFlag,2*currDim+1,1);	  
	    nextBoundaryChange[currDim] += currOffset;
	    break;
	  case 2: //on the negative boundary, moving off of it	  
	    setBit(boundaryFlag,2*currDim+1,0);
	    nextBoundaryChange[currDim] += (dim.vec[currDim]-2)*currOffset;
	    break;
	  default:
	    break;
	  }     
	}
    }
    delete [] nextBoundaryChange;
    nextBoundaryChange = NULL;
  }

  /** setup the offset table */
  void ImageSpaceSystem::setupOffsetTable(){
    int i;
    // initialize the offset table (for moving one index in a 
    // particular dimension)  
    if(offsetTable!= NULL) delete [] offsetTable;
    offsetTable= new long[dimension];
    offsetTable[0]= 1;
    for( i= 1; i < dimension; i++)
      offsetTable[i] = offsetTable[i-1]*dim.vec[i-1];
  }
} /* namespace */

#endif /* IMAGESPACESYSTEMS_IMAGESPACESYSTEM_C */

