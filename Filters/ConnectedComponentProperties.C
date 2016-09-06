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

#ifndef CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_C
#define CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_C

#include "ConnectedComponentProperties.hh"
#include "minmax.h"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


  /** default constructor. Takes the connected components image where the connected component shall be in channel 0 */
  template<class T>
  ConnectedComponentProperties<T>::ConnectedComponentProperties(Array<T> *connComp)
  {
      this->connComp = connComp;

      // Get size and dimensions of Array
      srcDim = connComp->getDimension();

      numEntries= 1;
      for( int l= srcDim.vec.size()-1 ; l>= 0 ; l-- )
	numEntries*= srcDim.vec[l];

      // Initialize property variables
      computedProps = NoProp;

      // Set all members to initial values
      numLabels = -1;
      bboxes = new std::vector<LongRangeList*>();
      area = new std::vector<unsigned long>();
      padding = new std::vector<unsigned int>();
      segments = new std::vector< Array<T>* >();
      eulerNum = new std::vector<double>();
      regionEntries = new std::vector< unsigned long >();
        
      return;
  }

  /** default destructor */
  template<class T>
  ConnectedComponentProperties<T>::~ConnectedComponentProperties()
  {
      delete bboxes;
      delete area;
      delete segments;
      delete eulerNum;
      delete padding;
      delete regionEntries;
  }

  /** returns bool if the specified property is computed. */
  template<class T>
  bool ConnectedComponentProperties<T>::getPropertyComputed( unsigned int property )
  {
      if( property == (unsigned)NoProp && computedProps == (unsigned)NoProp )
	  return true;

      return (property & computedProps);
  }

  /** computes the number of labels  */
  template<class T>
  void ConnectedComponentProperties<T>::computeNumLabels()
  {
      T result = 0.0;
      
      // read data scanline, and apply reduction ops to all pixels
      unsigned long scanLineSize = srcDim.vec[0];
      unsigned long numScanlines = numEntries / scanLineSize;

      T* dataScanline= new T[scanLineSize ];
      double* inScanline= new double[scanLineSize ];
      for( unsigned long i= 0 ; i< numScanlines ; i++ )
      {
	// Read next scanline and convert
	connComp->getScanline(i, dataScanline);

	// apply reduction operator to all pixels in scanline
	for( unsigned long j= 0 ; j< scanLineSize ; j++ )
	{
	    result = max( result, dataScanline[j] );
	}
      }

      numLabels = (long)(ceil(result)); 

      //Set bitvector
      computedProps |= (unsigned)NumLabels;

      return;
  }

  /** returns the number of labels if computed, else return -1 */
  template<class T>
  int ConnectedComponentProperties<T>::getNumLabels()
  {
      if( ! getPropertyComputed( NumLabels) )
      {
	  std::cerr << "Error in ConnectedComponentProperties: numLabels is queried"
		    << "but has not been computed." << std::endl;
      }
      return numLabels;
  }

  /** computes the boundingboxes and / or area */
  template<class T>
  void ConnectedComponentProperties<T>::computeIndependentProps( bool computeBBox, bool computeArea)
  {
      //Exit if number of labels not computed
      if( ! getPropertyComputed( NumLabels) )
      {
	  std::cerr << "Error in ConnectedComponentProperties: computeIndependentProps() "
		    << "needs to have numLabels computed in advance." << std::endl;
	  return;
      }

      //####Setup for Boundingboxes
      CoordinateVector srcPos;
      if( computeBBox )
      {
	for(  unsigned int i = 0 ; i < srcDim.vec.size() ; i++ )
	  srcPos.vec.push_back( 0 );

	//Setup boundingboxes
	bboxes->clear();
	for( unsigned long i = 0; i <= numLabels; i++ )
	{
	  bboxes->push_back( new LongRangeList( srcDim.vec.size() ));
	}
      }

      //####Setup for Area
      if( computeArea )
      {
	//Set all areas to zero
	area->clear();
	for( unsigned long i = 0; i <= numLabels; i++ )
	{
	  area->push_back( 0 );
	}
      }

      T* data = &((*((*connComp)[0]))[0]);
      long currRegion;
      long mini, maxi;
      for(  unsigned long i = 0 ; i < numEntries ; i++ )
      {
	//Typeconvert
	currRegion = (long)(data[i]);

	//####Process all regions for BBox computation
	if( computeBBox )
	{
	  if( currRegion > 0 && currRegion <= numLabels )
	  {
	      //Update bbox for current region
	      for( unsigned int j= 0 ; j< srcDim.vec.size() ; j++ )
	      {
		  mini = ((bboxes->at(currRegion))->vec[j]).val.first;
		  maxi = ((bboxes->at(currRegion))->vec[j]).val.second;
		  //Compare side on axis j
		  if( srcPos.vec[j] + 1 > maxi )
		      maxi = srcPos.vec[j] + 1;
		  if( srcPos.vec[j] < mini )
		      mini = srcPos.vec[j];

		  (bboxes->at(currRegion))->vec.at(j) = LongRange( mini, maxi);
	      } 
	  }

	  // update position in source array
	  for( unsigned int j= 0 ; j< srcDim.vec.size() ; j++ )
	  {
	    srcPos.vec[j]++;
	    if( srcPos.vec[j] >= srcDim.vec[j] )
		srcPos.vec[j] = 0;
	    else
	      break;
	  }
	}
	
	//####Process all regions for Area computation
	if( computeArea )
	{
	  if( currRegion > 0 && currRegion <= numLabels )
	    area->at(currRegion)++;
	}
      }

      //Set bitvector
      if( computeBBox )
	  computedProps |= (unsigned)BoundingBox;

      if( computeArea )
	  computedProps |= (unsigned)Area;

      return;
  }

  /** puts every boundingbox region in a segment array with a symetrical padding if in regions vector.
      if regions are not set, every boundingbox region is put into a segment array.*/
  template<class T>
  void ConnectedComponentProperties<T>::computeSegments( int pad, std::vector<unsigned int>* regions )
  {
       //Exit if bboxes not computed
      if( !  getPropertyComputed( BoundingBox)  || ! getPropertyComputed( Area ) )
      { 
	  std::cerr << "Error in ConnectedComponentProperties: computeSegmentsAreaRange() "
		    << "needs to have boundingbox and area computed in advance." << std::endl;
	  return;
      }

      //####Initialize segment arrays
      segments->clear();
      segments->push_back( new Array<T>() );
      padding->clear();
      padding->push_back( 0 );
      CoordinateVector segmentSize;
      Array<T> * currArray;

      regionEntries->clear();
      for( int i = 0; i <= numLabels; i++ )
	  regionEntries->push_back(0);

      bool noSegment;
      bool found;
      for( int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  //Retrieve segment size
	  segmentSize.vec.clear();
	  noSegment = false;
	  regionEntries->at(currRegion) = 1;

	  // Try to find currRegion in regions vector
	  found = false;
	  for( unsigned int j = 0; j < regions->size(); j++)
	  {
	      if( regions->at(j) == currRegion )
	      {
		  found = true;
		  break;
	      }
	  }
	  if( !found )
	      noSegment = true;

	  long bbox_side;
	  for(  unsigned int j = 0 ; j < srcDim.vec.size() && !noSegment ; j++ )
	  {
	      if( ((bboxes->at(currRegion))->vec[j]).val.first >= ((bboxes->at(currRegion))->vec[j]).val.second )
	      {
		  noSegment = true;
		  break;
	      }
	      bbox_side = ((bboxes->at(currRegion))->vec[j]).val.second - ((bboxes->at(currRegion))->vec[j]).val.first;
	      segmentSize.vec.push_back(  bbox_side + pad*2);
	      regionEntries->at(currRegion) *= segmentSize.vec.back();
	  }
	  
	  //Add segment if area not okay else add NULL
	  if( !noSegment)
	  {
	      currArray = new Array<T>( segmentSize );
	      currArray->addChannel();
	      padding->push_back( abs(pad) );
	  }else
	  {
	      currArray = NULL;
	      padding->push_back( 0 );
	  }

	  segments->push_back( currArray );
      }
      
      //####Copy data to segmented arrays
      unsigned long srcArrayPos;
      T* data = &((*((*connComp)[0]))[0]);

      CoordinateVector srcDimProd;
      srcDimProd.vec.push_back(1);
      for( unsigned int i = 1; i < srcDim.vec.size(); i++)
       	srcDimProd.vec.push_back(srcDim.vec[i-1] * srcDimProd.vec[i-1]);

      //Position in the region
      CoordinateVector regionPos;
      for(  unsigned int i = 0 ; i < srcDim.vec.size() ; i++ )
	  regionPos.vec.push_back( 0 );

      //Iterate over all labels
      T* regionData;
      double pixel;
      CoordinateVector regDimProd;
      for( int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  if( segments->at(currRegion) == NULL )
	      continue;

	  //Fetch regiondata
	  regionData = &((*((*(segments->at(currRegion)))[0]))[0]);
	  
	  //Compute padding and indices
	  for(  unsigned int j = 0 ; j < srcDim.vec.size() ; j++ )
	      regionPos.vec[j] = 0;

	  //Iterate over the region
	  long bbox_side;
	  for( unsigned long j = 0; j < regionEntries->at(currRegion); j++ )
	  {
	      // calculate position in source array
	      srcArrayPos = 0;
	      for( int k= srcDim.vec.size() - 1 ; k >=0 ; k-- )
		 srcArrayPos += (((bboxes->at(currRegion))->vec[k]).val.first - pad + regionPos.vec[k])*srcDimProd.vec[k] ;

	      // test if out of range due to padding
	      if( srcArrayPos < 0 || srcArrayPos >= numEntries )
	      {
		  regionData[j] = 0.0;
	      }else
	      {
		  // copy pixel value
		  pixel = (double)(data[srcArrayPos]);
		  regionData[j] = ((pixel == (double)currRegion) ? 1.0 : 0.0);
	      }	      

	      // update position in region array
	      for( unsigned int k= 0 ; k< srcDim.vec.size() ; k++ )
	      {
		regionPos.vec[k]++;
		bbox_side = ((bboxes->at(currRegion))->vec[k]).val.second - ((bboxes->at(currRegion))->vec[k]).val.first;
		if( regionPos.vec[k] >= bbox_side + 2*pad)
		    regionPos.vec[k] = 0;
		else
		    break;
	      }
	  }
      }

      //Set bitvector
      computedProps |= (unsigned) Segments;

      return;
  }
   
  /** reassembles the image from a the segments with indices in a vector*/
  template<class T>
  Array<T>* ConnectedComponentProperties<T>::reassemble( std::vector<unsigned int>* regions )
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Segments shall be reassembled "
		    << "but segments have not been computed in advance." << std::endl;
	  return NULL;
      }

      //Create array with same dimensions as connected components image
      Array<T>* result = new Array<T>( srcDim );
      result->addChannel();
      long srcArrayPos;
      T* data = &((*((*result)[0]))[0]);

      CoordinateVector srcDimProd;
      srcDimProd.vec.push_back(1);
      for( unsigned int i = 1; i < srcDim.vec.size(); i++)
	srcDimProd.vec[i]=srcDim.vec[i-1] * srcDimProd.vec[i-1];

      //###Iterate over all regions
      //Position in the region
      CoordinateVector regionPos;
      for(  unsigned int i = 0 ; i < srcDim.vec.size() ; i++ )
	  regionPos.vec.push_back( 0 );

      //Iterate over all labels
      T* regionData;
      double oldValue;
      CoordinateVector regDimProd;
      bool found;
      for( unsigned int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  // Try to find currRegion in regions vector
	  found = false;
	  
	  if( regions != NULL )
	  {
	    for( unsigned int j = 0; j < regions->size(); j++)
	    {
		if( regions->at(j) == currRegion )
		{
		    found = true;
		    break;
		}
	    }
	  }else
	    found = true;

	  // Next Region if not found or not valid
	  if( !found || (segments->at(currRegion) == NULL) )
	      continue;

	  //Fetch regiondata
	  regionData = &((*((*(segments->at(currRegion)))[0]))[0]);
	  
	  //Compute indices
	  for(  unsigned int j = 0 ; j < srcDim.vec.size() ; j++ )
	    regionPos.vec[j] = 0;

	  //Iterate over the region
	  long bbox_side;
	  for( unsigned long j = 0; j < regionEntries->at(currRegion); j++ )
	  {
	      // calculate position in source array
	      srcArrayPos = 0;
	      for( int k= srcDim.vec.size() - 1 ; k >=0 ; k-- )
		 srcArrayPos += (((bboxes->at(currRegion))->vec[k]).val.first - padding->at(currRegion) + regionPos.vec[k])*srcDimProd.vec[k] ;

	      // test if out of range due to padding
	      if( srcArrayPos < 0 || srcArrayPos >= numEntries )
		continue;

	      // copy pixel value
	      oldValue = (double)(data[srcArrayPos]);
	      data[srcArrayPos] = oldValue + (double)(regionData[j]);

	      // update position in region array
	      for( unsigned int k= 0 ; k< srcDim.vec.size() ; k++ )
	      {
		regionPos.vec[k]++;
		bbox_side = ((bboxes->at(currRegion))->vec[k]).val.second - ((bboxes->at(currRegion))->vec[k]).val.first;
		if( regionPos.vec[k] >= bbox_side + 2*padding->at(currRegion))
		    regionPos.vec[k] = 0;
		else
		    break;
	      }
	  }

      }
 
      return result;
  }

  /** compute eulernumbers of all segments  */
  template<class T>
  void ConnectedComponentProperties<T>::computeEulernumbers()
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Eulernumber shall be computed "
		    << "but need segments to be computed in advance." << std::endl;
	  return;
      }

      // Setup eulernumber container
      eulerNum->clear();
      for( int i = 0; i <= numLabels; i++ )
	  eulerNum->push_back(0.0);

      // Create the eulernumberfilter
      Filter<T> *filter;
      FilterFactory<T> factory;
      FilterType filterType= EulerNumberFiltering; 
      filter= factory.create( filterType );
      
      // Setup filter behaviour for 2D
      BoundaryMethod boundary= Clamp;

      ChannelList channels;
      channels.vec.push_back(0);

      AxisList axes;
      axes.vec.push_back(0);
      axes.vec.push_back(1);
      
      // Temporal Copy of a segment
      Array<T> * currSegment;
      //Loop variables
      double euler;
      CoordinateVector segmentDim;
      unsigned int segmentEntries;
      unsigned int segmentScanlines;
      T* segmentDataScanline;
      unsigned int scanlineLength;

      for( int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  if( segments->at(currRegion) == NULL )
	      continue;

	  // Copy to temporal array
	  currSegment = new Array<T>( (segments->at(currRegion))->getDimension() );
	  currSegment->addChannel();
	  
	  segmentDim = (segments->at(currRegion))->getDimension();
	  segmentEntries= 1;
	  for( int l= segmentDim.vec.size()-1 ; l>= 0 ; l-- )
	    segmentEntries*= segmentDim.vec[l];
	  
	  segmentScanlines = segmentEntries / segmentDim.vec[0];
	  scanlineLength = segmentDim.vec[0];
	  segmentDataScanline= new T[ scanlineLength ];
	  for( unsigned long i= 0 ; i< segmentScanlines ; i++ )
	  {
	    // Read next scanline and copy into new Array
	    (segments->at(currRegion))->getScanline(i, segmentDataScanline);
	    currSegment->updateScanline(i, segmentDataScanline);
	  }
      
	  // apply the filter
	  filter->apply( *currSegment, boundary, channels, axes );

	  //Compute actual eulernumber
	  euler = 0.0;
	  for( unsigned long i= 0 ; i< segmentScanlines ; i++ )
	  {
	    // Read next scanline and copy into new Array
	    currSegment->getScanline(i, segmentDataScanline);
	    // apply reduction operator to all pixels in scanline
	    for( unsigned long j= 0 ; j< scanlineLength ; j++ )
		euler += (double)(segmentDataScanline[j]);
	  }
	  delete[] segmentDataScanline;
	  
	  //Clean up temporal array
	  delete currSegment;

	  //Insert computed Eulernumber
	  eulerNum->at(currRegion) = euler;
      }

      // clean up filter
      delete filter;

      //Set bitvector
      computedProps |= EulerNumber;
      return;
  }

  /** computes the filled region of all valid segments if padded at least 1 pixel  */
  template<class T>
  void ConnectedComponentProperties<T>::computeFilledRegion()
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Filled region shall be computed "
		    << "but need segments to be computed in advance." << std::endl;
	  return;
      }

      // Create the floodfill filter
      Filter<T> *filter = new FloodFill<T>(1.0, 1.0);
      
      // Setup filter behaviour for 2D
      BoundaryMethod boundary= Clamp;

      ChannelList channels;

      AxisList axes;
      axes.vec.push_back(0);
      axes.vec.push_back(1);
      
      //Iterate over all regions
      Array<T>* currSegment;
      T* regionData;
      CoordinateVector segmentDim;
      unsigned int segmentEntries;
      for( int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  if( segments->at(currRegion) == NULL )
	      continue;

	  // No convex hull if padding not enabled
	  if( padding->at(currRegion) < 1 )
	  {	
	      std::cerr << "Could not compute convex hull of region " << currRegion
			<< ". No padding was done." << std::endl;
	      continue;
	  }

	  // Add channel for seed
	  currSegment = (segments->at(currRegion));
	  currSegment->addChannel();

	  segmentDim = currSegment->getDimension();
	  segmentEntries= 1;
	  for( int l= segmentDim.vec.size()-1 ; l>= 0 ; l-- )
	    segmentEntries*= segmentDim.vec[l];

	  regionData = &((*((*currSegment)[1]))[0]);
	  regionData[0] = 1.0;
	  for( unsigned long i= 1 ; i< segmentEntries ; i++ )
	    regionData[i] = 0.0;
	  currSegment->swapChannels(0,1);
      
	  // apply the filter
	  channels.vec.clear();
	  channels.vec.push_back(0);
	  channels.vec.push_back(1);
	  filter->apply( *currSegment, boundary, channels, axes );

	  //invert image
	  regionData = &((*((*currSegment)[0]))[0]);

	  for( unsigned long i= 0 ; i< segmentEntries ; i++ )
	    regionData[i] = (double)(regionData[i]) == 1.0? 0.0 : 1.0;
      }

      // Delete all added channels
      for( int currRegion = 1; currRegion <= numLabels; currRegion++ )
      {
	  if( segments->at(currRegion) == NULL || padding->at(currRegion) < 1 )
	      continue;
	  else
	      (segments->at(currRegion))->deleteChannel(1);
      }

      // clean up filter
      delete filter;

      //Set bitvector
      computedProps |= FilledRegion;
      return;
  }

  /** computes the outline (silhouette edge) of all valid segments in 2D if padded at least 1 pixel */
  template<class T>
  void ConnectedComponentProperties<T>::computeOutline()
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Outline shall be computed "
		    << "but need segments to be computed in advance." << std::endl;
	  return;
      }

      if( srcDim.vec.size() != 2 )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Oultline shall be computed "
		    << "but can only be computed on 2D images." << std::endl;
	  return;
      }
    
      //Set up circle offset vector
      IntRangeList circleOffset;
      circleOffset.vec.push_back(IntRange(-1, 0));
      circleOffset.vec.push_back(IntRange(-1, 1));
      circleOffset.vec.push_back(IntRange( 0, 1));
      circleOffset.vec.push_back(IntRange( 1, 1));
      circleOffset.vec.push_back(IntRange( 1, 0));
      circleOffset.vec.push_back(IntRange( 1,-1));
      circleOffset.vec.push_back(IntRange( 0,-1));
      circleOffset.vec.push_back(IntRange(-1,-1));

      CoordinateVector pPos(2);
      CoordinateVector qPos(2);
      CoordinateVector startPos(2);
      int circleStart;

      // Iterate over all regions
      // Temporal Copy of a segment
      Array<T> * resultSegment;
      T* resultData;
      unsigned int srcPos;

      Array<T>* currSegment;
      T* regionData;
      for( long currRegion = 1 ; currRegion <= numLabels; currRegion++ )
      {
	  if( segments->at(currRegion) == NULL )
	      continue;

	  // No outline if padding not enabled
	  if( padding->at(currRegion) < 1 )
	  {	
	      std::cerr << "Could not compute outline of region " << currRegion
			<< ". No padding was done." << std::endl;
	      continue;
	  }

	  //Create resulting segment and clear it
	  resultSegment = new Array<T>( (segments->at(currRegion))->getDimension() );
	  resultSegment->addChannel();
	  resultData = &((*((*resultSegment)[0]))[0]);
	  for( unsigned long j = 0; j < regionEntries->at(currRegion); j++ )
	      resultData[j] = 0.0;

	  //Fetch regiondata
	  currSegment = segments->at(currRegion);
	  regionData = &((*((*currSegment)[0]))[0]);
	  
	  //First, try to find an edge constellation
	  qPos.vec[0] = 0;
	  qPos.vec[1] = 0;
	  for( unsigned long j = 0; j < regionEntries->at(currRegion); j++ )
	  {
	      if( (double)(regionData[j]) > 0.5 )
		  break;

	      // update position in region array
	      qPos.vec[0]++;
	      if( qPos.vec[0] >= ((segments->at(currRegion))->getDimension()).vec[0] )
	      {
		  qPos.vec[1]++;
		  qPos.vec[0] = 0;
	      }
	  }
      
	  //First p is now on the left of qPos as padding is at least one pixel
	  pPos.vec[0] = qPos.vec[0] - 1;
	  pPos.vec[1] = qPos.vec[1];

	  //Save start position
	  startPos.vec[0] = qPos.vec[0];
	  startPos.vec[1] = qPos.vec[1];

	  //Keep track of samples for deadlock detection
	  Vector newVec( 2, 0.0 ), zeroVec( 2, 0.0 );
	  SampleVector outline;
	  //Vector of number of visits of outline points
	  vector<long> visited;

	  //Repeat until outline completed or until upper limit is achieved
	  for( unsigned long lineLength = 0 ; lineLength < ULONG_MAX; lineLength++ )
	  {
      
	      //End if line closed or deadlock
	      if( lineLength >= 1 && startPos.vec[0] == qPos.vec[0] && startPos.vec[1] == qPos.vec[1] )
		  break;

	      //Position in start array
	      srcPos = pPos.vec[1] * ((segments->at(currRegion))->getDimension()).vec[0] + pPos.vec[0];

	      //Create new sample for outline
	      newVec.copy( zeroVec );
	      newVec.set( 0, (double)(pPos.vec[0])); newVec.set( 1, (double)(pPos.vec[1]));

	      //Try to trace back already found outline samples
	      bool found = false;
	      bool deadlock = false;
	      for( long i = outline.data.size() - 1; i > 0; i-- )
	      {
		  if( (outline.data.at(i).get(0) == newVec.get(0)) &&
		      (outline.data.at(i).get(1) == newVec.get(1)) )
		  {
		      visited.at(i)++;
		      if( visited.at(i) >= 3 )
			  deadlock = true;

		      found = true;
		      break;
		  }
	      }

	      if( deadlock )
		  break;
	      
	      //Add to results if not found already
	      if( !found )
	      {
		  //Else mark new P 
		  resultData[ srcPos ] = 1.0;
	    
		  outline.data.push_back(newVec);
		  visited.push_back(1);
	      }

	      //Try to find start for circle
	      for(circleStart = 0; circleStart < 8; circleStart++)
		  if(( qPos.vec[0] == pPos.vec[0] +  (circleOffset.vec[circleStart]).val.first) &&
		     ( qPos.vec[1] == pPos.vec[1] +  (circleOffset.vec[circleStart]).val.second) )
		      break;

	      //Try to find new P and Q
	      for( int q = circleStart ; q < circleStart + 8; q++)
	      {
		  srcPos = (pPos.vec[1] + (circleOffset.vec[q % 8]).val.second) * ((segments->at(currRegion))->getDimension()).vec[0] 
			  + pPos.vec[0] + (circleOffset.vec[q % 8]).val.first;

		  if( (double)regionData[ srcPos ] == 0.0 )
		  {
		      //New P and Q found
		      pPos.vec[0] +=  (circleOffset.vec[q % 8]).val.first;
		      pPos.vec[1] +=  (circleOffset.vec[q % 8]).val.second;      
		      qPos.vec[0] = pPos.vec[0] + (circleOffset.vec[(q-1) % 8]).val.first;
		      qPos.vec[1] = pPos.vec[1] + (circleOffset.vec[(q-1) % 8]).val.second;
		      break;
		  }
	      }
	  }

	  //Finally delete original image and set to outline array
	  segments->at(currRegion) = resultSegment;
	  delete currSegment;
      }

      //Set bitvector
      computedProps |= Outline;
      return;
  }

  /** computes samples outline (silhouette edge) of all 2D segments in regions if padded at least 1 pixel 
      sampledistance is the distance on the outline between the distinct samples */
  template<class T>
  vector<SampleVector*>* ConnectedComponentProperties<T>::extractOutlineSamples( unsigned int sampleDist, std::vector<unsigned int>* regions  )
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Outline shall be computed "
		    << "but need segments to be computed in advance." << std::endl;
	  return NULL;
      }

      if( srcDim.vec.size() != 2 )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Oultline shall be computed "
		    << "but can only be computed on 2D images." << std::endl;
	  return NULL;
      }
	
      if( regions == NULL || regions->size() == 0 )
	  return NULL;

      //Create result array for the samples
      vector<SampleVector*>* result = new vector<SampleVector*>();

      //Fill result with empty pointers
      result->clear();
      for( unsigned int i = 0; i < regions->size(); i++ )
	  result->push_back(NULL);
    
      //Set up circle offset vector
      IntRangeList circleOffset;
      circleOffset.vec.push_back(IntRange(-1, 0));
      circleOffset.vec.push_back(IntRange(-1, 1));
      circleOffset.vec.push_back(IntRange( 0, 1));
      circleOffset.vec.push_back(IntRange( 1, 1));
      circleOffset.vec.push_back(IntRange( 1, 0));
      circleOffset.vec.push_back(IntRange( 1,-1));
      circleOffset.vec.push_back(IntRange( 0,-1));
      circleOffset.vec.push_back(IntRange(-1,-1));

      CoordinateVector pPos(2);
      CoordinateVector qPos(2);
      CoordinateVector startPos(2);
      int circleStart;

      // Iterate over all regions
      unsigned int srcPos;

      Array<T>* currSegment;
      T* regionData;
      for( long currRegion = 0 ; currRegion < regions->size(); currRegion++ )
      {
	  if( segments->at(regions->at(currRegion)) == NULL )
	      continue;

	  // No outline if padding not enabled
	  if( padding->at(regions->at(currRegion)) < 1 )
	  {	
	      std::cerr << "Could not compute outline of region " << currRegion
			<< ". No padding was done." << std::endl;
	      continue;
	  }

	  //Create sample vector for current region
	  result->at(currRegion) = new SampleVector();

	  //Fetch regiondata
	  currSegment = segments->at(regions->at(currRegion));
	  regionData = &((*((*currSegment)[0]))[0]);
	  
	  //First, try to find an edge constellation
	  qPos.vec[0] = 0;
	  qPos.vec[1] = 0;
	  for( unsigned long j = 0; j < regionEntries->at(regions->at(currRegion)); j++ )
	  {
	      if( (double)(regionData[j]) > 0.5 )
		  break;

	      // update position in region array
	      qPos.vec[0]++;
	      if( qPos.vec[0] >= ((segments->at(regions->at(currRegion)))->getDimension()).vec[0] )
	      {
		  qPos.vec[1]++;
		  qPos.vec[0] = 0;
	      }
	  }
      
	  //First p is now on the left of qPos as padding is at least one pixel
	  pPos.vec[0] = qPos.vec[0] - 1;
	  pPos.vec[1] = qPos.vec[1];

	  //Save start position
	  startPos.vec[0] = qPos.vec[0];
	  startPos.vec[1] = qPos.vec[1];

	  //Keep track of samples for deadlock detection
	  Vector newVec( 2, 0.0 ), zeroVec( 2, 0.0 );
	  SampleVector outline;
	  //Vector of number of visits of outline points
	  vector<long> visited;

	  //Repeat until outline completed or until upper limit is achieved
	  unsigned long sample = 0;
	  for( unsigned long lineLength = 0 ; lineLength < ULONG_MAX; lineLength++ )
	  {
      
	      //End if line closed or deadlock
	      if( lineLength >= 1 && startPos.vec[0] == qPos.vec[0] && startPos.vec[1] == qPos.vec[1] )
		  break;

	      //Position in start array
	      srcPos = pPos.vec[1] * ((segments->at(regions->at(currRegion)))->getDimension()).vec[0] + pPos.vec[0];

	      //Create new sample for outline
	      newVec.copy( zeroVec );
	      newVec.set( 0, (double)(pPos.vec[0])); newVec.set( 1, (double)(pPos.vec[1]));

	      //Try to trace back already found outline samples
	      bool found = false;
	      bool deadlock = false;
	      for( long i = outline.data.size() - 1; i > 0; i-- )
	      {
		  if( (outline.data.at(i).get(0) == newVec.get(0)) &&
		      (outline.data.at(i).get(1) == newVec.get(1)) )
		  {
		      visited.at(i)++;
		      if( visited.at(i) >= 3 )
			  deadlock = true;

		      found = true;
		      break;
		  }
	      }

	      if( deadlock )
		  break;
	      
	      //Add to results if not found already
	      if( !found )
	      {
		  //Else mark new P 
		  if( sampleDist == 0 || sample % (sampleDist + 1) == 0 )
		  {
		      ((result->at(currRegion))->data).push_back( newVec );
		      //Reset newVec for outline re-traversing
		      newVec.copy( zeroVec );
		      newVec.set( 0, (double)(pPos.vec[0])); newVec.set( 1, (double)(pPos.vec[1]));
		      sample=1;
		  }else
		      sample++;
		  
	    
		  outline.data.push_back(newVec);
		  visited.push_back(1);
	      }

	      //Try to find start for circle
	      for(circleStart = 0; circleStart < 8; circleStart++)
		  if(( qPos.vec[0] == pPos.vec[0] +  (circleOffset.vec[circleStart]).val.first) &&
		     ( qPos.vec[1] == pPos.vec[1] +  (circleOffset.vec[circleStart]).val.second) )
		      break;

	      //Try to find new P and Q
	      for( int q = circleStart ; q < circleStart + 8; q++)
	      {
		  srcPos = (pPos.vec[1] + (circleOffset.vec[q % 8]).val.second) * ((segments->at(regions->at(currRegion)))->getDimension()).vec[0] 
			  + pPos.vec[0] + (circleOffset.vec[q % 8]).val.first;

		  if( (double)regionData[ srcPos ] == 0.0 )
		  {
		      //New P and Q found
		      pPos.vec[0] +=  (circleOffset.vec[q % 8]).val.first;
		      pPos.vec[1] +=  (circleOffset.vec[q % 8]).val.second;      
 		      qPos.vec[0] = pPos.vec[0] + (circleOffset.vec[(q+8-1) % 8]).val.first;
 		      qPos.vec[1] = pPos.vec[1] + (circleOffset.vec[(q+8-1) % 8]).val.second;
		      break;
		  }
	      }
	  }

      }

      return result;
  }


  /** returns the boundingbox of the defined region*/
  template<class T>
  LongRangeList ConnectedComponentProperties<T>::getBoundingbox( unsigned int region )
  {
      LongRangeList result( srcDim.vec.size() );
      if( ! getPropertyComputed( BoundingBox) || region <= 0 || region > numLabels )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Boundingbox is queried "
		    << "but has not been computed or queried region out of range." << std::endl;
	  return result;
      }

      for( int j = 0; j < srcDim.vec.size(); j++)
	  result.vec[j] = LongRange( ((bboxes->at(region))->vec[j]).val.first, ((bboxes->at(region))->vec[j]).val.second );
      
      return result;
  }

  /** returns the area of the defined region*/
  template<class T>
  unsigned long ConnectedComponentProperties<T>::getArea( unsigned int region )
  {
      if( ! getPropertyComputed( Area ) || region <= 0 || region > numLabels )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Area is queried "
		    << "but has not been computed or queried region out of range." << std::endl;
	  return 0;
      }

      return area->at( region );
  }

  /** returns the extracted segment of the defined region*/
  template<class T>
  Array<T>* ConnectedComponentProperties<T>::getSegment( unsigned int region )
  {
      if( ! getPropertyComputed( Segments ) || region <= 0 || region > numLabels )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Segment is queried "
		    << "but has not been computed or queried region out of range." << std::endl;
	  return NULL;
      }

      if( segments->at(region) == NULL )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Segment is queried "
		    << "but is not valid" << std::endl;
	  return NULL;
      }

      return segments->at(region);
  }

  /** returns the symmetrical padding of segment of the defined region*/
  template<class T>
  unsigned int ConnectedComponentProperties<T>::getPadding( unsigned int region )
  {
      if( ! getPropertyComputed( Segments ) || region <= 0 || region > numLabels )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Padding is queried "
		    << "but segments have not been computed or queried region out of range." << std::endl;
	  return 0;
      }

      return padding->at(region);
  }

  /** returns the eulernumber of the defined region*/
  template<class T>
  double ConnectedComponentProperties<T>::getEulerNumber( unsigned int region )
  {
      if( ! getPropertyComputed( EulerNumber ) || region <= 0 || region > numLabels )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Eulernumber is queried "
		    << "but has not been computed or queried region out of range." << std::endl;
	  return 0.0;
      }

      if( segments->at(region) == NULL )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Eulernumber is queried "
		    << "but is not valid" << std::endl;
	  return 0.0;
      }

      return eulerNum->at(region);
  }

  /** returns if the segment is valid or not*/
  template<class T>
  bool ConnectedComponentProperties<T>::isValid( unsigned int region )
  {
      if( ! getPropertyComputed( Segments ) || region <= 0 || region > numLabels  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: Validation of segment " << region <<" is queried "
		    << "but segments hasve not been computed or queried region out of range." << std::endl;
	  return false;
      }

      return ((segments->at(region)) != NULL);
  }

  /** returns the number of segments*/
  template<class T>
  unsigned int ConnectedComponentProperties<T>::getNumSegments()
  {
      if( ! getPropertyComputed( Segments )  )
      {
	  std::cerr << "Error in ConnectedComponentProperties: NumSegments is queried "
		    << "but has not been computed." << std::endl;
	  return 0;
      }
	
      unsigned int numValidSegments = 0;
      for( unsigned int i = 1; i <  segments->size(); i++)
	  if( segments->at(i) != NULL )
	      numValidSegments++;

      return numValidSegments;
  }

// template instantiation code
template class ConnectedComponentProperties<float>;
template class ConnectedComponentProperties<double>;

} /* namespace */

#endif /* CONNCOMPPROPS_CONNECTEDCOMPONENTPROPERTIES_C */

