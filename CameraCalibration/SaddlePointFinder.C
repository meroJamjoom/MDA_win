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

#ifndef CONNCOMPPROPS_SADDLEPOINTFINDER_C
#define CONNCOMPPROPS_SADDLEPOINTFINDER_C

#include "SaddlePointFinder.hh"
#include <math.h>

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;


  /** default constructor, mask size is wintx * winty, a  central zero zone of the size  wx2, wy2 can be added */
  template<class T>
  SaddlePointFinder<T>::SaddlePointFinder( Array<T> *image, int wintx, int winty, int wx2, int wy2)
  {

      this->image = image;

      // Get size and dimensions of Array
      imgDim = image->getDimension();
      if( imgDim.vec.size() != 2 )
      {
	  cerr << " Error: SaddlePointFinder only supports 2D input images." << endl;
	  return;
      }

      //Create mask
      this->wintx = wintx;
      this->winty = winty;
      mask = new double[ (wintx * 2 + 1) * (winty * 2 + 1) ];

      //Fill mask
      int index = 0;
      for( double y = -winty; y <=winty; y+=1.0 )
      {
	for( double x = -wintx; x <=wintx; x+=1.0 )
	{
	    mask[index] =  exp( - pow( x / ((double) wintx) , 2.0) ) * exp( - pow( y / ((double) winty) , 2.0 ) ); 
	    index++;
	}
      }

      if( (wx2>0) && (wy2>0) )
      {
	if( ((wintx - wx2)>=2)&&((winty - wy2)>=2) )
	{
	    for( int y = winty-wy2; y <= winty+wy2; y++ )
	    {
		for( int x = wintx-wx2; x <= wintx+wx2; x++ )
		{
		    mask[y* ( 2 * wintx + 1 ) + x] = 0.0;
		}
	    }
	}
      }

      //Set algorithms default parameters
      resolution = 0.005;
      MaxIter = 10;

      return;
  }

  /** default destructor */
  template<class T>
  SaddlePointFinder<T>::~SaddlePointFinder()
  {
      delete[] mask;
  }

  /** finds the saddlepoints in the image defined in the constructor from initial guess positions
      the results are returned in the guesses SampleVector. The boolean vector of the same size as SampleVector
      holds information about which points diverged */
  template<class T>
  void SaddlePointFinder<T>::findSaddles( SampleVector *guesses,  vector<bool>* divergedGuesses )
  {
      //Check for dimensions
      if( imgDim.vec.size() != 2 )
      {
	  cerr << " Error: SaddlePointFinder only supports 2D input images." << endl;
	  return;
      }
    
      //Check for sanity of input parameters
      if( guesses == NULL || divergedGuesses == NULL )
	  return;

      //Initialize divergedGuesses Vector
      divergedGuesses->clear();
      for( unsigned long i=0; i < (guesses->data).size(); i++ )
	  divergedGuesses->push_back(false);
  
      //Initialize counter variables for loop
      unsigned long nx  = imgDim.vec[0];
      unsigned long ny = imgDim.vec[1]; 
      double v_extra;
      long compt, xmin, xmax, ymin, ymax;
      double cIx, cIy, crIx, crIy;
      Vector originalGuess( 2, 0.0 );
      //Pointer to original image data
      T* srcData = &((*((*image)[0]))[0]);
      //Extracted and interpolated region
      CoordinateVector extractedSamplesSize;
#if defined (_WIN32) || defined (_WIN64)
	  Array<T>* extractedSamples = NULL;
	  T* extractedSampleData = NULL;
#else
      Array<T>* extractedSamples;
      T* extractedSampleData;
#endif
      //Iterate over all saddles
      for( unsigned long i=0; i < (guesses->data).size(); i++ )
      {
	v_extra = resolution + 1.0; 			// just larger than resolution
	
	compt = 0; 					// no iteration yet

	originalGuess.assign( (guesses->data).at(i) );	// save original guess for divergence check
	
	// #### Start gradient zero finding iteration ####
	bool zeroFindingIteration = true;
	while( zeroFindingIteration )
	{
	    // End iteration if v_extra under the desired resolution or max resolution achieved
	    // After last iteration only extract and interpolate region and do final averaging
	    if( !( (v_extra > resolution) && (compt < MaxIter) ) )
		zeroFindingIteration = false;

	    //Samplepoint position
	    cIx = ((guesses->data).at(i)).get( 0 );
	    cIy = ((guesses->data).at(i)).get( 1 ); 
	    // on the initial image
	    crIx = round(cIx);
	    crIy = round(cIy);  
	    
	    // What if the sub image is not in?
	    if( crIx-wintx-2 < 0)
	    {
		xmin = 0; 
		xmax = 2 * wintx + 4;
	    }
	    else if( crIx+wintx+2 >= nx)
	    {
		xmax = (nx - 1); 
		xmin = (nx - 1) - 2*wintx - 4;
	    }
	    else
	    {
		xmin = crIx - wintx - 2;
		xmax = crIx + wintx + 2;
	    }
	    
	    if( crIy-winty-2 < 0)
	    {
		ymin = 0; 
		ymax = 2 * winty + 4;
	    }
	    else if(crIy+winty+2 >= ny)
	    {
		ymax = (ny - 1);
		ymin = (ny - 1) - 2*winty - 4;
	    }
	    else
	    {
		ymin = crIy - winty - 2;
		ymax = crIy + winty + 2;
	    }

	    
	    // Get necessary neighborhood
	    extractedSamplesSize.vec.clear();
	    extractedSamplesSize.vec.push_back( xmax - xmin + 1 );
	    extractedSamplesSize.vec.push_back( ymax - ymin + 1 );
	    extractedSamples = new Array<T>( extractedSamplesSize );
	    extractedSamples->addChannel();
	    extractedSampleData = &((*((*extractedSamples)[0]))[0]);
	    for( long y = 0; y < extractedSamplesSize.vec[1]; y++ )
	    {
		for( long x = 0; x < extractedSamplesSize.vec[0]; x++ )
		{
		    extractedSampleData[ y * extractedSamplesSize.vec[0] + x ] = srcData[ (y + ymin)*nx + x + xmin ];
		}
	    }

	    //Interpolate on subpixel level

	    // Coefficients to compute the sub pixel accuracy.
	    double itIx = cIx - crIx;
	    double itIy = cIy - crIy;
	    double* vIx = new double[3];
	    double* vIy = new double[3];

	    if( itIx > 0.0 )
	    {
		vIx[0] = 0.0;
		vIx[1] = 1.0 - itIx;
		vIx[2] = itIx;
	    }
	    else
	    {
		vIx[0] = -itIx;
		vIx[1] = 1.0 + itIx;
		vIx[2] = 0.0 ;
	    }

	    if( itIy > 0.0)
	    {
		vIy[0] = 0.0;
		vIy[1] = 1.0 - itIy;
		vIy[2] = itIy;
	    }
	    else
	    {
		vIy[0] = -itIy;
		vIy[1] = 1.0 + itIy;
		vIy[2] = 0.0 ;
	    }

	    //X-Interpolation
	    Linear1DFilter<T> interpolateX(1);
	    interpolateX.setFilter(1, vIx);
	    BoundaryMethod boundary= Clamp;
	    unsigned long interpolateXChannel= extractedSamples->addChannel();
	    interpolateX.apply( *extractedSamples, boundary, 0, 0, interpolateXChannel );

	    //Y-Interpolation
	    Linear1DFilter<T> interpolateY(1);
	    interpolateY.setFilter(1, vIy);
	    interpolateY.apply( *extractedSamples, boundary, 1, interpolateXChannel, 0 );

	    //Cleanup
	    delete[] vIx;
	    delete[] vIy;

	    //Exit iteration to final averaging
	    if( !zeroFindingIteration )
		break;

	    //Take only subarea of that
	    CoordinateVector neighborSamplesSize;
	    neighborSamplesSize.vec.push_back( 2*wintx+2 + 1 );
	    neighborSamplesSize.vec.push_back( 2*winty+2 + 1 );
	    Array<T>* neighborSamples = new Array<T>( neighborSamplesSize );
	    neighborSamples->addChannel();
	    T* sampleData = &((*((*neighborSamples)[0]))[0]);
	    for( long x = 1; x < 2*wintx+4; x++ )
	    {
		for( long y = 1; y < 2*winty+4; y++ )
		{
		    sampleData[ (y-1) * neighborSamplesSize.vec[0] + (x - 1) ] = extractedSampleData[ y * extractedSamplesSize.vec[0] + x ];
		}
	    }

	    //Cleanup
	    delete extractedSamples;

	    //Calculate gradient of the extracted neighbor samples
	    FirstDerivative1D<T> derivative;
	    // need two temporary channels:
	    unsigned long derivX= neighborSamples->addChannel();
	    unsigned long derivY= neighborSamples->addChannel();

	    // Compute derivative along x and y axis
	    derivative.apply( *neighborSamples, boundary, 0, 0, derivX );
	    derivative.apply( *neighborSamples, boundary, 1, 0, derivY );

	    // Extraction of the useful parts only of the gradients
	    CoordinateVector gradientImageSize;
	    gradientImageSize.vec.push_back( 2*wintx + 1 );
	    gradientImageSize.vec.push_back( 2*winty + 1 );

	    T* gradientX = new T[ (gradientImageSize.vec[0]) * (gradientImageSize.vec[1]) ];
	    T* gradientY = new T[ (gradientImageSize.vec[0]) * (gradientImageSize.vec[1]) ];
	    T* sampleDerivX = &((*((*neighborSamples)[derivX]))[0]);
	    T* sampleDerivY = &((*((*neighborSamples)[derivY]))[0]);
	    for( long y = 1; y < 2*winty+2; y++ )
	    {
		for( long x = 1; x < 2*wintx+2; x++ )
		{
		    gradientX[ (y-1) * gradientImageSize.vec[0] + x - 1 ] = (sampleDerivX[ y * neighborSamplesSize.vec[0] + x])/2.0;
		    gradientY[ (y-1) * gradientImageSize.vec[0] + x - 1 ] = (sampleDerivY[ y * neighborSamplesSize.vec[0] + x])/2.0;
		}
	    }

	    //Cleanup
	    delete neighborSamples;


	    //Find new gradient zero crossing
	    Vector newSaddle( 2, 0.0);
	    Vector tempSaddle( 2, 0.0);
	    T gxx, gyy, gxy;
	    unsigned long currIndex;
	    T a = 0.0;
	    T b = 0.0;
	    T c = 0.0;
	    for( long y = 0; y < gradientImageSize.vec[1]; y++ )
	    {
		for( long x = 0; x < gradientImageSize.vec[0]; x++ )
		{
		    currIndex = y * gradientImageSize.vec[0] + x;
		    gxx = (gradientX[currIndex]) * (gradientX[currIndex]) * mask[currIndex];
		    gyy = (gradientY[currIndex]) * (gradientY[currIndex]) * mask[currIndex];
		    gxy = (gradientX[currIndex]) * (gradientY[currIndex]) * mask[currIndex];
		    tempSaddle[0] += gxx * ( (T)x - (T)wintx + cIx ) + gxy * ( (T)y - (T)winty + cIy );
		    tempSaddle[1] += gxy * ( (T)x - (T)wintx + cIx ) + gyy * ( (T)y - (T)winty + cIy );

		    a += gxx;
		    b += gxy;
		    c += gyy;
		}
	    }
	    //Cleanup
	    delete[] gradientX;
	    delete[] gradientY;
	    
	    T dt = a*c - b*b;
          
	    newSaddle[0] = c * (tempSaddle[0]) - b * (tempSaddle[1]);
	    newSaddle[1] = a * (tempSaddle[1]) - b * (tempSaddle[0]);
	    newSaddle /= dt;

	    //Shifting vector length
	    v_extra = sqrt( pow( newSaddle[0] - cIx  , 2.0) + pow( newSaddle[1] - cIy , 2.0) );

	    //Assign new saddle position
	    ((guesses->data).at(i)).assign( newSaddle );

	    // Next zero crossing finding iteration
	    compt++;

	 } //End of zero crossing finding iteration


	  //#### Final averaging step, extractedSamples are extracted and interpolated  ####

	  //Compute Gauss-weighted image vector
	  Vector gaussVec( ( 2*wintx + 1 ) * ( 2*winty + 1 ), 0.0 );
	  long vecIndex = 0;
	  for( long x = 2; x < 2*wintx+2 + 1; x++ )
	  {
	      for( long y = 2; y < 2*winty+2 + 1; y++ )
	      {
		  gaussVec.set( vecIndex, (extractedSampleData[ y * extractedSamplesSize.vec[0] + x ]) * (mask[vecIndex])  );
		  vecIndex++;
	      }
	  }

	  //Cleanup
	  delete extractedSamples;


	  //Assemble new Matrix A for averaging
	  Matrix A( (2*wintx+1)*(2*winty+1), 6 );
	  double value;
	  double px;
	  double py;
	  for( long col = 0; col < A.getNumColumns(); col++ )
	  {
	      for( long row = 0; row < A.getNumRows(); row++ )
	      {
		  value = mask[row];
		  px = (double)( row % (2*wintx+1) ) - (double)wintx + cIx;
		  py = (double)( (row - row % (2*wintx+1)) / (2*wintx+1) ) - (double)winty + cIy;

		  switch( col )
		  {
		      case 0:
			    value *= px*px;
			    break;
		      case 1:
			    value *= px*py;
			    break;
		      case 2:
			    value *= py*py;
			    break;
		      case 3:
			    value *= px;
			    break;
		      case 4:
			    value *= py;
			    break;
		 }

		  A.set( row, col, value );
	      }
	  }

	  //Average and extract final point coordinates
	  Matrix A_transpose = A.getTranspose();
	  Matrix AA( 6, 6 );
	  multMatrixMatrix( A_transpose, A, AA );
	  Matrix AA_inv( 6, 6 );
	  AA_inv = inverse( AA );
	  Matrix AAA( 6, (2*wintx+1)*(2*winty+1) );
	  multMatrixMatrix( AA_inv, A_transpose, AAA );

	  //Multiply with gaussvec
	  Vector pointVec(6, 0.0);
	  AAA.rightMultiply( gaussVec, pointVec );
	  Matrix pointMat(2, 2 );
	  pointMat.set( 0, 0, 2.0 * pointVec.get(0) );
	  pointMat.set( 0, 1, pointVec.get(1) );
	  pointMat.set( 1, 0, pointVec.get(1) );
	  pointMat.set( 1, 1, 2.0 * pointVec.get(2) );
	  Matrix pointMat_inv(2, 2 );
	  pointMat_inv = inverse( pointMat );
	  pointMat_inv *= -1.0;
	  Vector pointCoeff( 2, 0.0);
	  pointCoeff.set( 0, pointVec.get( 3 ) );
	  pointCoeff.set( 1, pointVec.get( 4 ) );
	  Vector finalPoint( 2, 0.0 );
	  pointMat_inv.rightMultiply( pointCoeff, finalPoint );
	  
	  //Assign final vector
	  ((guesses->data).at(i)).assign( finalPoint );

	  // #### Check for points that diverge:
	  if( fabs( originalGuess.get(0) - finalPoint.get(0) ) > wintx || fabs( originalGuess.get(1) - finalPoint.get(1)) > winty )
	  {
		// For the diverged points, keep the original guesses
		((guesses->data).at(i)).assign( originalGuess ); 
		divergedGuesses->at(i) = true;
	  }


      }  //End of for-loop over all saddles

      return;
  }

  // template instantiation code
  template class SaddlePointFinder<float>;
  template class SaddlePointFinder<double>;

} /* namespace */

#endif /* CONNCOMPPROPS_SADDLEPOINTFINDER_C */

 
