// ==========================================================================
// $Id: CALTagPattern.C 757 2010-09-22 21:47:34Z bradleyd $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2010-, UBC
//
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef CAMERACALIBRATION_CALTAGPATTERN_C
#define CAMERACALIBRATION_CALTAGPATTERN_C

#include <time.h>
#include <fstream>

#include <MDA/Base/BitsAndBytes.hh>
#include <MDA/Base/CRCCode.hh>

#include "CALTagPattern.hh"
#include "CALTagPostscript.hh"

#if defined (_WIN32)|| defined(_WIN64)
#define srand48(num) srand(num)
#define lrand48() rand() 
#endif
/////////////////////// includes from Felix's code ////////////////////////////////////////////////////////
#include <algorithm>
#include <Magick++.h>
using namespace Magick;
///////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace MDA
{

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  /** read configuration from MetaData */
  bool
  CALTagPattern::configure( MetaData &md, const char *path ) {
    string p( path );

    bool success= true;

    // layout defaults to 1
    layout= 1;
    get( md, (p+":layout").c_str(), layout );

    // scale defaults to 1.0
    scale= 1.0;
    get( md, (p+":scale").c_str(), scale );

    // width and height are mandatory...
    rows= columns= 0u;
    success&= warnCond( get( md, (p+":rows").c_str(), rows ),
      "no rows attribute" );
    success&= warnCond( get( md, (p+":columns").c_str(), columns ),
      "no cols attribute" );

    // number of IDs has to at least match number of markers (extra ones
    // are ignored in the pipeline)
    success&= get( md, path, ids );
    warnCond( rows> 0 && columns> 0 && ids.vec.size()>= rows*columns,
      "number of marker IDs does not match grid dimensions" );

    return success;
  }

  /** write configuration to MetaData (pure virtual ) */
  void
  CALTagPattern::writeConfig( MetaData &md, const char *path ) {
    string p( path );

    // write attributes
    set( md, (p+":layout").c_str(), layout );
    set( md, (p+":scale").c_str(), scale );
    set( md, (p+":rows").c_str(), rows );
    set( md, (p+":columns").c_str(), columns );

    // write IDs
    set( md, path, ids );
  }

  /** change layout (0 or less for any parameter means to keep the current value)
      \param r: new # rows (0 means keep current)
      \param c: new # columns (0 means keep current)
      \param s: new scale factor (<= 0 means keep current)
      \param l: new layout (0 means keep current)
  */
  void
  CALTagPattern::changeLayout( unsigned r, unsigned c, double s, unsigned l ) {
    if( r> 0 )
      rows= r;
    if( c> 0 )
      columns= c;
    if( s> 0.0 )
      scale= s;
    if( l> 0 )
      layout= l;
  }

  /** generate the ID list for the current geometric layout, either
      from scratch, or using existing IDs for the first markers */
  bool
  CALTagPattern::makeIDs( bool fromScratch ) {
    unsigned i, j;

    // if we can use enough existing IDs, we are done
    if( !fromScratch && ids.vec.size()>= rows*columns )
      return true;

    // generate code table if not already existing
    if( codes== NULL || valid== NULL )
      computeCodes();

    // skip first "firstIndex-1" valid markers
    unsigned count= 0;
    unsigned numEntries= 1 << (bits-cBits);
    for( i= 0 ; i< numEntries && count< firstIndex ; i++ )
      if( valid[i] )
        count++;

    // determine how many IDs are missing, and copy them from the code table
    unsigned missing;
    if( fromScratch ) {
      ids.vec.clear();
      missing= rows*columns;
    }
    else
      missing= rows*columns - ids.vec.size();
    for( count= 0 ; count< missing && i< numEntries; i++ ) {
      if( valid[i] ) {
        ids.vec.push_back( codes[i] );
        count++;
      //	cerr << "appending code[i] == " << codes[i] << endl;
      //	CRCCode codec( idBits, cBits, genPoly );
      //	cerr << "   encoded would be: " << codec.encode(codes[i]) << endl;
      }
    }

    // check if we got enough codes...
    return warnCond( count== missing, "not enough valid codes to fill grid..." );
  }

  /** set IDs to specified list (either with or without CRC portion)  */
  bool
  CALTagPattern::setIDs( CoordinateVector &_ids, bool makeChecksum ) {
    // copy IDs
    ids= _ids;

    // create checksum if desired
    if( makeChecksum ) {
      errorCond( cBits< bits, "too many bits in CRC" );
      unsigned dBits= bits-cBits;
      CRCCode codec( dBits, cBits, genPoly );

      for( unsigned i= 0 ; i< ids.vec.size() ; i++ )
        ids.vec[i]= codec.encode( ids.vec[i] );
    }

    return warnCond( ids.vec.size()>= rows*columns,
      "insufficient number of marker IDs" );
  }

  //[brad] moved below
  ///** run the detection algorithm for the grid */
  //bool
  //CALTagPattern::detect()
  //{
  //	//!! TODO
  //}

  /** create a PostScript file with the calibration pattern */
  bool
  CALTagPattern::makePattern( const char *fileName, const PaperSize &paper ) {
    ofstream psFile( fileName );

    // get time (for creation date)
    time_t rawtime;
    time( &rawtime );

    // write postscript header with correct orientation
    bool landscape= rows<columns;
    psFile << "%!PS-Adobe-2.0" << endl
      << "%%Title: CALTag calibration pattern" << endl
      << "%%Pages: 1" << endl
      << "%%BoundingBox: 0 0 "
      << paper.width*72 << ' ' << paper.height*72 << endl
      << "%%Orientation: " << (landscape?"Landscape":"Portrait") << endl
      << "%%Creator: MDA tools - CALTagPattern class" << endl
      << "%%CreationDate: " << ctime( &rawtime )
      << "%%EndComments" << endl << endl
      << "%--------------------- config block ------------------" << endl
      << endl;

    // write configuration section
    psFile << "/layout { " << layout << " } bind def" << endl << endl
      << "/spacing { " << scale << " } bind def" << endl << endl
      << "/columns { " << columns << " } bind def" << endl << endl
      << "/rows { " << rows << " } bind def" << endl << endl
      << "/drawcropmarkers { true } bind def" << endl << endl
      << "/markerids { " << endl << "  [";
    for( unsigned i= 0 ; i< rows ; i++ ) {
      for( unsigned j= 0 ; j< columns ; j++ )
        psFile << ' ' << ids.vec[i*columns+j];
      if( i< rows-1 )
        psFile << endl << "   ";
    }
    psFile << " ]" << endl << "} bind def" << endl << endl
      << "%------------------- end config block ----------------" << endl
      << endl;

    // in landscape mode, we need to offset and rotate the coordinate frame
    if( landscape )
      psFile << endl << "% rotation to landscape mode" << endl
        << paper.width*72 << " 0 translate" << endl
        << "90 rotate" << endl;

    // write procedural elements for generating the pattern
    psFile << CALTAG_PS_CODE << endl;

    psFile.close();

    // see if this all fits onto the page...
    return warnCond( scale*((landscape?rows:columns)+1.5)<= paper.width &&
      scale*((landscape?columns:rows)+1.5)<= paper.height,
      "pattern exceeding paper dimensions - output clipped" );
  }

  /* purely local helper function for marking codes that are too close */
  void
  markTooClose( unsigned id, unsigned numElem, unsigned *codes, bool *valid ) {
    unsigned i;
    for( i= 0 ; i< numElem ; i++ )
      valid[i]&= hammingDist( id, codes[i] )> 1;
  }

  /** precompute the full code array and which codes are valid */
  void
  CALTagPattern::computeCodes() {
    unsigned i, j;

    // set up CRC codec
    errorCond( cBits< bits, "too many bits in CRC" );
    unsigned dBits= bits-cBits;
    CRCCode codec( dBits, cBits, genPoly );
    unsigned numEntries= 1 << dBits;

    // re-allocate memory
    // (future proof: in case the split between cBits and dBits has changed)
    if( codes )
      delete [] codes;
    codes= new unsigned[numEntries];
    if( valid )
      delete [] valid;
    valid= new bool[numEntries];

    // generate the code array and initial validity (based on bit counts)
    for( i= 0 ; i< numEntries ; i++ ) {
//cerr << "precomputing code:  i == " << i << ",  crc == " << codec.encode(i) << endl;
      codes[i]= codec.encode( i );
      j= hammingDist( codes[i] );
      valid[i]= j>= bits/4 && j<= 3*bits/4;
    }

    // now update validity based on rotational symmetry constraints
    for( i= 0 ; i< numEntries ; i++ )
    if( valid[i] ) {
      j= rotate( codes[i] );
      markTooClose( j, numEntries, codes, valid );
      j= rotate( j );
      markTooClose( j, numEntries, codes, valid );
      j= rotate( j );
      markTooClose( j, numEntries, codes, valid );
    }
  }

  /** code rotation - hard coded to 4x4 configurations right now */
  unsigned
  CALTagPattern::rotate( unsigned c ) {
    return
      ((c&0x0001)?0x1000:0) | ((c&0x0002)?0x0100:0) | ((c&0x0004)?0x0010:0) | ((c&0x0008)?0x0001:0) |
      ((c&0x0010)?0x2000:0) | ((c&0x0020)?0x0200:0) | ((c&0x0040)?0x0020:0) | ((c&0x0080)?0x0002:0) |
      ((c&0x0100)?0x4000:0) | ((c&0x0200)?0x0400:0) | ((c&0x0400)?0x0040:0) | ((c&0x0800)?0x0004:0) |
      ((c&0x1000)?0x8000:0) | ((c&0x2000)?0x0800:0) | ((c&0x4000)?0x0080:0) | ((c&0x8000)?0x0008:0);
  }

  //////////// Added by Brad /////////////////////////////////////////////////////////////////////////////////////

  void CALTagPattern::constructCorrespondences() {
    if( !successfulTracked ) {
      cerr << "Must run detect() before constructing correspondences" << endl;
      return;
    }
    int index;
    for( int row=0; row<rows+1; row++ ) {
      for( int col=0; col<columns+1; col++ ) {
        index = row*(columns+1) + col;
        double img0 = (imageCalibPoints->data).at(index).get(0);
        double img1 = (imageCalibPoints->data).at(index).get(1);
        if( index<(imageCalibPoints->data).size()  &&  img0>=0.0  && img0<arrayDim.vec[0] &&  img1>=0.0 && img1<arrayDim.vec[1]) {
          Vector iPt( 2, 0.0 );
          iPt[0] = img0;
          iPt[1] = img1;
          //[brad] possibly have to transpose?
          Vector wPt( 3, 0.0 );
          wPt[0] = row*scale;
          wPt[1] = col*scale;
          wPt[2] = 0.0;
          PointCorrespondence pc( iPt, wPt );
          correspondences.push_back( pc );
        }
      }
    }
  }

  void CALTagPattern::writeCorrespondences( MetaData& md, const char* path, const char* dest ) {
    if( getCorrespondences().size() <= 0 ) {
      cerr << "Must run constructCorrespondences() before writing them" << endl;
      return;
    }
    errorCond( path!=NULL, "Invalid file specified" );
    string p(path);
    set( md, (p+".corr:width").c_str(), arrayDim.vec[0] );
    set( md, (p+".corr:height").c_str(), arrayDim.vec[1] );

    // now write the correspondences
    set( md, path, correspondences);
    
    md.write( dest );
  }

  //////////// Taken from Felix's code ///////////////////////////////////////////////////////////////////////////

  // comparison function for the indexed sorting
  bool CALTagPattern::compareIndexPairs ( Vector i, Vector j) { return ( i[0] < j[0] ); }

  /** default constructor */
  CALTagPattern::CodeBlobProps::CodeBlobProps() : orientation( Vector(2,0.0) ) {
    ID = -1;
    valid = false;
    firstCorner = 0;

    corners.data.clear();

    return;
  }

  /** default destructor */
  CALTagPattern::CodeBlobProps::~CodeBlobProps() {
    corners.data.clear();
    return;
  }

  void CALTagPattern::setImage( Array<T> *array ) {
    // Early exit if NULL
    if( array == NULL )
      return;

    //Prepare the image
    this->array = array;

    arrayDim = array->getDimension();

    // Early exit if wrong dimension
    if( arrayDim.vec.size() != 2 ) {
      cerr << ": Supports only 2D images!\n\n";
      return;
    }

    // Convert rgb to gray if three Channels are available
    numEntries = 1;
    for( int d= arrayDim.vec.size()-1 ; d>= 0 ; d-- )
      numEntries*= arrayDim.vec[d];

    return;
  }

  /** wellners gaussian adaptive 2D thresholding function, default mode is relative(mode = 0), mode = 1 is fixed.  */
  void CALTagPattern::adaptiveThreshold( Array<T> *inArray, double t, int fsize, int mode ) {
    if( inArray == NULL )
      return;

    CoordinateVector inDim = inArray->getDimension();

    // Copy and filter with gaussian and upsample
    Array<T> inArrayCopy( inDim );
    inArrayCopy.setNativeType( inArray->getNativeType() );
    inArrayCopy.addChannel();

    T* inDataScanline = new T[ inDim.vec[0] ];
    for(  unsigned long i = 0 ; i < inDim.vec[1]  ; i++ ) {
      inArray->getScanline(i, inDataScanline );
      inArrayCopy.updateScanline(i, inDataScanline );
    }
    delete[] inDataScanline;

    //Apply gaussian filtering
    if( fsize == -1 ) {
      long maxDim = 0;
      for( int d = 0; d < inDim.vec.size(); d++ )
        maxDim = (inDim.vec[d] > maxDim)? inDim.vec[d] : maxDim;

      fsize = floor( ((double)maxDim) / 20.0 );
    }

    FilterParams filterParams;
    filterParams.radius= fsize * 3.0;
    filterParams.sigma= fsize;
    filterParams.edgeStopSigma= 0.1;
    filterParams.normalize= false;
    FilterFactory<T> filterFactory( &filterParams );
    FilterType filterType= FastGaussFiltering;
    Filter<T> *filter= filterFactory.create( filterType );

    // Setup filter behaviour for 2D
    BoundaryMethod boundary= Clamp;
    ChannelList channels;
    channels.vec.push_back(0);

    AxisList axes;
    axes.vec.push_back(0);
    axes.vec.push_back(1);

    filter->apply( inArrayCopy, boundary, channels, axes );
    delete filter;

    // Finally apply the threshold
    T* dataCh = &((*((*inArray)[0]))[0]);
    T* filCh = &((*(inArrayCopy[0]))[0]);
    unsigned long numEntries = 1;
    for( int d = inDim.vec.size()-1 ; d>= 0 ; d-- )
      numEntries*= inDim.vec.at(d);

    // Different actions for relative or fixed thresholding
    double t_dec = 0.0;
    double comp = 1.0-(t/100.0);
    if( mode == 1 ) {
      t_dec = t;
      comp = 1.0;
    }
    for(  unsigned long i = 0 ; i < numEntries ; i++ )
      dataCh[i] = (dataCh[i] > (filCh[i])*comp - t_dec) ? 1.0 : 0.0;

    return;
  }

  /** rotates the code pattern by 90 degrees counterclockwise */
  void CALTagPattern::rotCode90() {
    vector<bool> tempCode( resCode*resCode );
    for(unsigned int i=0;i < resCode; i++)
      for(unsigned int j=0;j < resCode; j++)
        tempCode.at(i*resCode + j) = code->at(j*resCode  + (resCode - i - 1)) ;

    //Copy back
    for(unsigned int i=0;i < resCode*resCode; i++)
      code->at(i) = tempCode.at(i);
  }

  bool CALTagPattern::validateSaddle( Array<T> *in, Vector saddle, int minR, int maxR, int window, double alpha, double kappa, double delta ) {
    //Early exit
    int bound = maxR+window;
    if(  maxR <= minR ||
      !(saddle.get(0) > bound && saddle.get(0) < arrayDim.vec[0] - bound && saddle.get(1) > bound && saddle.get(1) < arrayDim.vec[1] - bound ) )
      return true;

    //Source data location
    T* srcData = &((*((*in)[0]))[0]);


    vector<T> currCircle;
    vector<T> currCircleDeriv;
    int numRadii = maxR - minR + 1;

    //Current center of the circle
    long centerX = (long)saddle.get(0);
    long centerY = (long)saddle.get(1);

    //Try to search at least one valid circle in a window of size "window"
    //This is the number of found valid circles (with four alternating regions) in the window
    int windowSaddles = 0;

    //C{x,y} is center of current circle in window
    for( long cy= centerY - window; cy <= centerY + window; cy++ ) {
      for( long cx= centerX - window; cx <= centerX + window  ; cx++ ) {

        //Now iterate of the radii
        double validCheckerCircles = 0.0;

        for( long r = minR; r<= maxR; r++) 
	{
		//Current circle
		currCircle.clear();
		//First derivative of circle
		currCircleDeriv.clear();

		//Rasterize a circle with modified bresenham circle
		int f = 1 - r;
		int ddF_x = 0;
		int ddF_y = -2 * r;
		int x = 0;
		int y = r;


		//Push first points of 0, 2, 4, 6 octand to queue
		currCircle.push_back( srcData[(cy + r) * arrayDim.vec[0] + cx] );
		currCircle.push_back( srcData[ cy * arrayDim.vec[0] + cx + r] );
		currCircle.push_back( srcData[(cy - r) * arrayDim.vec[0] + cx] );
		currCircle.push_back( srcData[cy * arrayDim.vec[0] + cx - r ] );

		int i = 1;
		while(x < y) 
		{
			if(f >= 0) 
			{
			      y--;
			      ddF_y += 2;
			      f += ddF_y;
			}
			x++;
			ddF_x += 2;
			f += ddF_x + 1;

			//Push pixel values for 8 octands of circle in sorted order to queue
			currCircle.insert( currCircle.begin() + i, srcData[(cy + y) * arrayDim.vec[0] + cx + x ] );
			currCircle.insert( currCircle.begin() + i + 1, srcData[(cy + x) * arrayDim.vec[0] + cx + y ] );
			currCircle.insert( currCircle.begin() + 3*i + 1, srcData[(cy - x) * arrayDim.vec[0] + cx + y ] );
			currCircle.insert( currCircle.begin() + 3*i + 2, srcData[(cy - y) * arrayDim.vec[0] + cx + x ] );
			currCircle.insert( currCircle.begin() + 5*i + 2, srcData[(cy - y ) * arrayDim.vec[0] + cx - x ] );
			currCircle.insert( currCircle.begin() + 5*i + 3, srcData[(cy - x)* arrayDim.vec[0] + cx - y  ] );
			currCircle.insert( currCircle.begin() + 7*i + 3, srcData[(cy + x) * arrayDim.vec[0] + cx - y  ] );
			currCircle.insert( currCircle.begin() + 7*i + 4, srcData[(cy + y) * arrayDim.vec[0] + cx - x ] );

			i++;
		}

		//Circle is filled by now. Do adaptive Thresholding
		double circleThresh = 0.0;
		for(int i = 0; i < currCircle.size(); i++)
			circleThresh += currCircle.at(i);

		circleThresh /= (double)(currCircle.size());

		//Now do the actual thresholding
		for(int i = 0; i < currCircle.size(); i++)
			currCircle.at(i) = (currCircle.at(i) > circleThresh) ? 1.0 : 0.0;

		//Now count the differing opposite pixels
		int halfCircle = currCircle.size()/2;
		int differingOppositePixels = 0;
		for(int i = 0; i < halfCircle; i++)
			differingOppositePixels += (currCircle.at(i) != currCircle.at( halfCircle + i ) ) ? 1 : 0;

		//Debug
		//cerr << "size: " << halfCircle << " diff: " << differingOppositePixels << endl;

		//No valid circle if more than delta percent of opposite circle pixels are different
		if( ((double)differingOppositePixels) > delta * ((double)halfCircle) )
			continue;

		//Calculate the absolute derivatives
		currCircleDeriv.push_back( fabs( currCircle.at(0) - currCircle.back() ) );
		int lastRegionBegin = -1;
		if( currCircle.at(0) - currCircle.back() > 0 )
			lastRegionBegin = 0;

		for(int i = 1; i < currCircle.size(); i++) 
		{
			if( currCircle.at(i) -  currCircle.at(i-1) > 0 )
			  lastRegionBegin = i;

			currCircleDeriv.push_back( fabs( currCircle.at(i) -  currCircle.at(i-1) ) );

			if( currCircle.at(i) -  currCircle.at(i-1) < 0 && lastRegionBegin != -1
			&& (double)(i - lastRegionBegin) < kappa* (double)(currCircle.size()) ) 
			{
			    currCircleDeriv.at(i ) = 0.0;
			    currCircleDeriv.at(lastRegionBegin) = 0.0;
			    lastRegionBegin = -1;
			}

		}

		 //Calculate the sum of derivs
		double sumDeriv = 0.0;
		for(int i = 1; i < currCircleDeriv.size(); i++)
			sumDeriv += currCircleDeriv.at(i);

		if( sumDeriv == 4.0 && halfCircle)
			validCheckerCircles += 1.0;

        }

        if( validCheckerCircles >= alpha * (double)numRadii )
		windowSaddles++;

      }
    }

              //At least one saddle in the window
    return windowSaddles > 0 ;
  }

  /** computes the 2d homography h from the points in to the points out. H needs to be a 3x3 matrix */
  void CALTagPattern::homography2D( const SampleVector* in, const SampleVector *out, Matrix* h  ) {
    //Early exit if wrong size
    if( in == NULL || out == NULL || (in->data).size() != (out->data).size()
      || h == NULL || h->getNumRows() != 3 || h->getNumColumns() != 3 )
      return;

    //Copy both samplevectors
    SampleVector inSamples;
    SampleVector outSamples;
    Vector tempVec(3, 0.0), zeroVec(3, 0.0);
    for( unsigned long n = 0; n < (in->data).size() ; n++ ) {
      if( (in->data).at(n).getSize() != 3  ||  (out->data).at(n).getSize() != 3 ) {
        cerr << "Wrong size for homography 2D!" << endl;
        return;
      }
      tempVec.copy( zeroVec );
      tempVec.assign( (in->data).at(n) ); inSamples.data.push_back( tempVec );
      tempVec.copy( zeroVec );
      tempVec.assign( (out->data).at(n) ); outSamples.data.push_back( tempVec );
    }

    //Normalize input and output points
    Matrix t1(3,3);
    Matrix t2(3,3);
    normalizePoints(&inSamples , &t1);
    normalizePoints(&outSamples , &t2);

    //Initialize and fill point-correspondence Matrix
    Matrix A( 3* inSamples.data.size() , 9 );

    Vector X(3, 0.0), Y(3, 0.0);
    double xs, ys, ws;
    for( unsigned long n = 0; n < inSamples.data.size() ; n++ ) {

      X.assign( inSamples.data.at(n) );
      Y.assign( outSamples.data.at(n) );
      xs = Y[0]; ys = Y[1]; ws = Y[2];

      //Set whole 3x3 submatrix for points n
      A.set( 3*n, 0, 0.0 );
      A.set( 3*n, 1, 0.0 );
      A.set( 3*n, 2, 0.0 );
      A.set( 3*n, 3, -ws * X[0] );
      A.set( 3*n, 4, -ws * X[1] );
      A.set( 3*n, 5, -ws * X[2] );
      A.set( 3*n, 6, ys * X[0] );
      A.set( 3*n, 7, ys * X[1] );
      A.set( 3*n, 8, ys * X[2] );

      A.set( 3*n + 1, 0, ws * X[0] );
      A.set( 3*n + 1, 1, ws * X[1] );
      A.set( 3*n + 1, 2, ws * X[2] );
      A.set( 3*n + 1, 3, 0.0 );
      A.set( 3*n + 1, 4, 0.0 );
      A.set( 3*n + 1, 5, 0.0 );
      A.set( 3*n + 1, 6, -xs * X[0] );
      A.set( 3*n + 1, 7, -xs * X[1] );
      A.set( 3*n + 1, 8, -xs * X[2] );

      A.set( 3*n + 2, 0, -ys * X[0] );
      A.set( 3*n + 2, 1, -ys * X[1] );
      A.set( 3*n + 2, 2, -ys * X[2] );
      A.set( 3*n + 2, 3, xs * X[0] );
      A.set( 3*n + 2, 4, xs * X[1] );
      A.set( 3*n + 2, 5, xs * X[2] );
      A.set( 3*n + 2, 6, 0.0 );
      A.set( 3*n + 2, 7, 0.0 );
      A.set( 3*n + 2, 8, 0.0 );
    }

    // Find eigenvalues for matrix (A^T * A)
    Matrix AtA( 9, 9 );
    multMatrixMatrix( A.getTranspose(), A, AtA );

    // Solve Eigenvalue problem
    JacobiRotation solver;
    Vector eigValues( 9 );
    Matrix eigVectors( 9, 9 );
    solver.solve( AtA, eigValues, eigVectors );

    // Find the eigenVector with the least eigenValue
    unsigned long leastEigenVector = 0;
    for( unsigned long i = 1; i < 9; i++ ) {
      if( eigValues[i] < eigValues[leastEigenVector] )
        leastEigenVector = i;
    }

    // Assemble h matrix
    for( unsigned int i = 0; i < 9; i++) {
      h->set( (i- i%3)/3, i % 3,  eigVectors[leastEigenVector][i]);
    }

    //denormalize
    Matrix t2_inv( 3, 3 );
    t2_inv = inverse( t2 );
    Matrix h_temp( 3, 3 );
    multMatrixMatrix( t2_inv, (*h), h_temp );
    multMatrixMatrix( h_temp, t1, (*h) );

    return;
  }

  /** normalizes a set of points ˆXi = T*Xi such that the centroid of the new points
      has the coordinate (0, 0) and their average distance from the origin is sqrt(2). */
  void CALTagPattern::normalizePoints( SampleVector* points, Matrix* t ) {
    //Early exit if wrong size
    if( points == NULL || (points->data).size() == 0 ||
      t == NULL || t->getNumRows() != 3 || t->getNumColumns() != 3 )
      return;

    //Initialize transformation matrix
    t->zero();
    t->set(2,2, 1.0 );

    //Compute mean of all points
    double xmean = 0.0, ymean = 0.0;
    for( unsigned long i = 0; i < (points->data).size() ; i++ ) {
      if( (points->data).at(i).getSize() != 3  ) {
        cerr << "Wrong size for homography 2D!" << endl;
        return;
      }
      //First normalize
      ((points->data).at(i))[0] /= ((points->data).at(i))[2];
      ((points->data).at(i))[1] /= ((points->data).at(i))[2];
      ((points->data).at(i))[2] = 1.0;

      xmean += ((points->data).at(i)).get(0);
      ymean += ((points->data).at(i)).get(1);
    }
    xmean /= (double)((points->data).size());
    ymean /= (double)((points->data).size());

    //Compute s
    double s = 0.0;
    for( unsigned long i = 0; i < (points->data).size() ; i++ ) {
      s += sqrt( pow( ((points->data).at(i)).get(0) - xmean , 2.0) + pow( ((points->data).at(i)).get(1) - ymean, 2.0) );
    }
    s /= (double)((points->data).size());
    s= (sqrt(2.0))/s;

    //Save transformation matrix
    t->set( 0, 0, s );
    t->set( 1, 1, s );
    t->set( 0, 2, -s*xmean);
    t->set( 1, 2, -s*ymean);

    //Apply transformation to all points
    double x, y;
    Vector resVect( 3, 0.0);
    for( unsigned long i = 0; i < (points->data).size() ; i++ ) {
      t->rightMultiply( (points->data).at(i), resVect );
      ((points->data).at(i)).assign( resVect );
    }
    return;

  }

  /** applys a 2D homography h to a bunch of points */
  void CALTagPattern::homotrans( const SampleVector *points, Matrix* h, SampleVector *result  ) {
    //Early exit if wrong size
    if( points == NULL || (points->data).size() == 0 ||
      h == NULL || h->getNumRows() != 3 || h->getNumColumns() != 3 )
      return;

    Vector resVect( 3, 0.0);
    Vector newVect( 3, 0.0);
    (result->data).clear();
    for( unsigned long i = 0; i < (points->data).size(); i++ ) {
      if( ((points->data).at(i)).getSize() != 3 ) {
        cerr << "Transformation matrix and point dimensions do not match! " << endl;
        return;
      }
      resVect.copy( newVect );
      h->rightMultiply( (points->data).at(i), resVect );

      //Normalize
      for( int j = 0; j < 2; j++ )
        resVect[j] /= resVect[2];

      resVect[2] = 1.0;
      (result->data).push_back( resVect );
    }
    return;
  }

  // Ransac functions:
  /** Function to determine if a set of 4 2D points give rise
      to a degeneracy in the calculation of a homography as needed by RANSAC.
      This involves testing whether any 3 of the 4 points the set is collinear. */
  bool CALTagPattern::isDegenerate( SampleVector* points, double threshold ) {
    //Early exit if wrong size
    if( points == NULL || (points->data).size() != 4 )
      return false;

    int i,j,k;
    Vector p1(2, 0.0), p2(2, 0.0), p3(2, 0.0);
    bool linear = false;
    // check for each 3 points combination
    for(i=0; i < (points->data).size() - 2; i++) {
      p1.set(0, (points->data).at(i).get(0) );
      p1.set(1, (points->data).at(i).get(1) );
      for(j=i+1; j < (points->data).size() - 1; j++) {
        p2.set(0, (points->data).at(j).get(0) );
        p2.set(1, (points->data).at(j).get(1) );
        for(k=j+1; k < (points->data).size() ; k++) {
          p3.set(0, (points->data).at(k).get(0) );
          p3.set(1, (points->data).at(k).get(1) );
              //Crossproduct == 0.0 ?
          linear = ( fabs((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1])) < threshold );
          if(linear)
            break;
        }
        if(linear)
          break;
      }
      if(linear)
        break;
    }

    return linear;
  }

  /** Function to evaluate the symmetric transfer error of a homography with
    respect to a set of matched points as needed by RANSAC. */
  void CALTagPattern::homogdist2d(const SampleVector *pointsIn, const SampleVector *pointsOut, Matrix* H, double t, vector<unsigned long>* inliers ) {
    //Early exit if wrong size
    if( pointsIn == NULL || pointsOut == NULL || (pointsIn->data).size() != (pointsOut->data).size() || inliers == NULL
      || H == NULL || H->getNumRows() != 3 || H->getNumColumns() != 3 )
      return;

    //Precompute Inverse
    Matrix invH( 3, 3 );
    invH.assign(*H);
    invH = inverse( invH );

    //Initialize inliers
    inliers->clear();

    Vector p1(3, 0.0), p2(3, 0.0), Hp1(3, 0.0), invHp2(3, 0.0);
    for( unsigned long i = 0; i < (pointsIn->data).size(); i++ ) {
      if( ((pointsIn->data).at(i)).getSize() != 3 || ((pointsOut->data).at(i)).getSize() != 3 ) {
        cerr << "Transformation matrix and point dimensions do not match! " << endl;
        return;
      }
      p1.assign( (pointsIn->data).at(i) );
      p2.assign( (pointsOut->data).at(i) );
      H->rightMultiply( p1, Hp1 );
      invH.rightMultiply( p2, invHp2 );

      //Normalize
      for( int j = 0; j < 2; j++ ) {
        p1[j] /= p1[2];
        p2[j] /= p2[2];
        Hp1[j] /= Hp1[2];
        invHp2[j] /= invHp2[2];
      }

      p1[2] = 1.0;
      p2[2] = 1.0;
      Hp1[2] = 1.0;
      invHp2[2] = 1.0;

      p1 -= invHp2;
      p2 -= Hp1;

      if( fabs( dot( p1, p1 ) + dot( p2, p2 ) ) < t )
        inliers->push_back( i );
    }

    return;
  }

  /** fits 2D homography using the RANSAC algorithm. Returns true if successful*/
  bool CALTagPattern::ransacFitHomography(const SampleVector *pointsIn, const SampleVector *pointsOut, Matrix* H,
    vector<unsigned long>* inliers, double t, unsigned long s, unsigned long feedback,
  unsigned long maxDataTrials, unsigned long maxTrials ) {
    //Early exit if wrong size
    if( pointsIn == NULL || pointsOut == NULL || (pointsIn->data).size() != (pointsOut->data).size() || inliers == NULL
      || H == NULL || H->getNumRows() != 3 || H->getNumColumns() != 3 )
      return false;

    // Number of sampleIn<->sampleOut correspondences
    unsigned long npts = (pointsIn->data).size();

    // No inliers at beginning
    inliers->clear();

    //Check for minimum size
    if( npts < 4 ) {
      cerr << "Must have at least 4 points to fit homography!" << endl;
      return false;
    }

    //Copy both samplevectors
    SampleVector inSamples;
    SampleVector outSamples;
    Vector tempVec(3, 0.0), zeroVec(3, 0.0);
    for( unsigned long n = 0; n < npts ; n++ ) {
      if( (pointsIn->data).at(n).getSize() != 3  ||  (pointsOut->data).at(n).getSize() != 3 ) {
        cerr << "Wrong size for ransacFitHomography 2D!" << endl;
        return false;
      }
      tempVec.copy( zeroVec );
      tempVec.assign( (pointsIn->data).at(n) ); inSamples.data.push_back( tempVec );
      tempVec.copy( zeroVec );
      tempVec.assign( (pointsOut->data).at(n) ); outSamples.data.push_back( tempVec );
    }

    // Normalise each set of points so that the origin is at centroid and
    // mean distance from origin is sqrt(2).  normalise2dpts also ensures the
    // scale parameter is 1.
    Matrix t1(3,3);
    Matrix t2(3,3);
    normalizePoints(&inSamples , &t1);
    normalizePoints(&outSamples , &t2);

    //Initialize counting variables
    // Desired probability of choosing at least one sample free from outliers
    double p = 0.99;

    unsigned long trialcount = 0;
    unsigned long bestscore =  0;
    // Dummy initialisation for number of trials.
    long N = 1;

    //Initialize random number generator
    srand48(time(0));

    bool degenerate;
    unsigned long count;
    SampleVector randomSamplesIn,randomSamplesOut;
    vector<long> randNumbers;
    Matrix M(3,3,true), bestM(3,3,true);
    unsigned long ninliers;
    vector<unsigned long> bestinliers;
    bestinliers.clear();
    //Begin iteration
    while( N > trialcount ) {
      // Select at random s datapoints to form a trial model, M.
      // In selecting these points we have to check that they are not in
      // a degenerate configuration.
      degenerate = true;
      count = 1;
      while( degenerate ) {
              // Generate s random different(!) indicies in the range 0..npts-1
        randomSamplesIn.data.clear();
        randomSamplesOut.data.clear();
        randNumbers.clear();

        long currRand;
        for( unsigned long i = 0; i < s ; i++ ) {
          if( (npts - 1 - i) <= 0 )
            currRand = 0;
          else
            currRand = lrand48() % (npts - 1 - i);

              //Order sample into indiceslist
          unsigned long j = 0;
          for( j = 0; j < i; j++ ) {
            if( currRand >= randNumbers.at(j) )
              currRand++;
            else
              break;
          }
          randNumbers.insert( randNumbers.begin() + j, currRand );
              /*Debug output
          cerr << "CurrRand: " << currRand << " ";*/

              //Insert according Sample to randomSamples
          tempVec.copy( zeroVec );
          tempVec.assign( inSamples.data.at(currRand) );
          randomSamplesIn.data.push_back( tempVec );
          tempVec.copy( zeroVec );
          tempVec.assign( outSamples.data.at(currRand) );
          randomSamplesOut.data.push_back( tempVec );
        }

              // Test that these points are not a degenerate configuration.
        degenerate = isDegenerate( &randomSamplesIn ) || isDegenerate( &randomSamplesOut );
              /*Debug output
        cerr << " Degenerate: " << degenerate << endl;*/

        if( !degenerate ) {
              // Fit model M to this random selection of data points.
          homography2D(&randomSamplesIn,  &randomSamplesOut, &M);
        }

              // Safeguard against being stuck in this loop forever
        count = count + 1;
        if( count > maxDataTrials ) {
          cerr << "Unable to select a nondegenerate data set!" << endl;
          return false;
        }
      }

      // Once we are out here we should have some kind of model...
      // Evaluate distances between points and model returning the indices
      // of elements in x that are inliers.
      inliers->clear();
      homogdist2d(&inSamples, &outSamples, &M, t, inliers );

      // Find the number of inliers to this model.
      ninliers = inliers->size();

      // Largest set of inliers so far...
      if( ninliers > bestscore ) {
              // Record data for this model
        bestscore = ninliers;
        bestinliers.clear();
        for( unsigned long i = 0; i < inliers->size(); i++ )
          bestinliers.push_back( inliers->at(i) );

        bestM.copy( M );

              // Update estimate of N, the number of trials to ensure we pick,
              // with probability p, a data set with no outliers.
        double fracinliers =  (double)ninliers / (double)npts;
        double pNoOutliers = 1.0 -  pow(fracinliers, s);

        N = (long)( log(1.0 - p) / log( pNoOutliers ) );
      }

      trialcount = trialcount+1;
      /* Debug output
      cerr << "Trial " << trialcount << " out of " << N+1 << endl;*/

      // Safeguard against being stuck in this loop forever
      if( trialcount > maxTrials ) {
        cerr << "Ransac reached the maximum number of " << maxTrials << " trials." << endl;
        return false;
      }
    }

    // We got a solution
    if( bestinliers.size() > 0 ) {
      H->copy(bestM);
      inliers->clear();
      for( unsigned long i = 0; i < bestinliers.size(); i++ )
        inliers->push_back( bestinliers.at(i) );
    }
    else {
      H->identity();
      inliers->clear();
      cerr << "Ransac was unable to find a useful solution" << endl;
      return false;
    }

    // Now do a final least squares fit on the data points considered to be inliers.
    SampleVector inliersIn, inliersOut;

    for( unsigned long i = 0; i < inliers->size(); i++ ) {
      tempVec.copy( zeroVec );
      tempVec.assign( inSamples.data.at(inliers->at(i)) );
      inliersIn.data.push_back( tempVec );
      tempVec.copy( zeroVec );
      tempVec.assign( outSamples.data.at(inliers->at(i)) );
      inliersOut.data.push_back( tempVec );
    }
    homography2D(&inliersIn,  &inliersOut, H);

    // Denormalise
    Matrix t2_inv( 3, 3 );
    t2_inv = inverse( t2 );
    Matrix h_temp( 3, 3 );
    multMatrixMatrix( t2_inv, (*H), h_temp );
    multMatrixMatrix( h_temp, t1, (*H) );

    return true;
  }

/** corrects distortion for a image coordinate vector */
  void CALTagPattern::correctDistortion( Vector *point, Vector *imgCenter, double halfDiagonal, double r)
  {
	*point -= *imgCenter;
	
	//Normalized distance from center to point that corners have the distance 1.0
	//where halfDiagonal is the length of the half diagonal of the image
	double normShiftVecLength = point->norm(2) / halfDiagonal;

	//Compute correction factor
	double correctionFactor = pow((1.0 - normShiftVecLength), r );
	
	//Scale shift vector and compute corrected point coords
	*point *= correctionFactor;
	*point += *imgCenter;

	return;
  }

  /** corrects distortion for a couple of points */
  void CALTagPattern::correctDistortion( SampleVector *points, Vector *imgCenter, double halfDiagonal, double r)
  {
	for(unsigned int i = 0; i < (points->data).size(); i++)
	      correctDistortion( &((points->data).at(i)), imgCenter, halfDiagonal, r);
	return;
  }

  /** find an r which makes the points collinear */
  double CALTagPattern::findUndistortR(const SampleVector *points, Vector *imgCenter, double halfDiagonal, double precision )
  {
	//Precision of method
	long maxIterations = 1000;

	//Increment definitions
	double increment = 1.0;
	double oldDist = computeDistanceToBestFittedLine( points );
	double oldR = 0.0;
	
	if( oldDist < precision )
	    return oldR;
	
	//Search for an r in both directions
	double newDist = oldDist;
	double newR = oldR;
	SampleVector currLine; Vector tempVec(2, 0.0), zeroVec(2, 0.0);
	for( int iter = 0; iter < maxIterations; iter++ )
	{
	    //Check both directions
	    for( double dir = 1.0; dir > -2.0; dir -= 2.0 )
	    {
		do
		{
		    oldDist = newDist;
		    oldR = newR;
		    newR += increment * dir;

		    //Copy to currLine
		    currLine.data.clear();
		    for( unsigned long n = 0; n < (points->data).size() ; n++ )
		    {
			tempVec.copy( zeroVec ); 
			tempVec.assign( (points->data).at(n) ); currLine.data.push_back( tempVec );
		    }
		    
		    //Compute new collinearity rating
		    correctDistortion( &currLine, imgCenter, halfDiagonal, newR);
		    newDist = computeDistanceToBestFittedLine( &currLine );

		}while( newDist < oldDist );
	    }
	    
	    if( oldDist < precision )
		return oldR;
	    else
		increment /= 2.0;
	}

	return oldR;
  }

  /** find an r which maps the pointsIn to pointsOut by minimizing the distance between mapped points and pointsOut */
  double CALTagPattern::findInversionR( const SampleVector *pointsIn, const SampleVector *pointsOut, Vector *imgCenter, double halfDiagonal, double precision )
  {
	//Precision of method
	long maxIterations = 1000;

	//Increment definitions
	double increment = 1.0;
	double oldDist = computePointSetDifference( pointsIn, pointsOut );
	double oldR = 0.0;
	
	if( oldDist < precision )
	    return oldR;
	
	//Search for an r in both directions
	double newDist = oldDist;
	double newR = oldR;
	SampleVector guessedPointsOut; Vector tempVec(2, 0.0), zeroVec(2, 0.0);
	for( int iter = 0; iter < maxIterations; iter++ )
	{
	    //Check both directions
	    for( double dir = 1.0; dir > -2.0; dir -= 2.0 )
	    {
		do
		{
		    oldDist = newDist;
		    oldR = newR;
		    newR += increment * dir;

		    //Copy pointsIn to guessedPointsOut
		    guessedPointsOut.data.clear();
		    for( unsigned long n = 0; n < (pointsIn->data).size() ; n++ )
		    {
			tempVec.copy( zeroVec ); 
			tempVec.assign( (pointsIn->data).at(n) ); guessedPointsOut.data.push_back( tempVec );
		    }
		    
		    //Fit model with newR and compute new point set diff
		    correctDistortion( &guessedPointsOut, imgCenter, halfDiagonal, newR);
		    newDist = computePointSetDifference( &guessedPointsOut, pointsOut );

		}while( newDist < oldDist );
	    }
	    
	    if( oldDist < precision )
		return oldR;
	    else
		increment /= 2.0;
	}

	return oldR;
  }

  /** compute summed distance between points */
  double CALTagPattern::computePointSetDifference( const SampleVector *pointsIn, const SampleVector *pointsOut  )
  {
	if( (pointsIn->data).size() != (pointsOut->data).size() )
	{
	    cerr << "Error computing point set difference, assuming same number of compared SampleVectors" << endl;
	    return 0.0;
	}

	double totalDist = 0.0;
	for( unsigned int i = 0; i < (pointsIn->data).size(); i++ )
	    totalDist += dist(  (pointsIn->data).at(i), (pointsOut->data).at(i) );

	return totalDist;
  }

  /** rate collinearity */
  double CALTagPattern::computeDistanceToBestFittedLine( const SampleVector *points )
  {
	Line bestFitLine(2);
	fitToDataPoints( bestFitLine, *points );

	double totalDist = 0.0;
	for( unsigned int i = 0; i < (points->data).size(); i++ )
	    totalDist += dist( bestFitLine, (points->data).at(i) );

	return totalDist;
  }


  /*[brad] based on Felix's "track" function, but using the CalibrationPattern::detect().
 Computes the tracking for both checkerboard layouts. optionally generates debug output in specified fig-file
      debug = 0 ( no debugging information output )
      debug = 1 ( minimal debugging information ( Numbers of found markers, and if successful tracked )
      debug = 2 ( all available debugging information ) 
      debug = 3 ( output fig-file with sampled code positions, found saddles and coordinate system axis ) */
  bool CALTagPattern::detect() {
    //[brad] this stuff came from the CALTag constructor, but since now the CALTagPattern object does
    //stuff other than just detecting, it doesn't make sense to put it in the constructor, so we'll
    //just do it here
    debug = 0;
debugOutputFilename = "debug";
debugScaleFactor = 1.0;
//cerr << "debug:1" << endl;
    numMarkers = rows * columns;
    idBits = bits - cBits;
    successfulTracked = false;
    unitCalibPoints = new SampleVector();
    imageCalibPoints = new SampleVector();

    // width of border in units around the code dots in each square
    borderWidth = (double)(resSquare - resCode) / 2.0;

    // Compute unit sized calibration grid
    ((this->unitCalibPoints)->data).clear();
    Vector tmpPnt(3, 0.0);
    Vector oneVec(3, 1.0);
    for( unsigned int row = 0; row < rows+1 ; row++ ) {
      for( unsigned int col = 0; col < columns+1 ; col++ ) {
        tmpPnt.copy( oneVec );
        tmpPnt.set( 0, (double)row + 1.0 );
        tmpPnt.set( 1, (double)col + 1.0);
        ((this->unitCalibPoints)->data).push_back( tmpPnt );
      }
    }

    //Prepare matrices for homography
    // create a unit square:
    //[brad] original caltag coordinate system
    //unitSquare = new SampleVector();
    //tmpPnt.copy( oneVec );
    //tmpPnt[0] = 0.0; tmpPnt[1] = 0.0; (unitSquare->data).push_back(tmpPnt);
    //tmpPnt.copy( oneVec );
    //tmpPnt[0] = 0.0; tmpPnt[1] = 1.0; (unitSquare->data).push_back(tmpPnt);
    //tmpPnt.copy( oneVec );
    //tmpPnt[0] = 1.0; tmpPnt[1] = 1.0; (unitSquare->data).push_back(tmpPnt);
    //tmpPnt.copy( oneVec );
    //tmpPnt[0] = 1.0; tmpPnt[1] = 0.0; (unitSquare->data).push_back(tmpPnt);
    //[brad] changing now to the postscript coordinate system
    unitSquare = new SampleVector();
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 0.0; tmpPnt[1] = 1.0; (unitSquare->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 0.0; tmpPnt[1] = 0.0; (unitSquare->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 1.0; tmpPnt[1] = 0.0; (unitSquare->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 1.0; tmpPnt[1] = 1.0; (unitSquare->data).push_back(tmpPnt);

    //create bowtie intersection positions:
    bowtieIntersections = new SampleVector();
    tmpPnt.copy( oneVec );
    tmpPnt[0] = -0.5; tmpPnt[1] = 0.5; (bowtieIntersections->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 0.5; tmpPnt[1] = -0.5; (bowtieIntersections->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 1.5; tmpPnt[1] = 0.5; (bowtieIntersections->data).push_back(tmpPnt);
    tmpPnt.copy( oneVec );
    tmpPnt[0] = 0.5; tmpPnt[1] = 1.5; (bowtieIntersections->data).push_back(tmpPnt);

    // specify the points at which we want to sample (in unit square space)
    codeBlobSamplePos = new SampleVector();
    Vector blobPos(3, 0.0);

    if( layout == 1 ) //Original layout
    {
        for( int y = 2*borderWidth+1; y <= resSquare*2 - (2*borderWidth+1); y+=2 ) {
          for( int x = 2*borderWidth+1; x <= resSquare*2 - (2*borderWidth+1); x+=2 ) {
            blobPos.copy( oneVec );
            blobPos[0] = (double)x/ (double)(resSquare*2);
            blobPos[1] = (double)y/ (double)(resSquare*2);
            (codeBlobSamplePos->data).push_back(blobPos);
          }
        }
    } else if ( layout == 2 ) //Rotated layout
    {
        for( int y = 0; y < resCode; y++ )
        {
            for( int x = 0; x < resCode; x++ )
            {
                  blobPos.copy( oneVec );
		  blobPos[0] = ( borderWidth * 2.0 + 1.0 + (double)y + (double)x )/ (double)(resSquare*2);
		  blobPos[1] = ((double)(resSquare) - (double)x + (double)y )/ (double)(resSquare*2);
                  (codeBlobSamplePos->data).push_back(blobPos);
            }
        }
    }

    // current code pattern
    code = new vector<bool>(resCode * resCode);


    //[brad] end of stuff moved from the constructor

//cerr << "debug:2" << endl;
    // Early exit if wrong dimension
    if( arrayDim.vec.size() != 2 ) {
      cerr << ": Supports only 2D images!\n\n";
      //if( successfulTracked != NULL )
      successfulTracked = false;
      return false;
    }

    // RGB to grey conversion and copy to new array, since old greyscale array is needed later
    // for saddlepoint finding.
    Array<T> connCompArray( arrayDim );
    connCompArray.setNativeType((this->array)->getNativeType() );
    connCompArray.addChannel();
    T* connCompData = &((*(connCompArray[0]))[0]);

    T* dataCh0 = &((*((*(this->array))[0]))[0]);
    if( (this->array)->getNumChannels() >= 3 ) {
      unsigned int resultCh = (this->array)->addChannel();

      T* dataCh1 = &((*((*(this->array))[1]))[0]);
      T* dataCh2 = &((*((*(this->array))[2]))[0]);
      T* dataResultCh = &((*((*(this->array))[resultCh]))[0]);

      for(  unsigned long i = 0 ; i < numEntries ; i++ ) {
        dataResultCh[i] = 0.2989 * dataCh0[i] + 0.5870 * dataCh1[i] + 0.1140 * dataCh2[i];
        connCompData[i] = dataResultCh[i];
      }

      //Store result channel as first channel
      (this->array)->swapChannels( 0, resultCh );
    }else
    {
      for(  unsigned long i = 0 ; i < numEntries ; i++ )
        connCompData[i] = dataCh0[i];
    }

    //First channel is used, all other channels are deleted
    for( unsigned int i = (this->array)->getNumChannels() - 1; i > 0 ; i--)
      (this->array)->deleteChannel(i);

    // Output array as png for debug
    string greyscaleImageFilename;
    if( debug == 2 || debug == 3 ) {
      // setup ImageMagick image
      ostringstream geometry;
      geometry << arrayDim.vec[0] << 'x' << arrayDim.vec[1];
      Image outImage( geometry.str().c_str(), "black" );
#if defined (_WIN32) || defined (_WIN64)
	  outImage.type(GrayscaleType);
#else
      outImage.type( GrayscaleMatteType );
#endif
      outImage.modifyImage();
      Pixels view( outImage );

      T *scanlineIn = new T[ arrayDim.vec[0] ];
      unsigned short *scanlineOut= new unsigned short[arrayDim.vec[0]];

      // copy image data
      for( unsigned long i= 0 ; i< arrayDim.vec[1] ; i++ ) {
        // read scanline and convert to unsigned short
        (this->array)->getScanline(i, scanlineIn);
        typeConvert( scanlineIn, Double, scanlineOut, UShort,
          arrayDim.vec[0]);
#if defined (_WIN32) || defined (_WIN64)
		Quantum* pixels = view.get(0, i, arrayDim.vec[0], 1);
		for (unsigned long j = 0; j < arrayDim.vec[0]; j++, pixels++){
			//*pixels = (scanlineOut[j] << 24) | (scanlineOut[j] << 16) | (scanlineOut[j] << 8) + (0);
			*pixels++ = scanlineOut[j];
			*pixels++ = scanlineOut[j];
			*pixels++ = scanlineOut[j];
			*pixels++ = 0;
		}
#else
        // get scanline from the ImageMagick image for write access
        PixelPacket *pixels= view.get( 0, i, arrayDim.vec[0], 1 );

        // now copy all the channels for all pixels in the scanline
        for( unsigned long j= 0 ; j< arrayDim.vec[0] ; j++, pixels++ ) {

          pixels->red= scanlineOut[j];
          pixels->green= scanlineOut[j];
          pixels->blue= scanlineOut[j];
          pixels->opacity= 0;
        }
#endif

        // update the actual image with the scanline we just copied
        view.sync();
      }
      //[brad] changed debugOutputFilename from char* to std::string to get rid of
      //g++ deprecated conversion warning
      greyscaleImageFilename.assign( debugOutputFilename.c_str() );
      greyscaleImageFilename.append( ".png" );
      outImage.write( greyscaleImageFilename );
    }

    FilterParams filterParams;
    filterParams.radius= 3.0; //Default parameters
    filterParams.sigma= 1.0;
    filterParams.edgeStopSigma= 0.1;
    filterParams.normalize= false;
    FilterType filterType;
    FilterFactory<T> filterFactory( &filterParams );
    Filter<T> *filter;
    BoundaryMethod boundary;
    ChannelList channels;
    AxisList axes;
    axes.vec.push_back(0);
    axes.vec.push_back(1);
    if( layout == 1 ) //Original layout
    {
      // Apply a sobel operator as edge filter
      // gives similar good results as log in matlab version
      // the result is stored in a seperate Array, as the original scaled greyscale array
      // is used for later saddlepoint finding
      filterParams.radius= 5.0;
      filterType= SobelFiltering;
      filter= filterFactory.create( filterType );

      // Setup filter behaviour for 2D
      boundary= Clamp;
      channels.vec.clear();;
      channels.vec.push_back(0);

      connCompArray.setNativeType( Double );
      filter->apply( connCompArray, boundary, channels, axes );
      delete filter;

      // compute radius for gaussian mean finder as length of minimal marker size
      double minMarkerArea= round( ((double)(resSquare*resSquare))*2.0*2.0 );

      // Converting to a binary edge / non-edge mask and inversion
      // Use adaptive threshold to find edges
      // Try to estimate mean of whole checkerboard region
      // Use very hard adaptive threshold ( 95 percent ) to determine clear black and white difference
      adaptiveThreshold( &connCompArray, 5, sqrt(minMarkerArea ) * ((double)(min(rows, columns))) / 2.0 );
      
      // Thinning filter
      filterParams.radius = 0.0;
      filterType = ThinningFiltering;
      filter= filterFactory.create( filterType );
      boundary= Clamp;
      channels.vec.clear();
      channels.vec.push_back(0);

      axes.vec.clear();
      axes.vec.push_back(0);
      axes.vec.push_back(1);
      filter->apply( connCompArray, boundary, channels, axes );
      delete filter;
    }else if( layout == 2 ) //Changed layout
    {
      // Adaptive threshold to account for shadows and highlights
      adaptiveThreshold( &connCompArray, 15 );
    }

    // Prepare for connected-components filtering ( edges-binary image is the result )
    connCompData = &((*(connCompArray[0]))[0]);

    for(  unsigned long i = 0 ; i < numEntries ; i++ )
      connCompData[i] = (connCompData[i] == 1.0) ? 0.0 : 0.5;

    // #### Debug output
    if( debug == 2 ) {
      string foundEdgesFilename;
      foundEdgesFilename.assign( debugOutputFilename.c_str() );
      foundEdgesFilename.append( "_foundEdges.mda" );
      //Output edge image
      connCompArray.write( (char*)foundEdgesFilename.c_str() ,Double);
    }

    // Label grid ( for connected component filter bgr = 0.0, 0.0 < foreground < 1.0 )
    filterParams.radius= 1.0;
    filterType = ConnectedComponentFiltering;
    filter= filterFactory.create( filterType );
    boundary= Clamp;
    channels.vec.clear();
    channels.vec.push_back(0);

    axes.vec.clear();
    axes.vec.push_back(0);
    axes.vec.push_back(1);
    filter->apply( connCompArray, boundary, channels, axes );
    delete filter;

//cerr << "debug:5" << endl;
    // Compute Number of Labels
    ConnectedComponentProperties<T> ccp( &connCompArray );
    ccp.computeNumLabels();

    // Debug output
    if( debug == 2 )
      cerr << "Numlabels: " << ccp.getNumLabels() << std::endl;

    ccp.computeIndependentProps(true, true);

    // Sort out regions with area < minArea or area > 1/8 of image size
    std::vector<unsigned int>* regions = new std::vector<unsigned int>();

    //Compute min and max Area
    double maxArea = ( (double)(arrayDim.vec[0] * arrayDim.vec[1]) ) / 8.0;
    // minimum (unfilled) square area assumes at least 2x2 pixels per code dot
    double minArea = round( ((double)(resSquare*resSquare))*4.0 - ((double)(resCode*resCode))*4.0 );

    // Sort out iteration:
    double area;
    for( int i = 2; i <= ccp.getNumLabels() ; i++) {
      area = (double)(ccp.getArea(i));
      if( area >= minArea && area < maxArea )
        regions->push_back( i );
    }

    // Extract segments
    ccp.computeSegments( 5, regions );
    if( debug == 2 )
      cerr << "Found " << ccp.getNumSegments() << " segments." << std::endl;

    // Compute eulernumbers and sort out segments with wrong eulernumber
    ccp.computeEulernumbers();
    regions->clear();
    double currEulerNum;
    for( int i = 1; i <= ccp.getNumLabels() ; i++) {
      if( !ccp.isValid(i) )
        continue;

      currEulerNum = round( ccp.getEulerNumber( i ) );
      if(  currEulerNum <= 0.0 &&  currEulerNum > (double)dotThresh)
        regions->push_back(i);
    }

//cerr << "debug:6" << endl;

    // #### Debug output
    if( debug == 2 ) {
      //Deassemble all extracted regions and output
      Array<T>* result = ccp.reassemble( regions );

      string foundMarkersFilename;
      foundMarkersFilename.assign( debugOutputFilename.c_str() );
      foundMarkersFilename.append( "_foundMarkers.mda" );
      //Output edge image
      result->write( (char*)foundMarkersFilename.c_str() ,Double);
      delete result;
    }

    // early exit if we detect too few potential markers
    if( (double)(regions->size()) < ((double)numMarkers)/8.0 ) {
      if( debug == 2 )
	cerr << "Not enough markers were detected !" << endl;
      delete regions;
      //if( successfulTracked != NULL )
      successfulTracked = false;
      return false;
    }
    
    if( debug == 2 )
      cerr << "Detected " << regions->size() << " markers." << std::endl;

    // #### Debug output:
    XFigWriter writer;
    Vector drawOffset( 2, 0.0);
    double circleSize = 9.0;
    Vector debugPoint(2, 0.0);
    if( debug == 2 || debug == 3 ) {
      string outputFigFilename;
      outputFigFilename.assign( debugOutputFilename.c_str() );
      outputFigFilename.append( ".fig" );

      // Setup debug output writer for following fig output
      writer.open( (char*)outputFigFilename.c_str() );
      writer.getStyle().thickness= 6.0;

      // Writer offset and scaling
      drawOffset[0] = 0.5;
      drawOffset[1] = 0.5;

      //Output corner points in green
      writer.getStyle().penColor = 2;
      writer.getStyle().depth= 69;

    }

    /* #### Outline output
    if( debug == 2 )
    {
        //Deassemble all extracted regions and output
        ccp.computeOutline();
        Array<T>* result = ccp.reassemble( regions );
        result->write( (char*)"outlines.mda" ,Double);
        delete result;

        if( successfulTracked != NULL )
            successfulTracked = false;
        return;
    }*/

//cerr << "debug:7" << endl;
    // Extract outline samples of all segments
    vector<SampleVector*>* outlineSamples = ccp.extractOutlineSamples( sampleDist, regions );

    // Iterate over all valid regions, fit Quadliteral and reassemble corners
    SampleVector* cornerSamples = new SampleVector();
    vector<long> cornerSampleRegion;
    SampleVector* currOutline;
    lloydQuadFinder quadFinder;
    SampleVector markerCorners;

    //Ordering of outlineSamples is the same as regions
    Vector currSample( 2, 0.0 ), zeroVec( 2, 0.0);

    for( unsigned long i = 0; i < regions->size(); i++ ) {
      currOutline = outlineSamples->at(i);

      if( currOutline == NULL )
        continue;

      //Find lloyds corners
      markerCorners.data.clear();
      quadFinder.computeCorners( currOutline, &markerCorners );

      //Add to cornerSamples if not outside padded bbox area /bbox neigborhood
      Vector offset( 2, 0.0 ), padVec( 2, 0.0);
      LongRangeList bbox;
      vector<long> bboxSide(2, 0);

      for( unsigned long j = 0; j < markerCorners.data.size(); j++ ) {
        currSample.copy( markerCorners.data.at(j) );

              //Directly reject all samples more away from a 3*3 window of adjacent boundingboxes
        bbox = ccp.getBoundingbox( regions->at(i) );
        for( long dim = 0; dim < 2; dim++)
          bboxSide.at(dim) = (bbox.vec[dim]).val.second - (bbox.vec[dim]).val.first ;

        padVec.set( 0, ccp.getPadding( regions->at(i) ) ); padVec.set( 1, ccp.getPadding( regions->at(i) ) );
        currSample -= padVec;
        offset.set(0, (bbox.vec[0]).val.first ); offset.set(1, (bbox.vec[1]).val.first );
        currSample += offset;

              //Reject corners
        if( isnan(currSample[0]) || isnan(currSample[1]) ||
          currSample[0] < (bbox.vec[0]).val.first - bboxSide.at(0) || currSample[1] < (bbox.vec[1]).val.first - bboxSide.at(1) ||
          currSample[0] > (bbox.vec[0]).val.second + bboxSide.at(0) || currSample[1] > (bbox.vec[1]).val.second + bboxSide.at(1))
          continue;

              //Add to global samplelist and store region index in seperate vector
        (cornerSamples->data).push_back( currSample );
        cornerSampleRegion.push_back( i );


              //Debug output
        if( debug == 2 ) {
          debugPoint.assign( currSample );
          debugPoint += drawOffset;
          debugPoint*= debugScaleFactor;
          writer.circle( debugPoint, circleSize / 2.0);
        }
      }
    }

//cerr << "debug:8" << endl;
    // group corners within sqrt(minArea) pixels of each other into clusters, since each saddle
    // point could have up to 4 nearby guessed points coming from each of the
    //adjacent squares.
    ClusterModifier clusterSolver;
    vector<long>* cornerClusters = new vector<long>( (cornerSamples->data).size() );
    long clusterNum = clusterSolver.clusterCutoff( sqrt( (double) minArea ), cornerSamples, cornerClusters);

    if( debug == 2 )
      cerr << "Found " << clusterNum << " initial corner-guesses."<<  endl;

    SampleVector* markerCornerGuesses = new SampleVector();
    clusterSolver.computeClusterMeans( cornerSamples, cornerClusters, markerCornerGuesses);

    //Now create for each cluster a list of marker regions it corresponds to
    vector< vector<long> > meansToSamplesMapping( (markerCornerGuesses->data).size() );
    for( long i = 0; i < cornerClusters->size(); i++ )
      (meansToSamplesMapping.at( cornerClusters->at(i) )).push_back( cornerSampleRegion.at(i) );

    //Cleanup
    cornerClusters->clear();
    delete cornerClusters;
    (cornerSamples->data).clear();
    delete cornerSamples;

    //Store points in region-indexed vector
    vector<CodeBlobProps*> regionCodeProps;
    for( unsigned long i = 0; i < regions->size(); i++ )
        regionCodeProps.push_back( new CodeBlobProps() );

    if( layout == 2 )
    {
      //Iterate over size of all found points and assign to according regionCodeProps
      for( unsigned long i = 0; i < (markerCornerGuesses->data).size(); i++ )
      {
          for( long j = 0; j < meansToSamplesMapping.at( i ).size(); j++ )
          {
              currSample.copy( (markerCornerGuesses->data).at(i) );
              ((regionCodeProps.at( (meansToSamplesMapping.at( i )).at(j) ))->corners).data.push_back( currSample );
          }
      }
    }

    //Debug output
    if( debug == 2 ) {
      //Output markerCornerGuesses in blue
      writer.getStyle().depth= 68;
      writer.getStyle().penColor = 10;
    }

    //Debug output
    if( debug == 2 ) {
      for( unsigned long i = 0; i < (markerCornerGuesses->data).size(); i++ ) {
        debugPoint.assign( (markerCornerGuesses->data).at(i) );
        debugPoint += drawOffset;
        debugPoint*= debugScaleFactor;
        writer.circle( debugPoint, circleSize / 1.5 );
      }
    }


    SampleVector bowtieIntersecGuesses;
    SampleVector bowtieIntersecSaddles;
    Vector tempVec(3, 1.0), currValueIndex(2, 0.0);   
    SaddlePointFinder<T> saddleFinder( this->array );
    vector<bool>* saddleDiverged = new vector<bool>();
    long goodMarkers = 0;
    long rejectedMarkers = 0;
    if( layout == 1 )
    {
      // find the precise saddle point coordinates in the original greyscale image
      saddleFinder.findSaddles( markerCornerGuesses, saddleDiverged);

      //Debug output of the final found saddles
      if( debug == 2 ) {
        for( unsigned long i = 0; i < (markerCornerGuesses->data).size(); i++ ) {
                //Continue if diverged
          if( saddleDiverged->at(i) )
            continue;

                //Get sample
          debugPoint.assign( (markerCornerGuesses->data).at(i)  );
          debugPoint += drawOffset;
          debugPoint*= debugScaleFactor;
          writer.circle( debugPoint, circleSize * 2.0 );
        }
      }

      //Iterate over size of all found points and assign to according regionCodeProps
      //Create saddlesGuesses without repetitions of same saddles
      SampleVector singleSaddles;
      for( unsigned long i = 0; i < (markerCornerGuesses->data).size(); i++ ) {
        //Continue if diverged
        if( saddleDiverged->at(i) )
          continue;

        for( long j = 0; j < meansToSamplesMapping.at( i ).size(); j++ ) {
          currSample.copy( (markerCornerGuesses->data).at(i) );
          ((regionCodeProps.at( (meansToSamplesMapping.at( i )).at(j) ))->corners).data.push_back( currSample );
        }

        //Add to saddleGuesses without repetitions
        currSample.copy( (markerCornerGuesses->data).at(i) );
        bool alreadyfound = false;
        for( long j = 0; j < singleSaddles.data.size(); j++ ) {
          if( dist( currSample, singleSaddles.data.at(j) ) < 0.1  ) {
            alreadyfound = true;
            break;
          }
        }
        if( !alreadyfound )
          singleSaddles.data.push_back( currSample );
      }

      //Destroy all regions without exactly 4 points
      //Simultaneously count good markers ( non-diverged samples >= 4 )
      for( unsigned int currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ ) {
        //Debug output
        //cerr << "Corners found in region: " << ((regionCodeProps.at( currRegion ))->corners).data.size() << endl;
        if( ((regionCodeProps.at( currRegion ))->corners).data.size() != 4  )
        {
          ((regionCodeProps.at( currRegion ))->corners).data.clear();
          rejectedMarkers++;
        }
        else
          goodMarkers++;
      }
    }else if( layout == 2 )
    {
      //Fit homography to find bowtie intersection guesses
      vector<Vector> currTheta; 
      SampleVector quadSquare;
      //Initialize
      for( unsigned int i = 0; i < 4; i++ )
      {
          tempVec.copy( oneVec ); quadSquare.data.push_back( tempVec );
          currValueIndex.copy( zeroVec ); currTheta.push_back( currValueIndex );
      }

      //Homography matrix
      Matrix H(3,3);

      //iterate over all regions
      //Destroy all regions without exactly 4 points
      //Simultaneously count good markers ( non-diverged samples >= 4 )
      for( unsigned long currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ )
      {
          //Continue if not enough saddles were found for homography calculation
          if( ((regionCodeProps.at( currRegion ))->corners).data.size() != 4 )
              continue;
          
          // sort cornerguesses into clockwise order
          Vector mean(2, 0.0);
          for( unsigned int i = 0; i < 4; i++ )
              mean += ((regionCodeProps.at( currRegion ))->corners).data.at(i);
          mean /= 4.0;

          for( unsigned int i = 0; i < 4; i++ )
          {
              (currTheta.at(i))[0] = atan2( (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(1) - mean.get(1),
                                       (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(0) - mean.get(0));
              if( (currTheta.at(i))[0] < 0 ) 
                    (currTheta.at(i))[0] += 2 * M_PI;

              (currTheta.at(i))[1] = (double)i; //Index for sorting
          }
          std::sort( currTheta.begin(), currTheta.end(), CALTagPattern::compareIndexPairs );

          //create a clockwise sorted quadsquare in homogeneous 2D coords
          int index;
          for( int i = 3; i >= 0; i-- )
          {
              index = (unsigned int)((currTheta.at(i)).get(1));
              quadSquare.data.at(i).set( 0, (((regionCodeProps.at( currRegion ))->corners).data.at(index)).get(0) );
              quadSquare.data.at(i).set( 1, (((regionCodeProps.at( currRegion ))->corners).data.at(index)).get(1) );
              quadSquare.data.at(i).set( 2, 1.0 );
          }

          //Compute homography
          bowtieIntersecGuesses.data.clear();
          homography2D(unitSquare, &quadSquare, &H);
          homotrans( bowtieIntersections, &H, &bowtieIntersecGuesses );

          //Debug output
          if( debug == 2 )
          {
              //Output triangle intersection positions in red
              writer.getStyle().penColor = 4;
              writer.getStyle().depth= 68;

              for( unsigned int i = 0; i < bowtieIntersecGuesses.data.size(); i++ )
              {
                  debugPoint.assign( bowtieIntersecGuesses.data.at(i).getSubVector(0,1) );
                  debugPoint += drawOffset;
                  debugPoint*= debugScaleFactor;
                  if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
                        writer.circle( debugPoint, circleSize );
              }
          }

          //get rid of third coordinate for saddlefinder
          bowtieIntersecSaddles.data.clear();
          for( unsigned int i = 0; i < bowtieIntersecGuesses.data.size(); i++ )
          {
              currSample.copy( zeroVec );
              currSample.assign( bowtieIntersecGuesses.data.at(i).getSubVector(0,1) );
              bowtieIntersecSaddles.data.push_back( currSample );
          }

          //Now try to find saddles and insert as new regionprop corners if none of them diverged
          saddleDiverged->clear();
          saddleFinder.findSaddles( &bowtieIntersecSaddles, saddleDiverged);

          //Refill regiondata with found saddles
          ((regionCodeProps.at( currRegion ))->corners).data.clear();
          for( unsigned int i = 0; i < bowtieIntersecSaddles.data.size(); i++ )
          {
              if( saddleDiverged->at(i) )
              {
                  ((regionCodeProps.at( currRegion ))->corners).data.clear();
                  rejectedMarkers++;
                  break;
              }
              currSample.copy( zeroVec );
              currSample.assign( bowtieIntersecSaddles.data.at(i) );
              ((regionCodeProps.at( currRegion ))->corners).data.push_back(currSample) ;
          }

          if( ((regionCodeProps.at( currRegion ))->corners).data.size() == 4 )
              goodMarkers++;

          //Debug output of the final found saddles
          if( debug == 2 )
          {   
              //Output final found saddles in blue
              writer.getStyle().depth= 67;
              writer.getStyle().penColor = 10;

              for( unsigned int i = 0; i < ((regionCodeProps.at( currRegion ))->corners).data.size(); i++ )
              {
                  //Get sample
                  debugPoint.assign( ((regionCodeProps.at( currRegion ))->corners).data.at(i)  );
                  debugPoint += drawOffset;
                  debugPoint*= debugScaleFactor;
                  if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
                      writer.circle( debugPoint, circleSize * 2.0 );
              }
          }
  
      }
    }

    //Cleanup
    (markerCornerGuesses->data).clear();
    delete markerCornerGuesses;
    saddleDiverged->clear();
    delete saddleDiverged;
    outlineSamples->clear();
    delete outlineSamples;

    if( debug == 1 || debug == 2 )
    {
        cerr << "Good markers found: " << goodMarkers << endl;
        cerr << "Markers rejected: " << rejectedMarkers << endl;
    }

//cerr << "debug:10" << endl;
    // early exit if there are not enough good markers
    if( (double)goodMarkers < ((double)numMarkers)/8.0) {
      if( debug == 2 )
        cerr << "Not enough markers were detected !" << endl;
      regionCodeProps.clear();
      if(debug == 2 || debug == 3) {
              // image in the background
        writer.getStyle().depth= 70;
        writer.image( greyscaleImageFilename.c_str(), arrayDim.vec[0]*debugScaleFactor, arrayDim.vec[1]*debugScaleFactor );

              // #### End of Debug output
        writer.close();
      }

      //[brad] changed successfulTracked from bool* to bool
      //if( successfulTracked != NULL )
      successfulTracked = false;
      return false;
    }

    //set up linear point sampler for code sample interpolation
    LinearPointSampler<T> linInterpolSampler( arrayDim );
    boundary = Clamp;

    //get code sample data from array
    T* arrayData = &((*((*this->array)[0]))[0]);
    unsigned long numArrayEntries = 1;
    for( int d= arrayDim.vec.size()-1 ; d>= 0 ; d-- )
      numArrayEntries*= arrayDim.vec[d];

    //Compute means for adaptive thresholding
    vector<double> meanIntensities( regions->size(), 0.0 );
    //Iterate over regions
    LongRangeList bbox;
    for( unsigned long currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ ) {
      //Continue if not enough saddles were found
      if( ((regionCodeProps.at( currRegion ))->corners).data.size() != 4 )
        continue;

      //Get boundingbox
      bbox = ccp.getBoundingbox( regions->at(currRegion) );
      vector<long> bboxSide(2, 0);
      for( long dim = 0; dim < 2; dim++)
        bboxSide.at(dim) = (bbox.vec[dim]).val.second - (bbox.vec[dim]).val.first ;

      //iterate over 9 "adjacent" markers in original greyscale image
      unsigned long srcPos;
      unsigned long setPixels = 0;
      for( unsigned long y = (bbox.vec[1]).val.first - bboxSide.at(1) ; y <= (bbox.vec[1]).val.second + bboxSide.at(1)  ; y++ ) {
        for( unsigned long x = (bbox.vec[0]).val.first - bboxSide.at(0) ; x <= (bbox.vec[0]).val.second + bboxSide.at(0) ; x++ ) {
          srcPos = y * arrayDim.vec[0] + x;
          if( srcPos >= 0 && srcPos < numArrayEntries ) {
            meanIntensities.at(currRegion) += arrayData[srcPos];
            setPixels++;
          }
        }
      }
      meanIntensities.at(currRegion) /= (double)setPixels;
    }

//cerr << "debug:11" << endl;
    //given the corner points of each marker, we can sample the interior to read the code
    //pattern. Adaptive thresholding is used to account for shadows.
    vector<Vector> currTheta;
    SampleVector quadSquare;

    //Initialize
    for( unsigned int i = 0; i < 4; i++ ) {
      tempVec.copy( oneVec ); quadSquare.data.push_back( tempVec );
      currValueIndex.copy( zeroVec ); currTheta.push_back( currValueIndex );
    }
    Matrix H(3,3);
    SampleVector foundCodeSamplePos;
    unsigned int validCodesFound = 0;

    //try to read all marker codes:
    for( unsigned long currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ ) {
      //Continue if not enough saddles were found for homography calculation
      if( ((regionCodeProps.at( currRegion ))->corners).data.size() != 4 )
        continue;

      if( layout == 1 ) //Original pattern
      {
        // sort saddles into clockwise order
        Vector mean(2, 0.0);
        for( unsigned int i = 0; i < 4; i++ )
          mean += ((regionCodeProps.at( currRegion ))->corners).data.at(i);
        mean /= 4.0;

        for( unsigned int i = 0; i < 4; i++ ) {
          (currTheta.at(i))[0] = atan2( (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(1) - mean.get(1),
            (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(0) - mean.get(0));
          if( (currTheta.at(i))[0] < 0 )
            (currTheta.at(i))[0] += 2 * M_PI;

                //Index for sorting
          (currTheta.at(i))[1] = (double)i;
        }
        std::sort( currTheta.begin(), currTheta.end(), CALTagPattern::compareIndexPairs );

        //create a clockwise sorted quadsquare in homogeneous 2D coords and also sort regionCodeProps vectors
        int index;
        for( int i = 3; i >= 0; i-- ) {
          index = (unsigned int)((currTheta.at(i)).get(1));
          quadSquare.data.at(i).set( 0, (((regionCodeProps.at( currRegion ))->corners).data.at(index)).get(0) );
          quadSquare.data.at(i).set( 1, (((regionCodeProps.at( currRegion ))->corners).data.at(index)).get(1) );
          quadSquare.data.at(i).set( 2, 1.0 );
        }

        for( unsigned int i = 0; i < 4; i++ ) {
          (((regionCodeProps.at( currRegion ))->corners).data.at(i)).set( 0, quadSquare.data.at(i).get(0) );
          (((regionCodeProps.at( currRegion ))->corners).data.at(i)).set( 1, quadSquare.data.at(i).get(1) );
        }
      }else if( layout == 2 ) //rotated pattern
      {
        //Compute homography ( regionCodeProps corners are allready in clockwise sorting )
        for( int i = 0; i < 4; i++ )
        {
            quadSquare.data.at(i).set( 0, (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(0) );
            quadSquare.data.at(i).set( 1, (((regionCodeProps.at( currRegion ))->corners).data.at(i)).get(1) );
            quadSquare.data.at(i).set( 2, 1.0 );
        }
      }

      homography2D(unitSquare, &quadSquare, &H);
      homotrans( codeBlobSamplePos, &H, &foundCodeSamplePos );

      //Save into binary code pattern and output debug
      //Debug output
      if( debug == 2 || debug == 3  ) {
        //Output codeblob positions
        writer.getStyle().depth= 66;
        writer.getStyle().penColor = 2;
        writer.getStyle().fillColor = 2;
        writer.getStyle().areaFill = 20;
      }
      for( unsigned int i = 0; i < foundCodeSamplePos.data.size(); i++ ) {
        linInterpolSampler.setSampleLocation(foundCodeSamplePos.data.at(i).getSubVector(0,1) , boundary);
        code->at(i) = ((linInterpolSampler.getSample( arrayData )) > meanIntensities.at(currRegion));

              //Debug output
        if( debug  == 2 || debug == 3 ) {
          debugPoint.assign( foundCodeSamplePos.data.at(i).getSubVector(0,1) );
          debugPoint += drawOffset;
          debugPoint*= debugScaleFactor;
          if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
            writer.circle( debugPoint, circleSize * 0.2 );
        }
      }

      //Debug output
      if( debug == 2 || debug == 3 ) {
              //Reset xfig output props
        writer.getStyle().penColor = 4;
        writer.getStyle().fillColor = -1;
        writer.getStyle().areaFill = -1;
      }

      //check each of the 4 possible orientations (ignoring reflections) to
      //see if they are valid codes
      unsigned int ID, mult, crc;
      unsigned int maxid = ( 2 << (idBits - 1) );
      Vector orientVec(2, 0.0);
//cerr << "Checking code: " << endl;
//cerr << (*code)[15] << (*code)[14] << (*code)[13] << (*code)[12] << endl
//		 << (*code)[11] << (*code)[10] << (*code)[9] << (*code)[8] << endl
//		 << (*code)[7] << (*code)[6] << (*code)[5] << (*code)[4] << endl
//		 << (*code)[3] << (*code)[2] << (*code)[1] << (*code)[0] << endl;
      for( int orient = 0; orient < 4; orient++ ) {
              // read the id from the image and convert to decimal. Note that
              // each binary 2D code pattern is constructed so that the code
              // string can be read out in columnwise ordering, first the ID
              // and then the CRC, each in most-significant-bit-first order.
        mult = 1; ID = 0;
        for( int i = idBits - 1; i >= 0; i-- ) {
          ID += (int)code->at(i) * mult;
          mult *= 2;
        }

        if( ( ID > 0 ) && ( ID <= maxid ) ) {
              // read the crc from the image and convert to decimal
          mult = 1; crc = 0;
          for( int i = code->size() - 1; i >= code->size() - cBits; i-- ) {
            crc += (int)code->at(i) * mult;
            mult *= 2;
          }

          // check if image crc matches the expected value
          //[brad] replacing precomputed crcTable with check via CRC class
          //if( crc == crcTable[ID-1] )
          CRCCode codec( idBits, cBits, genPoly );
//cerr << "Validating ID=" << ID << ", CRC="<<crc << endl;
          if( codec.validate( (ID<<cBits)|crc ) ) {
//cerr << "debug:validate==true" << endl;
            (regionCodeProps.at( currRegion ))->valid = true;
            (regionCodeProps.at( currRegion ))->ID = ID;
              // orientation as vector from first to second corner
            (regionCodeProps.at( currRegion ))->firstCorner = (orient + 1) % 4;
            orientVec.assign( ((regionCodeProps.at( currRegion ))->corners).data.at( (orient + 2) % 4 ) );
            orientVec -= ((regionCodeProps.at( currRegion ))->corners).data.at( (orient + 1) % 4 ) ;

              //Debug output of orientation
            if( debug == 2 ) {
              Vector orientVec1(2, 0.0), orientVec2(2, 0.0);
              orientVec1.assign( ((regionCodeProps.at( currRegion ))->corners).data.at( (orient + 1) % 4 ) );
              orientVec2.assign( orientVec ); orientVec2 *= 0.5;
              orientVec2 += orientVec1;
              orientVec1 += drawOffset;
              orientVec1*= debugScaleFactor;
              orientVec2 += drawOffset;
              orientVec2*= debugScaleFactor;

              writer.lineSegment( orientVec1, orientVec2 );
            }

            orientVec.normalize();
            ((regionCodeProps.at( currRegion ))->orientation).assign( orientVec );

              /*Debug output
            cerr << "Valid Code found: ID: " << ID << " CRC: " << crc << endl;*/
            validCodesFound++;
            break;
          }
          else {
            rotCode90();
          }
        }
      }
      /*Debug output ID and code pattern found
      cerr << "ID " << (regionCodeProps.at( currRegion ))->ID << endl;
      for( int y = 0; y < resCode; y++ )
      {
          for( int x = 0; x < resCode; x++ )
          {
              cerr << code->at( y* resCode +x) << " ";
          }
          cerr << endl;
      }*/
    }
    //Debug output
    if( debug == 1 || debug == 2 )
      cerr << "Found " << validCodesFound << " valid codes." << endl;

//cerr << "debug:12" << endl;
    //Cleanup
    delete regions;
    foundCodeSamplePos.data.clear();
    quadSquare.data.clear();
    currTheta.clear();


    // filter out any squares that have a very different orientation to the
    // majority (can occur on random background patterns, or else marker squares
    // where a few code bits are incorrectly read causing the entire pattern to
    // mistaken for a different code at another orientation).
    // First get median from sorting of angle to x-axis
    Vector medianOrient(2, 0.0);
    vector<Vector> orientations;
    for( unsigned int i = 0; i < regionCodeProps.size(); i++ ) {
      //Continue if no valid code was found
      if( !((regionCodeProps.at( i ))->valid ) )
        continue;

      currSample.copy( zeroVec );
      currSample.set(0, acos( ((regionCodeProps.at( i ))->orientation).get(0) ) );
      currSample.set(1, i);
      orientations.push_back(currSample);
    }

    //Sort out any region with angle between its orientation and medianOrient > angleThreshold (default 90 degrees)
    unsigned int outsortedOrients = 0;
              //Do outsorting only if something can be outsorted
    if( orientations.size() > 1 ) {
      std::sort( orientations.begin(), orientations.end(), compareIndexPairs );
      unsigned int medianIndex =  (unsigned int)(round((double)(orientations.size())/2.0));
      medianOrient.assign( (regionCodeProps.at( (unsigned int)(round((orientations.at(medianIndex)).get(1))) ))->orientation );

      //Sort out orientations
      for( unsigned int i = 0; i < regionCodeProps.size(); i++ ) {
              //Continue if no valid code was found
        if( !((regionCodeProps.at( i ))->valid ) )
          continue;

        if( acos( dot((regionCodeProps.at( i ))->orientation, medianOrient) ) > (angleThreshold/ 180.0 * M_PI) ) {
          (regionCodeProps.at( i ))->valid = false;
          ((regionCodeProps.at( i ))->corners).data.clear();
          outsortedOrients++;
        }
      }
    }
    if( debug == 2 )
      cerr << "Outsorted Orientations: " << outsortedOrients << endl;

    if( validCodesFound - outsortedOrients <= 0 ) {
      if(debug == 2 || debug == 3 ) {
              // image in the background
        writer.getStyle().depth= 70;
        writer.image( greyscaleImageFilename.c_str(), arrayDim.vec[0]*debugScaleFactor, arrayDim.vec[1]*debugScaleFactor );

              // #### End of Debug output
        writer.close();
      }

      //if( successfulTracked != NULL )
      successfulTracked = false;

      return false;
    }

//cerr << "debug:13" << endl;
    //Now try to detect the whole calibration grid (with the missed markers
    //A samplevector for the whole calibration grid, accessed as array[numRows][numCols]
    ((this->imageCalibPoints)->data).clear();
    Vector badVector(3, -1.0); badVector.set(2, 1.0);
    for( unsigned int i = 0; i < (rows+1)*(columns+1); i++ ) {
      tempVec.copy( badVector );
      ((this->imageCalibPoints)->data).push_back( tempVec );
    }

    //Get the coordinates of the grid corners from the data
    unsigned int currFirstCorner;
    for( unsigned long currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ ) {
      //Continue if no valid code was found
      if( !((regionCodeProps.at( currRegion ))->valid ) )
        continue;

      //Try to find position of ID in grid
      unsigned int row, col;
      bool found = false;
      for( row = 0; row < rows ; row++ ) {
        for( col = 0; col < columns ; col++ ) {
              //[brad] replace hardcoded codeTable[9][8] array with the ids property
              //which is a CoordinateVector which is a std::vector<int> which gets
              //filled during the call to configure(). Note that we may have to do some
              //nasty transposes and mirror flips soon to get the coordinate frames to
              //align. And note that ids includes the crc, but we only want the id here.
              //if( codeTable[row][col] == ((regionCodeProps.at( currRegion ))->ID) )
          if( (ids.vec[row*columns+col]>>cBits) == ((regionCodeProps.at( currRegion ))->ID) ) {
            found = true;
            break;
          }
        }
        if( found )
          break;
      }
      //Continue if not found
      if( !found )
        continue;

      //Save gridpositions
      currFirstCorner = (regionCodeProps.at( currRegion ))->firstCorner;

      /*Debug if marker points are scattered right into imageCalibPoints
      if( row == 3 && col == 4 )
      {
		writer.getStyle().depth= 30;
		writer.getStyle().penColor = 5;
		for( int i = 0; i < 4; i++ )
		{
			debugPoint.set(0, ((regionCodeProps.at( currRegion ))->corners).data.at((currFirstCorner + i) % 4).get(0) );
			debugPoint.set(1, ((regionCodeProps.at( currRegion ))->corners).data.at((currFirstCorner + i)% 4).get(1) );
			debugPoint += drawOffset;
			debugPoint*= debugScaleFactor;
			cerr << debugPoint.get(0) << " " << debugPoint.get(1) << endl;
			if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
			    writer.circle( debugPoint, circleSize * (1+i));
		}
		      cerr << "endmarker";
      }*/

      for( int i = 0; i < 2 ; i++ ) {
        (((this->imageCalibPoints)->data).at((row+1)*(columns+1) + col)).
          set( i, ((regionCodeProps.at( currRegion ))->corners).data.at(currFirstCorner % 4).get(i)  );

        (((this->imageCalibPoints)->data).at((row+1)*(columns+1) + (col+1))).
          set( i, ((regionCodeProps.at( currRegion ))->corners).data.at((currFirstCorner + 1) % 4).get(i)  );

        (((this->imageCalibPoints)->data).at(row*(columns+1) + (col+1))).
          set( i, ((regionCodeProps.at( currRegion ))->corners).data.at((currFirstCorner + 2) % 4).get(i)  );

        (((this->imageCalibPoints)->data).at(row*(columns+1) + col)).
          set( i, ((regionCodeProps.at( currRegion ))->corners).data.at((currFirstCorner + 3) % 4).get(i)  );
      }
    }


//cerr << "debug:14" << endl;
    // fill in missing points
    // Construct a ransac homography with all good points and
    // then compute positions of the missed points with that homography
    SampleVector goodGridPoints, goodUnitGridPoints, missedUnitPoints, distortionCorrectedPoints;
    unsigned int rowcolIndex;
    for( unsigned int row = 0; row < rows+1 ; row++ ) 
    {
      for( unsigned int col = 0; col < columns+1 ; col++ ) 
      {
        rowcolIndex = row*(columns+1) + col;
              // If bad vector, add missed unit gridpoint
        if( ((this->imageCalibPoints)->data).at( rowcolIndex ).get(0) < 0.0 || ((this->imageCalibPoints)->data).at( rowcolIndex ).get(0) >= arrayDim.vec[0] ||
	    ((this->imageCalibPoints)->data).at( rowcolIndex ).get(1) < 0.0  || ((this->imageCalibPoints)->data).at( rowcolIndex ).get(1) >= arrayDim.vec[1] ) 
	  {
              //Missed unit grid point
          tempVec.copy( oneVec );
          tempVec.set( 0, (double)row + 1.0 );
          tempVec.set( 1, (double)col + 1.0);
          missedUnitPoints.data.push_back( tempVec );
        }else
        {
          //Gridpoint
          currSample.copy( zeroVec );
	  currSample.set(0, ((this->imageCalibPoints)->data).at( rowcolIndex ).get(0) );
	  currSample.set(1, ((this->imageCalibPoints)->data).at( rowcolIndex ).get(1) );
          goodGridPoints.data.push_back( currSample );

	  //Distortion vector
	  currSample.copy( zeroVec ); 
	  currSample.set(0, ((this->imageCalibPoints)->data).at( rowcolIndex ).get(0) );
	  currSample.set(1, ((this->imageCalibPoints)->data).at( rowcolIndex ).get(1) );
	  distortionCorrectedPoints.data.push_back( currSample );

          //Good unit grid point
          tempVec.copy( oneVec );
          tempVec.set( 0, (double)row + 1.0 );
          tempVec.set( 1, (double)col + 1.0);
          goodUnitGridPoints.data.push_back( tempVec );
        }
      }
    }

//## Distortion correction:

      //Debug output of found saddles
      if( debug == 3 )
      {
	  writer.getStyle().depth= 50;
	  writer.getStyle().penColor = 4;

	  for(unsigned int i = 0; i < distortionCorrectedPoints.data.size(); i++)
	  {
	      debugPoint.set(0,  distortionCorrectedPoints.data.at( i ).get(0) );
	      debugPoint.set(1,  distortionCorrectedPoints.data.at( i ).get(1) );
	      debugPoint += drawOffset;
	      debugPoint*= debugScaleFactor;
	      if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
	      writer.circle( debugPoint, circleSize );
	  }
      }

      //####Handling of distortion
      //Parameters
      double precision = 0.0001;	//Precision of distortion correction
      double diagScale = 1.0;

      //Apply distortion correction
      Vector imgCenter(2, 0.0); 
      imgCenter.set(0, ((double)arrayDim.vec[0]) / 2.0 ); 
      imgCenter.set(1, ((double)arrayDim.vec[1]) / 2.0 );
      double halfDiagonal = imgCenter.norm(2) * diagScale;
      
      //Search for an r by computing r's over all rows and columns 
      //and finally averaging these r's
      SampleVector collinearSamples;
      vector<double> rVector; vector<double> rVecWeights;
      for( unsigned int row = 0, col = 0; (row < rows + 1) && (col < columns + 1); row++, col++)
      {
	      //Process Col
	      collinearSamples.data.clear();
	      for( unsigned int tmpRow = 0; tmpRow < rows+1 ; tmpRow++ )
	      {
			unsigned int index = tmpRow*(columns+1) + col;
			// If bad vector, add missed unit gridpoint
			if( ((this->imageCalibPoints)->data).at( index ).get(0) > 0.0 && ((this->imageCalibPoints)->data).at( index ).get(0) < arrayDim.vec[0] && ((this->imageCalibPoints)->data).at( index ).get(1) > 0.0  && ((this->imageCalibPoints)->data).at( index ).get(1) < arrayDim.vec[1] )
		{
			    currSample.copy( zeroVec ); 
			    currSample.set(0, ((this->imageCalibPoints)->data).at( index ).get(0) );
			    currSample.set(1, ((this->imageCalibPoints)->data).at( index ).get(1) );
			    collinearSamples.data.push_back( currSample );
			}
	      }
	      if( collinearSamples.data.size() > 4 ) //Only process lines with at least 4 points
			rVector.push_back( findUndistortR( &collinearSamples, &imgCenter, halfDiagonal, precision ) );

	      //Process Row
	      collinearSamples.data.clear();
	      for( unsigned int tmpCol = 0; tmpCol < columns+1 ; tmpCol++ )
	      {
			unsigned int index = row*(columns+1) + tmpCol;
			// If bad vector, add missed unit gridpoint
			if( ((this->imageCalibPoints)->data).at( index ).get(0) > 0.0 && ((this->imageCalibPoints)->data).at( index ).get(0) < arrayDim.vec[0] && ((this->imageCalibPoints)->data).at( index ).get(1) > 0.0  && ((this->imageCalibPoints)->data).at( index ).get(1) < arrayDim.vec[1] )
			{
			    currSample.copy( zeroVec ); 
			    currSample.set(0, ((this->imageCalibPoints)->data).at( index ).get(0) );
			    currSample.set(1, ((this->imageCalibPoints)->data).at( index ).get(1) );
			    collinearSamples.data.push_back( currSample );
			}
	      }
	      if( collinearSamples.data.size() > 4 ) //Only process lines with at least 4 points
			rVector.push_back( findUndistortR( &collinearSamples, &imgCenter, halfDiagonal, precision ) );
      }

      //Average r
      double r = 0.0;
      if( rVector.size() > 0 )
      {
	  for( unsigned int i = 0; i < rVector.size(); i++ )
		r += rVector.at(i);

	  r /= ((double)rVector.size());
      }

      if( debug == 2 )
	  cerr << "Undistortion R is: " << r << endl;

      //Do undistortion if an significant R was found
      bool applyDistortionCorrection = fabs(r) > precision; //if false do not apply correction at all
      if( applyDistortionCorrection ) 
      {
	  correctDistortion( &distortionCorrectedPoints, &imgCenter, halfDiagonal, r);

	  //Recompute the undistorted length of the diagonal
	  Vector diag(2, 0.0);
	  correctDistortion( &diag, &imgCenter, halfDiagonal, r );
	  halfDiagonal = diag.norm(2) * diagScale;
      }

      /*###DEBUGGING PURPOSE: Invert the undistortion operation for debugging of distortion correction
      if( applyDistortionCorrection ) 
      {
	    double r_inv = findInversionR( &distortionCorrectedPoints, &goodGridPoints, &imgCenter, halfDiagonal, precision );

	    //Now invert Undistortion for a copy of distortionCorrectedPoints
	    SampleVector debugDistortionCorrectedPoints;
	    for( unsigned int i = 0; i < distortionCorrectedPoints.data.size(); i++ )
	    {
		  currSample.copy( zeroVec ); 
		  currSample.assign( distortionCorrectedPoints.data.at( i ) );
		  debugDistortionCorrectedPoints.data.push_back( currSample );
	    }

	    correctDistortion( &debugDistortionCorrectedPoints, &imgCenter, halfDiagonal, r_inv);

	    //Debug output of debugDistortionCorrectedPoints
	    if( debug == 3 )
	    {
		writer.getStyle().depth= 49;
		writer.getStyle().penColor = 5;

		for(unsigned int i = 0; i < debugDistortionCorrectedPoints.data.size(); i++)
		{
		    debugPoint.assign( debugDistortionCorrectedPoints.data.at( i ) );
		    debugPoint += drawOffset;
		    debugPoint*= debugScaleFactor;
		    cerr << debugPoint.get(0) << " " << debugPoint.get(1) << endl;
		    if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
		    writer.circle( debugPoint, circleSize );
		}
	    }
      }*/

      //Prepare distpoints for homography
      SampleVector homogeneDistortionCorrectedPoints;
      for( unsigned int i = 0; i < distortionCorrectedPoints.data.size(); i++ )
      {
	    tempVec.copy( oneVec ); 
	    tempVec.set(0, distortionCorrectedPoints.data.at( i ).get(0) );
	    tempVec.set(1, distortionCorrectedPoints.data.at( i ).get(1) );
	    homogeneDistortionCorrectedPoints.data.push_back( tempVec );
      }

//## End of distortion correction

    //Fit ransac homography with all good (undistorted) points
    Matrix ransacH(3,3);
    vector<unsigned long>* inliers = new vector<unsigned long>();
    bool ransacOK = ransacFitHomography(&goodUnitGridPoints, &homogeneDistortionCorrectedPoints, &ransacH, inliers, 0.01 );

    //Debug output
    if( ! ransacOK ) {
      if(debug == 2  || debug == 3) {
        cerr << "Could not calculate Ransac homography." << endl;
              // image in the background
        writer.getStyle().depth= 70;
        writer.image( greyscaleImageFilename.c_str(), arrayDim.vec[0]*debugScaleFactor, arrayDim.vec[1]*debugScaleFactor );

              // #### End of Debug output
        writer.close();
      }

      //if( successfulTracked != NULL )
      successfulTracked = false;

      return false;
    }

    //Now refit outliers and originally bad points
    //First add oultiers
    unsigned int outliersNum = 0;
    for( unsigned int i = 0; i < goodUnitGridPoints.data.size() && ransacOK; i++ ) {
      bool foundInlier = false;
      for( unsigned int j = 0; j < inliers->size(); j++ ) {
        if( inliers->at(j) == i ) {
          foundInlier = true;
          break;
        }
      }
      
      if( ! foundInlier ) {
              //Add missed unit grid point
        tempVec.copy( goodUnitGridPoints.data.at(i) );
        missedUnitPoints.data.push_back( tempVec );
        outliersNum++;
      }
    }
    delete inliers;
    if( debug == 2 )
      cerr << "Outliers added: " << outliersNum << endl;

//cerr << "debug:16" << endl;
    //Find missed Samples guesses in image and remove samples outside of image border
    SampleVector missedSampleGuesses;
    if( ransacOK )
      homotrans( &missedUnitPoints, &ransacH, &missedSampleGuesses );

    SampleVector missedSampleGuesses2d, missedSampleOriginalPos;
    for( unsigned int i = 0; i < missedSampleGuesses.data.size() && ransacOK; i++ ) {
      // remove the points that map outside the image borders from consideration
      /* and keep track of original grid position
      if(  missedSampleGuesses.data.at(i).get(0) < 5   ||
        missedSampleGuesses.data.at(i).get(0) > arrayDim.vec[0]-5  ||
        missedSampleGuesses.data.at(i).get(1) < 5   ||
        missedSampleGuesses.data.at(i).get(1) > arrayDim.vec[1]-5 )
        continue;*/ //Has to be disabled for proper distortion handling.

      currSample.copy( zeroVec );
      currSample.set( 0, missedSampleGuesses.data.at(i).get(0) );
      currSample.set( 1, missedSampleGuesses.data.at(i).get(1) );
      missedSampleGuesses2d.data.push_back( currSample );
      currSample.copy( zeroVec );
      currSample.set( 0, missedUnitPoints.data.at(i).get(0) - 1.0 );
      currSample.set( 1, missedUnitPoints.data.at(i).get(1) - 1.0 );
      missedSampleOriginalPos.data.push_back( currSample );
    }

    //Invert the undistortion operation now
    if( applyDistortionCorrection ) 
    {
	  double r_inv = findInversionR( &distortionCorrectedPoints, &goodGridPoints, &imgCenter, halfDiagonal, precision );
	  if( debug == 2 )
	      cerr << "Inversion R is: " << r_inv << endl;
     
	  //Now invert Undistortion
	  correctDistortion( &missedSampleGuesses2d, &imgCenter, halfDiagonal, r_inv);
    }

    /*Debug output of saddleguesses
    if( debug == 3 )
    {
	//Guesses in purple
	writer.getStyle().penColor = 2;
	writer.getStyle().depth= 40;

	debugPoint.assign( imgCenter );
	debugPoint += drawOffset;
	debugPoint*= debugScaleFactor;
	writer.circle( debugPoint, circleSize* 2.0 );

	writer.getStyle().penColor = 5;

	for(unsigned int i = 0; i < missedSampleGuesses2d.data.size(); i++)
	{
	    debugPoint.assign( missedSampleGuesses2d.data.at( i ) );
	    debugPoint += drawOffset;
	    debugPoint*= debugScaleFactor;
	    if( ! isnan(debugPoint[0]) && ! isnan(debugPoint[1]) )
	    writer.circle( debugPoint, circleSize* 1.7 );
	}
    }*/

    // find the precise saddle points of the missed points
    vector<bool>* missedSaddleDiverged = new vector<bool>();
    saddleFinder.findSaddles( &missedSampleGuesses2d, missedSaddleDiverged);

//Felix: has to be tested more
//###Validate new missed saddlepoints by sampling circles around each possible saddle
    for( unsigned int i = 0; i < missedSampleGuesses2d.data.size(); i++) 
    {
	if( missedSaddleDiverged->at(i) )
	    continue;

	double maxRadius = 0.0;
	long numRadii = 0;
	for( unsigned long currRegion = 0; currRegion < regionCodeProps.size(); currRegion++ ) 
	{
	    //Continue if no valid code was found
	    if( !((regionCodeProps.at( currRegion ))->valid ) )
	      continue;

	    // Calc mean
	    Vector mean(2, 0.0);
	    for( unsigned int i = 0; i < 4; i++ )
		mean += ((regionCodeProps.at( currRegion ))->corners).data.at(i);
	    mean /= 4.0;

	    for( unsigned int i = 1; i < 4; i++ ) 
	    {
		maxRadius += dist( ((regionCodeProps.at( currRegion ))->corners).data.at(i), mean);
		numRadii++;
	    }
	}

	//Average
	if( numRadii > 0 ) 
	    maxRadius /= (double)numRadii;


	if( layout == 1 ) //Original pattern
	    maxRadius /= 5.0; //Seems to be good heuristic although  maxRadius /= 2.0 would be the maximal possible number for layout 1
	else if( layout == 2 ) //Rotated layout
	    maxRadius /= 2.0;

	if( maxRadius < 3.0 )
	    maxRadius = 3.0;

	//Debug output
	//cerr << "MaxRad: " << maxRadius << endl;

	//Now try to validate all saddle
	if( !validateSaddle(this->array, missedSampleGuesses2d.data.at(i), 2, maxRadius/2.0, maxRadius) )
	  {

	    missedSaddleDiverged->at(i) = true;

	  }

    }
//####
//cerr << "debug:17" << endl;

    //Debug output
    unsigned int recoveredSaddles = 0;
    for( unsigned int i = 0; i < missedSaddleDiverged->size(); i++)
      recoveredSaddles += (missedSaddleDiverged->at(i))? 0 : 1;

    if( debug == 2 )
      cerr << "Tryed to recover " << missedSampleGuesses2d.data.size() << " missing corners. "
        << recoveredSaddles << " recovered." << endl;

    //Store newly found saddles in original grid:
    for( unsigned int i = 0; i < missedSampleGuesses2d.data.size(); i++) {
      //If saddle points were found, add them to the points matrix
      if( missedSaddleDiverged->at(i) )
        continue;

      //Get index from missedSampleOriginalPos position vector
      rowcolIndex = (unsigned int)(missedSampleOriginalPos.data.at(i).get(0))  * (columns+1)
        + (unsigned int)(missedSampleOriginalPos.data.at(i).get(1));
      for( unsigned int j = 0; j < 2; j++ )
        ((this->imageCalibPoints)->data).at( rowcolIndex ).set(j, missedSampleGuesses2d.data.at(i).get(j) );
    }
    //Destroy regioncodeprobs finally
    regionCodeProps.clear();

    // #### Debug output:
    if( debug == 2 || debug == 3 ) {
      //Final point in red
      writer.getStyle().penColor = 4;
      writer.getStyle().depth= 66;
      for( unsigned int row = 0; row < rows+1 ; row++ ) {
        for( unsigned int col = 0; col < columns+1 ; col++ ) {
              // If bad vector, continue
          if( ((this->imageCalibPoints)->data).at( row*(columns+1) + col  ).get(0) < 0.0 || ((this->imageCalibPoints)->data).at( row*(columns+1) + col  ).get(0) >= arrayDim.vec[0] ||
	      ((this->imageCalibPoints)->data).at( row*(columns+1) + col  ).get(1) < 0.0 || ((this->imageCalibPoints)->data).at( row*(columns+1) + col  ).get(1) >= arrayDim.vec[1])
            continue;

          debugPoint.set(0, ((this->imageCalibPoints)->data).at( row*(columns+1) + col ).get(0) );
          debugPoint.set(1, ((this->imageCalibPoints)->data).at( row*(columns+1) + col ).get(1) );
          debugPoint += drawOffset;
          debugPoint*= debugScaleFactor;
          writer.circle( debugPoint, circleSize );
        }
      }
    }

//cerr << "debug:18" << endl;
    if( debug == 2 ) {
      //Missed new point in purple
      writer.getStyle().penColor = 5;
      writer.getStyle().depth= 65;
      for( unsigned int i = 0; i < missedSampleGuesses2d.data.size() ; i++ ) {
              // If diverged, continue
        if( ! (missedSaddleDiverged->at(i)) ) {
          debugPoint.assign( missedSampleGuesses2d.data.at(i) );
          debugPoint += drawOffset;
          debugPoint*= debugScaleFactor;
          writer.circle( debugPoint, circleSize );
        }
      }
    }

    if( debug == 2 || debug == 3 ) {
      //Now draw axes
      SampleVector xAxis;
      tempVec.copy( oneVec );
      tempVec[0] = 0.0; tempVec[1] = 0.0; xAxis.data.push_back( tempVec );
      tempVec.copy( oneVec );
      tempVec[0] = 5.0; tempVec[1] = 0.0; xAxis.data.push_back( tempVec );
      SampleVector imageXAxis;
      homotrans( &xAxis, &ransacH, &imageXAxis );

      SampleVector yAxis;
      tempVec.copy( oneVec );
      tempVec[0] = 0.0; tempVec[1] = 0.0; yAxis.data.push_back( tempVec );
      tempVec.copy( oneVec );
      tempVec[0] = 0.0; tempVec[1] = 5.0; yAxis.data.push_back( tempVec );
      SampleVector imageYAxis;
      homotrans( &yAxis, &ransacH, &imageYAxis );

      //Finally draw axes
      writer.getStyle().depth= 64;
      writer.getStyle().thickness = 5.0;
      Vector l1(2, 0.0), l2(2, 0.0);
      writer.getStyle().penColor = 4;
      l1.assign( imageXAxis.data.at(0).getSubVector(0,1) ); l1 += drawOffset; l1*= debugScaleFactor;
      l2.assign( imageXAxis.data.at(1).getSubVector(0,1) ); l2 += drawOffset; l2*= debugScaleFactor;
      writer.lineSegment( l1, l2 );
      writer.getStyle().penColor = 9;
      l1.assign( imageYAxis.data.at(0).getSubVector(0,1) ); l1 += drawOffset; l1*= debugScaleFactor;
      l2.assign( imageYAxis.data.at(1).getSubVector(0,1) ); l2 += drawOffset; l2*= debugScaleFactor;
      writer.lineSegment( l1, l2 );

      // image in the background
      writer.getStyle().depth= 70;
      writer.image( greyscaleImageFilename.c_str(), arrayDim.vec[0]*debugScaleFactor, arrayDim.vec[1]*debugScaleFactor );

      // #### End of Debug output
      writer.close();
    }

    //if( successfulTracked != NULL )
    successfulTracked = true;

    if( debug == 1 || debug == 2 )
      cerr << "Successfully tracked !" << endl << endl;

    delete missedSaddleDiverged;

    return successfulTracked;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

}							/* namespace */
#endif						/* CAMERACALIBRATION_CALTAGPATTERN_C */
