// ==========================================================================
// $Id: ImageSpaceSystem.hh 999 2014-05-28 15:07:31Z heidrich $
// baseclass for image-space linear systems
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

#ifndef IMAGESPACESYSTEMS_IMAGESPACESYSTEM_H
#define IMAGESPACESYSTEMS_IMAGESPACESYSTEM_H

/*! \file  ImageSpaceSystem.hh
    \brief baseclass for image-space linear systems
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <utility>
#include <vector>

#include "MDA/LinearAlgebra/LinAlg.hh"

namespace MDA {
 

  using namespace std;
  
  /**convenience method for setting a bit at a particular position (starting
     from the right) to a particular value*/
  inline void setBit(unsigned & flag, unsigned short index, bool value)
  {
    flag&= ~(1 << index); 
    flag|= (value ? 1 : 0) << index;
  }
  
  /**convenience method for obtaining a bit at a particular position (starting
     from the right)*/
  inline bool getBit(unsigned flag, unsigned index)
  {
    return (flag >> index & 1);
  } 
  
  /** \class ImageSpaceSystem ImageSpaceSystem.hh
      baseclass for image-space linear systems */
  
  class ImageSpaceSystem: public LinearOperator {
    
  protected:
    
    /** dimension of the space */
    unsigned dimension;
    
    /** array dimensions */
    CoordinateVector dim;

    /** total number of points in system = product_i(dim[i]) */
    unsigned long totalPts;
    
    /** determines offsets for moving a unit in a particular dimension*/ 
    long * offsetTable;
    
    /** neighbouring pixel configuration (of the original image,
	before constraints are applied) */
    unsigned * neighbourConfig;
    
    /** count of support size for given bit vector */
    unsigned char *supportSize;
    
    /** number of different configurations */
    unsigned numCases;
    
    /** offsets and weights for all possible non-zero element configurations */
    pair<long,double> **elems;
    
    /** array bitvectors describing, for each pixel in a given image
	system, which of the possible neighbors are active
     */
    unsigned *nonZeroes;

    /** setup the offset table */
    virtual void setupOffsetTable();

    /** setup the tables describing different pixel configuration
	cases */
    virtual void setupCaseTable();

    /** determine the neighbourhood configuration for each pixel */
    void computeNeighborCounts ( );

  public:
    
    /** default constructor */
    inline ImageSpaceSystem()
      : nonZeroes( NULL ), supportSize( NULL ), elems( NULL ), numCases( 0 ), 
	offsetTable(NULL), neighbourConfig(NULL)	
    {}
    
    /** destructor */
    virtual ~ImageSpaceSystem();
    
    /** multiply a vector with the image-space system (from the right) */
    virtual Vector &rightMultiply( const Vector &v, Vector &result ) const;
    
    /** left multiply for symmetric systems */
    inline virtual Vector &leftMultiply( const Vector &v, Vector &result ) const { return rightMultiply( v, result ); }
    
    inline virtual unsigned char getSupportSize( unsigned long row ) const { return supportSize[ nonZeroes[row] ]; }

    inline virtual unsigned getNonzeroConfig( unsigned long row ) const { return nonZeroes[row]; }

    inline virtual void getOffsetConfig( unsigned long row, pair<long,double>*& es ) const 
    {
      unsigned char i;
      unsigned nZ = nonZeroes[row];
      unsigned char size = supportSize[ nZ ];
      pair<long,double>* temp = elems[nZ];
      for(i=0; i<=size; i++)
      {
	es[i].first = temp[i].first;
	es[i].second = temp[i].second;
      }
    }    

    virtual CoordinateVector getDimensions() const { return dim; }
    
    /** returns number of rows in the image */
    virtual inline unsigned long getR(){ return dim.vec[1]; }
    
    /** returns number of columns in the image */
    virtual inline unsigned long getC(){ return dim.vec[0]; }

    /** returns number of rows of the poisson system. equals to
	row+col of the image */
    virtual inline unsigned long getNumRows() const { return totalPts; }
    
    /** returns number of columns of the poisson system. equals to
	row+col of the image */
    virtual inline unsigned long getNumColumns() const { return totalPts; }

  };




  /** \class MultigridSystem MultigridSystem.hh
      baseclass for multigrid linear systems */

  class MultigridSystem: public ImageSpaceSystem {
    
  protected:
    
    /** grid level */
    int gridLevel;

    /** points to the  half system needed for multigrid. */ 
    MultigridSystem* halfSystem;
   
  public:

    CoordinateVector getHalfDimension(){
      CoordinateVector v=dim;
      for (int i=0; i<v.vec.size(); i++)
	v.vec[i]=(v.vec[i]+1)/2; //ceil
      return v;
    }
    
    MultigridSystem(int _gridLevel):gridLevel(_gridLevel),halfSystem(NULL),ImageSpaceSystem(){}

    ~MultigridSystem(){
      if (halfSystem!=NULL) delete halfSystem;
    }
    
    /** multigrid prolongation operation */
    virtual inline void prolongate(const Vector &half, Vector &v)=0;
    
    /** multigrid restriction operation */
    virtual inline void restrict(const Vector &v, Vector &half)=0;
    
    /** returns grid level of this system */
    virtual int getGridLevel(){ return gridLevel; }
    
    /** returns the half system */
    virtual MultigridSystem* getHalfSystem(){ return halfSystem; }
    
  };

} /* namespace */

#endif /* IMAGESPACESYSTEMS_IMAGESPACESYSTEM_H */

