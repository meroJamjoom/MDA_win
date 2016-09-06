// ==========================================================================
// $Id: PoissonSystem.hh 849 2011-08-29 22:09:48Z joelaf $ 
// An image-space linear system describing a Poisson equation
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

#ifndef IMAGESPACESYSTEMS_POISSONSYSTEM_H
#define IMAGESPACESYSTEMS_POISSONSYSTEM_H

/*! \file  PoissonSystem.hh
  \brief An image-space linear system describing a Poisson equation
*/

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Array/Array.hh"
#include "ImageSpaceSystem.hh"

namespace MDA {

  /** \class PoissonSystem PoissonSystem.hh
      An image-space linear system describing a Poisson equation */
  template<class T>
  class PoissonSystem: public MultigridSystem {

  protected:

    /** reference to gradients */
    T ** gradPtr;

    /** reference to mask channel */
    T * maskPtr;

    /** reference to constraints channel */
    T *constraints;

    /** smaller array required for multigrid computation */
    Array<T> a2;

    /** main setup method; everything implemented here */
    bool _setup(Array<T> &a, 
		Vector & rhs, 
		ChannelList &gradientCh, 
		unsigned constraintCh, 
		unsigned maskCh, 
		bool computeDivergence, 
		unsigned sourceCh, 
		bool contract,
		bool multigrid=false);
 
    /** convenience method for computing the div(grad field) at a point*/
    inline T computeRHS(long index, bool computeDivergence, T * potential);

    /** rescaler function. returns: unconstrained count in the half matrix */
    inline int halfRescale(Array<T> &a, 
			   Array<T> &a2, 
			   int maskCh, 
			   int constraintCh,
			   ChannelList &gradientCh);

    /** convenience method for computing the contribution to the RHS vector
	from the boundary conditions*/
    inline T computeBCContribution(long index, T * constraints);

    /** function to create the half sized system; same as setup()
	except included computing "half" array  */  
    inline PoissonSystem<T>* createHalfSystem( Array<T> &a, 
					       //Vector & rhs, 
					       ChannelList &gradientCh, 
					       unsigned constraintCh, 
					       unsigned maskCh, 
					       bool computeDivergence,
					       unsigned sourceCh, 
					       bool cleanMatrix);

  public:

    enum { 
      valueMask = 0, 
      gradientMask = -1, 
      unconstrainedMask = 1 
    };

    /** default constructor */
    inline PoissonSystem(int gridLevel=0):gradPtr(NULL), 
					  maskPtr(NULL),
					  constraints(NULL),
					  MultigridSystem(gridLevel)
    {}
    
    /** destructor */
    inline ~PoissonSystem(){}
    
    /** setup - includes computing the RHS vector*/
    inline bool setup( Array<T> &a, 
		       Vector & rhs, 
		       ChannelList &gradientCh, 
		       unsigned constraintCh, 
		       unsigned maskCh, 
		       bool computeDivergence, 
		       unsigned sourceCh,
		       bool multigrid )
    {
      return _setup(a,rhs,gradientCh,constraintCh,maskCh,
		    computeDivergence,sourceCh, false, multigrid);
    }

    /** setup - computes divergence by default */
    inline bool setup( Array<T> &a, 
		       Vector &rhs, 
		       ChannelList &gradientCh, 
		       unsigned constraintCh, 
		       unsigned maskCh ) {
      return _setup(a,rhs,gradientCh,constraintCh,maskCh,true,0, false);
    } 

    /** setup - removes all identity entries from the matrix */
    inline bool cleanSetup(Array<T> &a, 
			   Vector & rhs, 
			   ChannelList &gradientCh, 
			   unsigned constraintCh, 
			   unsigned maskCh, 
			   bool computeDivergence, 
			   unsigned sourceCh ) {
      return _setup(a,rhs,gradientCh,constraintCh,
		    maskCh,computeDivergence,sourceCh, true);
    }

    /** setup - for creating an image space system without 
	any gradient information */
    bool setup( Array<T> &a, unsigned maskCh);

    /** recalculates RHS with new "source" data */
    void recalculateRHS(T * source, 
			T * constraints, 
			Vector &rhs, 
			bool computeDivergence=false);
		   
    /** rightmultiply is called by the solver */
    Vector &rightMultiply( const Vector &v, Vector &result ) const;
    
    /** multigrid prolongation operation */
    virtual inline void prolongate(const Vector &half, Vector &v);
    
    /** multigrid restriction operation */
    virtual inline void restrict(const Vector &v, Vector &half);
  };


} /* namespace */

#endif /* IMAGESPACESYSTEMS_POISSONSYSTEM_H */

