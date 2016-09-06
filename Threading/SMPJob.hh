// ==========================================================================
// $Id:$
// a single job in a symmetric multiprocessign algorithm (pure virtual class)
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

#ifndef THREADING_SMPJOB_H
#define THREADING_SMPJOB_H

/*! \file  SMPJob.hh
    \brief a single job in a symmetric multiprocessing algorithm (pure
    virtual class)
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>
#include <iostream>

#include <MDA/Config.hh>

namespace MDA {
  
  using namespace std;
  
  /** \class SMPJob SMPJob.hh
      a single job in a symmetric multiprocessing algorithm
      (virtual class that shoud always be subclassed) */
  class SMPJob {

  public:

    /** default constructor
     * \param timeEst: how expensive we expect the job to be in FLOPS
     * This helps the scheduler to decide whether or not to collate
     * multiple jobs.
     * Default is just high enough to not be collated with other jobs
     */
    SMPJob( double timeEst= SMP_COLLATION_THRESHOLD )
      : applyReduction( false ), timeEstimate( timeEst )
    {}
    
    /** destructor */
    virtual ~SMPJob() {}
    
    /** execute job (pure virtual)
     * \param threadID is an int that identifies individual threads
     * primarily for debugging
     */
    virtual void execute( int threadID )= 0;
    
    /** Reduction operator. Merges the result of this job computation
	with other results. Needs to be associative and
	commutative. The scheduler insures that only one reduction
	operator is running at any given time (very much unlike the
	main method) */
    virtual void reduce( int jobID )
    {}
    
    /** whether or not the reduction is applied */
    bool applyReduction;

    /** how time consuming the job is expected to be (in FLOPS) */
    double timeEstimate;
  };
  
  
  /** \typedef a job list */
  typedef list<SMPJob *> SMPJobList;

} /* namespace */

#endif /* THREADING_SMPJOB_H */

