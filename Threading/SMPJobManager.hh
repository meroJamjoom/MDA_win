// ==========================================================================
// $Id:$
// job manager for symmetric multiprocessing
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

#ifndef THREADING_SMPJOBMANAGER_H
#define THREADING_SMPJOBMANAGER_H

/*! \file  SMPJobManager.hh
    \brief job manager for symmetric multiprocessing
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <list>
#ifdef HAVE_PTHREADS
#include <pthread.h>
#endif

#include "SMPJob.hh"


namespace MDA {
  
  using namespace std;
  
  /** \class SMPJobManager SMPJobManager.hh
      job manager for symmetric multiprocessing */
  
  class SMPJobManager {
    
  private:

    /** constructor is private to ensure there is ever only one
	instance (i.e. the SMBJobManager is a singleton) */
    SMPJobManager();
    
  public:
    
    /** return pointer to the job manager object */
    inline static SMPJobManager *getJobManager()
    {
      if( jobManager== NULL )
	jobManager= new SMPJobManager;
      return jobManager;
    }
    
    /** batch a list of jobs for execution & block until termination
     *  Use serial execution if job list not previously empty The job
     *  manager removes all jobs from the JobList, and returns it
     *  empty.
     */
    void batch( SMPJobList &jobs );
    
    /** set the threshold (minimum process cost) for combining
	multiple jobs */
    inline void setCollationThreshold( double threshold )
    {
      collationThreshold= threshold;
    }
    
    /** returns the number of threads allocated by the manager */
    inline static int getNumThreads()
    {
      return numThreads;
    }

  protected:

    /** wait until job list is empty */
    void wait();
    
    /** serially execute a list of jobs */
    void serialExecute( SMPJobList &jobs );
    
    /** collation threshold for small jobs */
    static double collationThreshold;
    
    /** the global job manager object */
    static SMPJobManager *jobManager;
    
    /** current job list */
    static SMPJobList jobList;
    
    /** how many threads we should use on this system */
    static const int numThreads;
    
    /** ThreadingOption can set some of the configration values */
    friend class ThreadingOption;
    
#ifdef HAVE_PTHREADS
    /** number of busy threads (0 if all jobs complete) */
    static int numBusy;
    
    /** a thread that waits for/gets jobs and executes them */
    static void *runThread( void *threadData );
    
    /** mutex for job list */
    static pthread_mutex_t *jobMutex;
    
    /** mutex for reduction */
    static pthread_mutex_t *reduceMutex;
    
    /** condition varialble to wait for jobs in the job list */
    static pthread_cond_t *jobListReady;
    
    /** condition varialble to signal that all batch jobs are done */
    static pthread_cond_t *allThreadsIdle;
#endif
  };


} /* namespace */

#endif /* THREADING_SMPJOBMANAGER_H */

