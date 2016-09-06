// ==========================================================================
// $Id:$
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

#ifndef THREADING_SMPJOBMANAGER_C
#define THREADING_SMPJOBMANAGER_C

#include <iostream>

#include <MDA/Config.hh>
#include "SMPJobManager.hh"

namespace MDA {

  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  //
  // initializations od static members
  //

// the job list
SMPJobList SMPJobManager::jobList;
// no job manager object initially
SMPJobManager *SMPJobManager::jobManager= NULL;
// default threshold for merging jobs (in FLOPS)
double SMPJobManager::collationThreshold= SMP_COLLATION_THRESHOLD;
// multithreading off by default
const int SMPJobManager::numThreads= 1;


#ifdef HAVE_PTHREADS
pthread_mutex_t *SMPJobManager::jobMutex= NULL;
pthread_mutex_t *SMPJobManager::reduceMutex= NULL;
pthread_cond_t *SMPJobManager::jobListReady= NULL;
pthread_cond_t *SMPJobManager::allThreadsIdle= NULL;
int SMPJobManager::numBusy= 0;
#endif


/** constructor is private to ensure there is ever only one
    instance (i.e. the SMBJobManager is a singleton) */
SMPJobManager::SMPJobManager()
{
  jobList.clear();
  int returnVal;
  
  if( numThreads> 1 )
  {
#ifdef HAVE_PTHREADS
    // initially, set all threads to busy
    numBusy= numThreads;

    // create mutex, condition variable
    jobMutex= new pthread_mutex_t;
    pthread_mutex_init( jobMutex, NULL );
    reduceMutex= new pthread_mutex_t;
    pthread_mutex_init( reduceMutex, NULL );
    jobListReady= new pthread_cond_t;
    pthread_cond_init( jobListReady, NULL );
    allThreadsIdle= new pthread_cond_t;
    pthread_cond_init( allThreadsIdle, NULL );
    
    // create threads and make them joinable
    pthread_t threads[numThreads];
    pthread_attr_t joinable;
    pthread_attr_init( &joinable );
    pthread_attr_setdetachstate( &joinable, PTHREAD_CREATE_JOINABLE );
    for( long i= 0 ; i< numThreads ; i++ )
    {
      returnVal= pthread_create( &threads[i], &joinable, runThread, (void *)i );
      if( returnVal )
      {
	cerr << "Error creating thread " << i << endl;
	numBusy--; // threads that didn't spawn are not busy...
      }
    }
#endif
  }
}


#ifdef HAVE_PTHREADS
/** a thread that waits for/gets jobs and executes them */
void *
SMPJobManager::runThread( void *threadData )
{
  double jobCost;
  SMPJobList myJobs;
  
  while( true )
  {
    // get next job from job list (quit if job list empty)
    pthread_mutex_lock( jobMutex );
    while( jobList.empty() )
    {
#ifdef DEBUG_THREADS
      cerr << "Job " << (int)(threadData) << " waiting for jobs\n";
#endif
      if( --numBusy<= 0 )
      {
#ifdef DEBUG_THREADS
	cerr << "Signaling idle...\n";
#endif
	pthread_cond_signal( allThreadsIdle );
      }
      pthread_cond_wait( jobListReady, jobMutex );
      numBusy++;
    }
    // grab as many small jobs as needed to bring the total
    // computational cost above collationThreshold
    for( jobCost= 0.0 ; jobCost< collationThreshold && !jobList.empty() ; )
    {
      SMPJob *h= jobList.front();
      myJobs.push_back( h );
      jobCost+= h->timeEstimate;
      jobList.pop_front();
    }
    pthread_mutex_unlock( jobMutex );
    
    // execute the jobs
    while( !myJobs.empty() )
    {
      SMPJob *nextJob= myJobs.front();
      nextJob->execute( (long int)threadData );
      
      // use reduction operator if reduction isn't the null operator
      if( nextJob->applyReduction )
      {
	pthread_mutex_lock( reduceMutex );
	nextJob->reduce( (long int)threadData );
	pthread_mutex_unlock( reduceMutex );
      }
      
      // delete job after we are done
      delete nextJob;
      myJobs.pop_front();
    }
  }
  
  // never reached
  pthread_mutex_unlock( jobMutex );
  pthread_exit( NULL );
  return NULL;
}
#endif


/** serially execute a list of jobs */
void
SMPJobManager::serialExecute( SMPJobList &jobs )
{
  while( !jobs.empty() )
  {
    SMPJob *current= jobs.front();
    jobs.pop_front();
    current->execute( -1 );
    if( current->applyReduction )
      current->reduce( -1 );
    delete current;
  }
}


/** batch a list of jobs for execution & block until termination
 *  (use serial execution if job list not previously empty)
 */
void
SMPJobManager::batch( SMPJobList &jobs )
{
#ifdef HAVE_PTHREADS
  if( numThreads> 1)
  {
    pthread_mutex_lock( jobMutex );
    if( jobList.empty() )
    {
      // we have pthreads, and the threading is not already active
      // - thus, actually create the threads for parallel processing
      
      // remove all jobs from the provided list, and add them to jobList
      jobList.splice( jobList.begin(), jobs );
      pthread_cond_broadcast( jobListReady );
      pthread_mutex_unlock( jobMutex );
      
      // now wait for completion
      wait();
    }
    else
    {
      // we have pthreads, but the threads are already active working
      // on a job list -> serial execution of the new jobs
      pthread_mutex_unlock( jobMutex );
      serialExecute( jobs );
    }
  }
  else
#endif
  {
    // if threading is not compiled or enabled, use serial execution
    serialExecute( jobs );
  }
}

/** wait until job list is empty */
void
SMPJobManager::wait()
{
#ifdef HAVE_PTHREADS
  if( numThreads> 1 )
  {
    pthread_mutex_lock( jobMutex );
    pthread_cond_wait( allThreadsIdle, jobMutex );
    pthread_mutex_unlock( jobMutex );
  }
#endif
}
    

} /* namespace */

#endif /* THREADING_SMPJOBMANAGER_C */

