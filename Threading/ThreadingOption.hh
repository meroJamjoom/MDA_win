// ==========================================================================
// $Id:$
// command line options for configuring multithreading
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

#ifndef THREADING_THREADINGOPTION_H
#define THREADING_THREADINGOPTION_H

/*! \file  ThreadingOption.hh
    \brief command line options for configuring multithreading
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include "MDA/Base/CommandlineParser.hh"
#include "SMPJobManager.hh"


namespace MDA {

  /** \class ThreadingOption ThreadingOption.hh
      command line options for configuring multithreading */
  
  class ThreadingOption: public ScalarOption<int> {

  public:

    /** default constructor */
    ThreadingOption()
      : ScalarOption<int>(static_cast<int>(SMPJobManager::numThreads),
			   "\tnumber of threads to use\n",
			   "--threads" , NULL )
    {}
    
  protected:

  };


} /* namespace */

#endif /* THREADING_THREADINGOPTION_H */

