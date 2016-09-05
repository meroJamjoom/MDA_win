// ==========================================================================
// $Id: Errors.hh 199 2008-04-08 03:58:26Z heidrich $
// methods and macros for dealing with error messages
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2006-, UBC
// 
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef BASE_ERRORS_H
#define BASE_ERRORS_H

#include <stdlib.h>
#include <iostream>

/*! \file  Errors.hh
    \brief methods and macros for dealing with error messages and warnings
 */

#ifdef DEBUG

/* unconditional warning */
#define warning( msg ) __warning( msg, __FILE__, __LINE__, __FUNCTION__ )

/* if an error occurs in debug mode, warnCond outputs the error
   message followed by the location of the error. Result value is the
   value of the condition (false if error, true if everything is ok) */
#define warnCond( cond, msg ) __warnCond( cond, msg, __FILE__, __LINE__,\
					  __FUNCTION__ )

/* warning about allocating a large chunk of memory */
#define warnMem( size, thres ) __warnMem( size, thres, __FILE__, __LINE__,\
					  __FUNCTION__ )

/* if an error occurs in debug mode, errorCond outputs the error
   message followed by the location of the error. The program
   temrinates if there is an error */
#define errorCond( cond, msg ) __errorCond( cond, msg, #cond,\
					    __FILE__, __LINE__, __FUNCTION__ )

namespace MDA {

/** internal method for unconditional warning */
inline void
__warning( const char *msg, const char *file, unsigned line, const char* func )
{
  std::cerr << "WARNING: " << file << ", line " << line
	    << " (" << func << "):" << std::endl
	    << msg << std::endl;
}

/** internal method for issuing warnings (if condition is false, the
    warning is triggered) */
inline bool
__warnCond( bool cond, const char *msg,
	    const char *file, unsigned line, const char* func )
{
  if( cond )
    return true;
  std::cerr << "WARNING: " << file << ", line " << line
	    << " (" << func << "):" << std::endl
	    << msg << std::endl;
  return false;
}

/** internal method for issuing memory allocation warnings
    (triggered if mem is above thres) */
inline void
__warnMem( unsigned long mem, unsigned long thres,
	   const char *file, unsigned line, const char* func )
{
  static char prefixes[]= { 'k', 'M', 'G', 'T', 'P' };
  if( mem>= thres )
  {
    int i= 0;
    do
      mem/= 1024ul;
    while( i< 5 && mem> 1024 );
    std::cerr << "WARNING: " << file << ", line" << line
	      << " (" << func << "):" << std::endl
	      << "Allocating more than " << mem << ' '
	      << prefixes[i] << "B of memory!" << std::endl;
  }
}
    
/** internal method for issuing errors (if condition is false, the
    error is triggered, and the program aborted) */
inline void
__errorCond( bool cond, const char *msg, const char *condStr,
	     const char *file, unsigned line, const char *func )
{
  if( cond )
    return;
  std::cerr << "ERROR: " << file << ", line " << line
	    << " (" << func << "):" << std::endl
	    << msg << std::endl
	    << "  " << condStr << std::endl;
  abort();
}

  

} /* namespace */


#else

/* in optimized mode, no warnings or errors */
#define warning( msg ) (true)
#define warnCond( cond, msg ) (cond)
#define warnMem( mem, thres ) (true)
#define errorCond( cond, msg ) ((void) 0)

#endif


#endif /* BASE_ERRORS_H */

