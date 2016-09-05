// ==========================================================================
// $Id: MemoryObject.hh 547 2010-01-06 01:42:18Z skoch $
// reference counter object for large memory regions
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

#ifndef LINEARALGEBRA_MEMORYOBJECT_H
#define LINEARALGEBRA_MEMORYOBJECT_H

/*! \file  MemoryObject.hh
    \brief reference counter object for large memory regions
 */

#ifdef _WIN32
// this header file must be included before all the others
#define NOMINMAX
#include <windows.h>
#endif

#include <stdlib.h>

#include "MDA/base/Errors.hh"

namespace MDA {
  
  using namespace std;
  
  /** \class MemoryObject MemoryObject.hh
      reference counter object for large memory regions */
  
  class MemoryObject {
    
    
  public:
    
    /** constructor */
    inline MemoryObject( double *_data= NULL, MemoryObject *_owner= NULL )
      : owner( NULL ), data( _data )
    {
      // if we are sharing data with another linear algebra object,
      // make sure that data stays around for our lifetime
      refCounter = 0;
      if( _owner!= NULL )
	owner= _owner->ref();
      // create initial reference to self
      ref();
    }
    
    /** new reference to the object */
    inline MemoryObject *ref()
    {
      refCounter++;
      return this;
    }
    
    /** dereference the object */
    inline void deref()
    {
      if( --refCounter== 0 )
	delete this;
    }
    
    /** pointer to raw memory */
    double *mem;
    
    /** the data vector */
    double *data;
    
    /** the "owner" of the data (equals this if we own it, NULL if no
	linear algebra object does) */
    MemoryObject *owner;
    
  protected:
    
    /** destructor */
    virtual ~MemoryObject()
    {
      if( owner== this )
      {
	if( mem!= NULL )
	  delete [] mem;
      }
      else
	if( owner!= NULL )
	  owner->deref();
    }
    
    /** create new memory for "size" double values. The data pointer
	is guaranteed to be aligned with 16 Byte boundaries. Also,
	data[-1] is a valid address, and can be used to store data
	such as the background color of a channel. */
    inline void allocate( unsigned long size )
    {
      errorCond( data== NULL, "  cannot change memory allocation!\n" );
      
      // create new memory
      if( size!= 0ul )
      {
	// Allocate 2 extra doubles to ensure alignent of the the data
	// pointer to 16 byte boundaries.
	mem= new double[size+2];
	data= (double *)(((size_t)mem+size_t(16)) & ~size_t(15));
	owner= this;
      }
    }
    
    /** reference counter */
    unsigned refCounter;
    
  };


} /* namespace */

#endif /* LINEARALGEBRA_MEMORYOBJECT_H */

