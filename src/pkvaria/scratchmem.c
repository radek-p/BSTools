
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pkvprivate.h"


boolean pkv_critical = false;
boolean pkv_signal = false;
void (*pkv_signal_handler)( void ) = NULL;
void (*pkv_register_memblock)( void *ptr, boolean alloc ) = NULL;
void (*pkv_malloc_wrapper)( char _case ) = NULL;

/* pointers to procedures, which by default point to simple  */
/* allocation/deallocation routines in a scratch memory pool */
/* which is a simple stack; in an application using threads  */
/* running in parallel, alternative routines (using          */
/* a separate stack for each thread) have to be assigned     */
void *(*pkv_GetScratchMem) ( size_t size ) = NULL;
void (*pkv_FreeScratchMem) ( size_t size ) = NULL;
void *(*pkv_GetScratchMemTop) ( void ) = NULL;
void (*pkv_SetScratchMemTop) ( void *p ) = NULL;
size_t (*pkv_ScratchMemAvail) ( void ) = NULL;
size_t (*pkv_MaxScratchTaken) ( void ) = NULL;

/* private variables */
char *ScratchPtr = NULL;
char *FreeScratch = NULL;
size_t ScratchSize = 0;
size_t FreeScratchSize = 0;
size_t MinFreeScratch;

/* ////////////////////////////////////////////////////////////////////////// */
static void *_pkv_GetScratchMem ( size_t size )
{
  char *p;
  
  if ( size <= FreeScratchSize )
  {
    p = FreeScratch;
    FreeScratch += size;
    FreeScratchSize -= size;
    if ( FreeScratchSize < MinFreeScratch )
      MinFreeScratch = FreeScratchSize;
  }
  else if ( !ScratchPtr ) {
    PKV_SIGNALERROR ( LIB_PKVARIA, ERRCODE_0, ERRMSG_0 );
    return NULL;
  }
  else
    p = NULL;
  return p;
} /*_pkv_GetScratchMem*/

static void _pkv_FreeScratchMem ( size_t size )
{
  FreeScratch -= size;
  FreeScratchSize += size;
} /*_pkv_FreeScratchMem*/

static void *_pkv_GetScratchMemTop ( void )
{
  return FreeScratch;
} /*_pkv_GetScratchMemTop*/

static void _pkv_SetScratchMemTop ( void *p )
{
  int fsm;

  fsm = (char*)p - FreeScratch;
  FreeScratch = (char*)p;
  FreeScratchSize -= fsm;
} /*_pkv_SetScratchMemTop*/

static size_t _pkv_ScratchMemAvail ( void )
{
  return FreeScratchSize;
} /*_pkv_ScratchMemAvail*/

static size_t _pkv_MaxScratchTaken ( void )
{
  size_t s;

  s = ScratchSize-MinFreeScratch;
  MinFreeScratch = FreeScratchSize;
  return s;
} /*_pkv_MaxScratchTaken*/

/* ////////////////////////////////////////////////////////////////////////// */
void _pkv_AssignDefaultScratchMemProc ( void )
{
  pkv_GetScratchMem    = _pkv_GetScratchMem;
  pkv_FreeScratchMem   = _pkv_FreeScratchMem;
  pkv_GetScratchMemTop = _pkv_GetScratchMemTop;
  pkv_SetScratchMemTop = _pkv_SetScratchMemTop;
  pkv_ScratchMemAvail  = _pkv_ScratchMemAvail;
  pkv_MaxScratchTaken  = _pkv_MaxScratchTaken;
} /*_pkv_AssignDefaultScratchMemProc*/

boolean pkv_InitScratchMem ( size_t size )
{
  _pkv_AssignDefaultScratchMemProc ();
  PKV_MALLOC ( ScratchPtr, size );
  if ( !ScratchPtr ) {
    PKV_SIGNALERROR ( LIB_PKVARIA, ERRCODE_1, ERRMSG_1 );
    return false;
  }
  FreeScratch = ScratchPtr;
  MinFreeScratch = FreeScratchSize = ScratchSize = size;
  return true;
} /*pkv_InitScratchMem*/

void pkv_DestroyScratchMem ( void )
{
  if ( ScratchPtr ) {
    PKV_FREE ( ScratchPtr );
    FreeScratchSize = ScratchSize = 0;
  }
} /*pkv_DestroyScratchMem*/

void PrintScratchMemData ( void )
{
#if __WORDSIZE == 64
  printf ( "ScratchPtr = 0x%lx, FreeScratch = 0x%lx, ScratchSize = %ld, FreeScratchSize = %ld\n",
   (size_t)ScratchPtr, (size_t)FreeScratch, ScratchSize, FreeScratchSize );
#else
  printf ( "ScratchPtr = 0x%x, FreeScratch = 0x%x, ScratchSize = %d, FreeScratchSize = %d\n",
   (size_t)ScratchPtr, (size_t)FreeScratch, ScratchSize, FreeScratchSize );
#endif
} /*PrintScratchMemData*/

