
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef PKVTHREADS_H
#define PKVTHREADS_H

#ifndef _PTHREAD_H
#include <pthread.h>
#endif

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef void*(*PKVThreadProc)(void*);
typedef boolean (*PKVThreadWorkToDo)(void*,int3*);

typedef struct {
    short int         next;        /* list link */
    pthread_t         thread;      /* pthread identifier */
    pthread_mutex_t   mutex;
    pthread_cond_t    cond;
    PKVThreadWorkToDo jobproc;
    void              *jobdata;
    void              *auxdata;
    int3              jobnum;
    boolean           valid, waiting, success, pad;
    char              *ScratchPtr, *FreeScratchPtr;
    size_t            ScratchSize, FreeScratchSize, MinFreeScratch;
  } pkv_thread;


extern pthread_mutex_t thread_mutex;
extern pthread_t       main_thread;

extern short int       max_threads;
extern pkv_thread      *pkvthread;

boolean pkv_InitPThreads ( short int maxthreads );
void pkv_DestroyPThreads ( void );

short int pkv_PThreadMyPos ( pthread_t *thr );
short int pkv_PThreadIPos ( pthread_t thr );

void pkv_CancelPThread ( short int pos );
void pkv_CancelPThreads ( void );

boolean pkv_SetPThreadsToWork ( int3 *jobsize, int npthreads,
                                size_t stacksize, size_t scratchmemsize,
                                void *usrdata,
                                PKVThreadWorkToDo jobproc,
                                void *extradata,
                                PKVThreadWorkToDo extrajob,
                                boolean *success );

short int pkv_NewJoinablePThread ( size_t stacksize, size_t scratchmemsize,
                                   PKVThreadWorkToDo jobproc, void *jobdata,
                                   int3 *jobnum, void *auxdata,
                                   pthread_t *thread );

int pkv_FindNCPU ( void );
#ifdef __cplusplus
}
#endif

#endif /*PKVTHREADS_H*/

