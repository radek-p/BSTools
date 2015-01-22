
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include "pkvaria.h"
#include "pkvthreads.h"

#include "pkvprivate.h"

/*#define DEBUG*/

/* ////////////////////////////////////////////////////////////////////////// */
/* I do not know about processors, whose stacks grow up (to higher addresses) */
/* but if such a thing exists, then before compilation for it, delete or */
/* comment out the line below */
#define CPU_STACK_GROWS_DOWN

typedef struct {
    short int pos;
    size_t    stack_var;
  } thrstack_desc;

/* ////////////////////////////////////////////////////////////////////////// */
boolean                pkv_threads_in_use = false;
pthread_mutex_t        thread_mutex;
pthread_t              main_thread;

short int              max_threads = 0;
pkv_thread             *pkvthread = NULL;

static short int       pkvt_busylist = -1,
                       pkvt_freelist = -1,
                       pkvt_firstfree = 0;

static short int       current_thrstacks = 0;
static pthread_mutex_t cpustack_mutex;
static thrstack_desc   *thrstacks = NULL;

/* ////////////////////////////////////////////////////////////////////////// */
/* as the thread identifiers (pthread_t) are opaque objects, linear search is */
/* used in the procedures. */
short int pkv_PThreadIPos ( pthread_t thr )
{
  int pos;

  pthread_mutex_lock ( &thread_mutex );
  for ( pos = pkvt_busylist; pos >= 0; pos = pkvthread[pos].next )
    if ( pkvthread[pos].valid &&
         pthread_equal ( thr, pkvthread[pos].thread )  )
      break;
  pthread_mutex_unlock ( &thread_mutex );
  return pos;
} /*pkv_PThreadIPos*/

short int _pkv_PThreadMyPos ( pthread_t *thr )
{
  pthread_t myself;

  myself = pthread_self ();
  if ( thr )
    *thr = myself;
  return pkv_PThreadIPos ( myself );
} /*_pkv_PThreadMyPos*/

/* ////////////////////////////////////////////////////////////////////////// */
static void RegisterThreadStack ( int pos, char *stack_var )
{
  int           i;
  thrstack_desc s;

  s.pos = pos;
  s.stack_var = (size_t)stack_var;
  pthread_mutex_lock ( &cpustack_mutex );
#ifdef DEBUG
printf ( "r: %d, ", pos );
#endif
        /* insertion sort */
#ifdef CPU_STACK_GROWS_DOWN
  for ( i = current_thrstacks; i > 0; i-- ) {
    if ( s.stack_var > thrstacks[i-1].stack_var )
      thrstacks[i] = thrstacks[i-1];
    else
      break;
  }
#else
  for ( i = current_thrstacks; i > 0; i-- ) {
    if ( s.stack_var < thrstacks[i-1].stack_var )
      thrstacks[i] = thrstacks[i-1];
    else
      break;
  }
#endif
  thrstacks[i] = s;
  current_thrstacks ++;
  pthread_mutex_unlock ( &cpustack_mutex );
} /*RegisterThreadStack*/

static void WithdrawThreadStack ( short int pos )
{
  int i;

  pthread_mutex_lock ( &cpustack_mutex );
#ifdef DEBUG
printf ( "w: %d, ", pos );
#endif
  for ( i = 0; i < current_thrstacks; i++ )
    if ( thrstacks[i].pos == pos ) {
      current_thrstacks --;
      if ( i < current_thrstacks )
        memmove ( &thrstacks[i], &thrstacks[i+1],
                  (current_thrstacks-i)*sizeof(thrstack_desc) );
      break;
    }
  pthread_mutex_unlock ( &cpustack_mutex );
} /*WithdrawThreadStack*/

short int pkv_PThreadMyPos ( void )
{
  size_t    qqaddr;
  short int a, b, c, pos;

  qqaddr = (size_t)&qqaddr;
  pthread_mutex_lock ( &cpustack_mutex );
      /* binary search */
  a = 0;  b = current_thrstacks;
#ifdef CPU_STACK_GROWS_DOWN
  while ( b-a > 1 ) {
    c = (a+b)/2;
    if ( thrstacks[c].stack_var > qqaddr )
      a = c;
    else
      b = c;
  }
#else
  while ( b-a > 1 ) {
    c = (a+b)/2;
    if ( thrstacks[c].stack_var < qqaddr )
      a = c;
    else
      b = c;
  }
#endif
  pos = thrstacks[a].pos;
        /* just in case */
  pthread_mutex_unlock ( &cpustack_mutex );
#ifdef DEBUG
{
  pthread_t myself;

  myself = pthread_self ();
  if ( (pos >= 0 && !pthread_equal ( myself, pkvthread[pos].thread )) ||
       (pos < 0 && !pthread_equal ( myself, main_thread )) ) {
    pos = _pkv_PThreadMyPos ( NULL );
    printf ( "Q: %d %d, ", thrstacks[a].pos, pos );
  }
}
#endif
  return pos;
} /*pkv_PThreadMyPos*/

/* ////////////////////////////////////////////////////////////////////////// */
static void *_pkv_PTHRGetScratchMem ( size_t size )
{
  short int my_pos;
  char      *p;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 ) {
    if ( !pkvthread[my_pos].ScratchPtr ) {
      PKV_SIGNALERROR ( LIB_PKVARIA, ERRCODE_0, ERRMSG_0 );
      return NULL;
    }
    else if ( size <= pkvthread[my_pos].FreeScratchSize ) {
      p = pkvthread[my_pos].FreeScratchPtr;
      pkvthread[my_pos].FreeScratchPtr += size;
      pkvthread[my_pos].FreeScratchSize -= size;
      if ( pkvthread[my_pos].FreeScratchSize < pkvthread[my_pos].MinFreeScratch )
        pkvthread[my_pos].MinFreeScratch = pkvthread[my_pos].FreeScratchSize;
    }
    else
      p = NULL;
  }
  else {
    if ( !ScratchPtr ) {
      PKV_SIGNALERROR ( LIB_PKVARIA, ERRCODE_0, ERRMSG_0 );
      return NULL;
    }
    else if ( size <= FreeScratchSize ) {
      p = FreeScratch;
      FreeScratch += size;
      FreeScratchSize -= size;
      if ( FreeScratchSize < MinFreeScratch )
        MinFreeScratch = FreeScratchSize;
    }
    else
      p = NULL;
  }
  return p;
} /*_pkv_PTHRGetScratchMem*/

static void _pkv_PTHRFreeScratchMem ( size_t size )
{
  short int my_pos;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 ) {
    pkvthread[my_pos].FreeScratchPtr -= size;
    pkvthread[my_pos].FreeScratchSize += size;
  }
  else {
    FreeScratch -= size;
    FreeScratchSize += size;
  }
} /*_pkv_PTHRFreeScratchMem*/

static void *_pkv_PTHRGetScratchMemTop ( void )
{
  short int my_pos;
  void      *p;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 )
    p = pkvthread[my_pos].FreeScratchPtr;
  else
    p = FreeScratch;
  return p;
} /*_pkv_PTHRGetScratchMemTop*/

static void _pkv_PTHRSetScratchMemTop ( void *p )
{
  short int my_pos;
  long int  fsm;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 ) {
    fsm = (char*)p - pkvthread[my_pos].FreeScratchPtr;
    pkvthread[my_pos].FreeScratchPtr = (char*)p;
    pkvthread[my_pos].FreeScratchSize -= fsm;
  }
  else {
    fsm = (char*)p - FreeScratch;
    FreeScratch = (char*)p;
    FreeScratchSize -= fsm;
  }
} /*_pkv_PTHRSetScratchMemTop*/

static size_t _pkv_PTHRScratchMemAvail ( void )
{
  short int my_pos;
  size_t    s;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 )
    s = pkvthread[my_pos].FreeScratchSize;
  else
    s = FreeScratchSize;
  return s;
} /*_pkv_PTHRScratchMemAvail*/

static size_t _pkv_PTHRMaxScratchTaken ( void )
{
  short int my_pos;
  size_t    s;

  my_pos = pkv_PThreadMyPos ();
  if ( my_pos >= 0 ) {
    s = pkvthread[my_pos].ScratchSize - pkvthread[my_pos].MinFreeScratch;
    pkvthread[my_pos].MinFreeScratch = pkvthread[my_pos].FreeScratchSize;
  }
  else {
    s = ScratchSize - FreeScratchSize;
    MinFreeScratch = FreeScratchSize;
  }
  return s;
} /*_pkv_PTHRMaxScratchTaken*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean pkv_InitPThreads ( short int maxthreads )
{
  int rc;

  main_thread = pthread_self ();
  PKV_MALLOC ( thrstacks, (maxthreads+1)*sizeof(thrstack_desc) );
  PKV_MALLOC ( pkvthread, maxthreads*sizeof(pkv_thread) );
  if ( !thrstacks || !pkvthread )
    goto failure;
  max_threads = maxthreads;
  current_thrstacks = 0;
  RegisterThreadStack ( -1, (char*)&rc );
  pkvt_busylist = pkvt_freelist = -1;
  pkvt_firstfree = 0;

  rc = pthread_mutex_init ( &cpustack_mutex, NULL );
  if ( rc )
    goto failure;
  rc = pthread_mutex_init ( &thread_mutex, NULL );
  if ( rc ) {
    pthread_mutex_destroy ( &cpustack_mutex );
    goto failure;
  }
        /* substitite the scratch memory procedures */
  pkv_GetScratchMem    = _pkv_PTHRGetScratchMem;
  pkv_FreeScratchMem   = _pkv_PTHRFreeScratchMem;
  pkv_GetScratchMemTop = _pkv_PTHRGetScratchMemTop;
  pkv_SetScratchMemTop = _pkv_PTHRSetScratchMemTop;
  pkv_ScratchMemAvail  = _pkv_PTHRScratchMemAvail;
  pkv_MaxScratchTaken  = _pkv_PTHRMaxScratchTaken;
  pkv_threads_in_use = true;
  return true;

failure:
  if ( thrstacks ) PKV_FREE ( thrstacks );
  if ( pkvthread ) PKV_FREE ( pkvthread );
  pkv_threads_in_use = false;
  return false;
} /*pkv_InitPThreads*/

/* to be called from the main thread */
void pkv_DestroyPThreads ( void )
{
  pkv_CancelPThreads ();
  pkv_threads_in_use = false;
  PKV_FREE ( pkvthread );
  PKV_FREE ( thrstacks );
  max_threads = 0;
  pthread_mutex_destroy ( &thread_mutex );
  pthread_mutex_destroy ( &cpustack_mutex );
        /* restore default scratch memory procedures */
  _pkv_AssignDefaultScratchMemProc ();
} /*pkv_DestroyPThreads*/

/* ////////////////////////////////////////////////////////////////////////// */
static short _alloc_thr ( void )
{
  short pos;

  pthread_mutex_lock ( &thread_mutex );
  if ( pkvt_freelist >= 0 ) {
        /* take element from the free list if nonempty */
    pos = pkvt_freelist;
    pkvt_freelist = pkvthread[pos].next;
  }
  else if ( pkvt_firstfree >= max_threads ) {
        /* too many threads at the same moment */
    pthread_mutex_unlock ( &thread_mutex );
    return -1;
  }
  else
    pos = pkvt_firstfree ++;
  memset ( &pkvthread[pos].thread, 0, sizeof(pthread_t) );
  pkvthread[pos].valid = false;
  pkvthread[pos].next = pkvt_busylist;
  pkvt_busylist = pos;
  pthread_mutex_unlock ( &thread_mutex );
  return pos;
} /*_alloc_thr*/

static void _free_thr ( short pos )
{
  short *nextp;

  pthread_mutex_lock ( &thread_mutex );
        /* extract the element from the busy list */
  for ( nextp = &pkvt_busylist; *nextp != pos; nextp = &pkvthread[*nextp].next )
    if ( *nextp < 0 ) {
      pthread_mutex_unlock ( &thread_mutex );
      return;
    }
  *nextp = pkvthread[pos].next;
        /* insert it to the free list */
  memset ( &pkvthread[pos].thread, 0, sizeof(pthread_t) );
  pkvthread[pos].valid = false;
  pkvthread[pos].next = pkvt_freelist;
  pkvt_freelist = pos;
  pthread_mutex_unlock ( &thread_mutex );
} /*_free_thr*/

/* ////////////////////////////////////////////////////////////////////////// */
static short int _pkv_NewThread ( int detachstate, PKVThreadProc startproc,
                                  size_t stacksize, size_t scratchmemsize,
                                  PKVThreadWorkToDo jobproc, void *jobdata,
                                  int4 *jobnum, void *auxdata,
                                  pthread_t *thread )
{
  pthread_attr_t attr;
  int            rc;
  short          pos;
  size_t         dstk;

  pos = _alloc_thr ();
  if ( pos >= 0 ) {
        /* create the scratch memory pool, if required */
    if ( scratchmemsize ) {
      PKV_MALLOC ( pkvthread[pos].ScratchPtr, scratchmemsize );
      if ( !pkvthread[pos].ScratchPtr ) {
        _free_thr ( pos );
        return -1;
      }
      pkvthread[pos].FreeScratchPtr = pkvthread[pos].ScratchPtr;
    }
    else
      pkvthread[pos].ScratchPtr = NULL;
    pkvthread[pos].ScratchSize = pkvthread[pos].FreeScratchSize =
    pkvthread[pos].MinFreeScratch = scratchmemsize;
        /* initialise the thread description structure */
    pkvthread[pos].jobproc = jobproc;
    pkvthread[pos].jobdata = jobdata;
    if ( jobnum )
      pkvthread[pos].jobnum = *jobnum;
    else
      memset ( &pkvthread[pos].jobnum, 0, sizeof(int4) );
    pkvthread[pos].auxdata = auxdata;
    pkvthread[pos].waiting = false;
    pkvthread[pos].success = true;  /* to be changed in case of failure */
    rc = pthread_cond_init ( &pkvthread[pos].cond, NULL );
    if ( rc )
      goto failure1;
    rc = pthread_mutex_init ( &pkvthread[pos].mutex, NULL );
    if ( rc )  /* cannot create the necessary mutex */
      goto failure2;
        /* create the thread */
    pthread_attr_init ( &attr );
    pthread_attr_setdetachstate ( &attr, detachstate );
    pthread_attr_getstacksize ( &attr, &dstk );
    if ( stacksize > dstk )
      pthread_attr_setstacksize ( &attr, stacksize );
    rc = pthread_create ( &pkvthread[pos].thread, &attr,
                          startproc, (void*)&pkvthread[pos] );
    if ( rc ) {  /* thread creation process failed */
      pthread_mutex_destroy ( &pkvthread[pos].mutex );
      pthread_attr_destroy ( &attr );
failure2:
      pthread_cond_destroy ( &pkvthread[pos].cond );
failure1:
      if ( pkvthread[pos].ScratchPtr )
        PKV_FREE ( pkvthread[pos].ScratchPtr );
      _free_thr ( pos );
      return -1;
    }
    pkvthread[pos].valid = true;
    if ( thread )
      *thread = pkvthread[pos].thread;
    pthread_attr_destroy ( &attr );
    return pos;
  }
  else
    return -1;
} /*_pkv_NewThread*/

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    pthread_mutex_t queue_mutex;     /* mutex and condition variable */
    pthread_cond_t  queue_cond;      /* used to wait for a free worker */
    pthread_mutex_t complete_mutex;  /* mutex and condition variable */
    pthread_cond_t  complete_cond;   /* used to wait until the last job */
                                     /* is complete */
    boolean         waiting1,        /* true when waiting for a free thread */
                    waiting2,        /* true when waiting for all threads to */
                                     /* complete their jobs */
                    waiting3;        /* true when waiting for thread */
                                     /* termination */
    int             totaljobs,       /* total number of jobs */
                    jobsdone;        /* number of jobs done */
    short int       npos, nthr;      /* queue capacity and number of threads */
                                     /* the actual queue */
    short int       fr, en;          /* front, end and capacity */
    short int       qtab[1];
  } PThreadsQueue;

/* this procedure is called by each working thread after its creation */
/* and after completing its job */
static void _pkv_EnqueuePThread ( PThreadsQueue *q, int pos )
{
  pthread_mutex_lock ( &q->queue_mutex );
  q->qtab[q->en++] = pos;
  if ( q->en > q->npos )
    q->en = 0;
  if ( q->waiting1 ) {
    pthread_cond_signal ( &q->queue_cond );
    q->waiting1 = false;
  }
  pkvthread[pos].waiting = true;
  pthread_mutex_unlock ( &q->queue_mutex );
} /*_pkv_EnqueuePThread*/

/* this procedure is called by the supervising procedure when it still has */
/* some job to assign and it needs a free working thread for it */
static short int _pkv_DequeuePThread ( PThreadsQueue *q )
{
  short int pos;

  pthread_mutex_lock ( &q->queue_mutex );
  if ( q->fr == q->en ) {
    q->waiting1 = true;
    pthread_cond_wait ( &q->queue_cond, &q->queue_mutex );
  }
  pos = q->qtab[q->fr++];
  if ( q->fr > q->npos )
    q->fr = 0;
  pthread_mutex_unlock ( &q->queue_mutex );
  return pos;
} /*_pkv_DequeuePThread*/

/* the start routine for detached threads, after thread creation it will wait */
/* for a job to do, after completing it, it will wait for another job etc.    */
static void* _pkv_DetachedSR ( void *usrdata )
{
  pkv_thread        *thr;
  pthread_t         myself;
  int               pos;
  PKVThreadWorkToDo jobproc;
  void              *jobdata;
  int4              jobnum;
  PThreadsQueue     *q;

  thr = (pkv_thread*)usrdata;
  thr->valid = true;
  pos = _pkv_PThreadMyPos ( &myself );
  RegisterThreadStack ( pos, (char*)&thr );
  q = (PThreadsQueue*)thr->auxdata;
  thr->jobproc = NULL;

        /* do subsequent jobs in a loop */
  for (;;) {
        /* stand in the queue */
    _pkv_EnqueuePThread ( q, pos );
          /* wait until a job to do arrives */
    pthread_mutex_lock ( &thr->mutex );
    if ( !thr->jobproc ) {
      thr->waiting = true;
      pthread_cond_wait ( &thr->cond, &thr->mutex );
    }
    pthread_mutex_unlock ( &thr->mutex );
    jobproc = thr->jobproc;
    thr->jobproc = NULL;
    jobdata = thr->jobdata;
    jobnum = thr->jobnum;
          /* do the job */
    thr->success = jobproc ( jobdata, &jobnum );
          /* signal the supervisor */
    pthread_mutex_lock ( &q->complete_mutex );
    q->jobsdone ++;
    if ( q->jobsdone == q->totaljobs && q->waiting2 ) {
      pthread_cond_signal ( &q->complete_cond );
      q->waiting2 = false;
    }
    pthread_mutex_unlock ( &q->complete_mutex );
  }
        /* this routine does not return */
  return NULL;
} /*_pkv_DetachedSR*/

boolean _pkv_BlackPill ( void *data, int4 *jobnum )
{
  PThreadsQueue *q;
  short int     pos;

  q = (PThreadsQueue*)data;
  pos = _pkv_PThreadMyPos ( NULL );
  WithdrawThreadStack ( pos );
  pthread_mutex_destroy ( &pkvthread[pos].mutex );
  pthread_cond_destroy ( &pkvthread[pos].cond );
  if ( pkvthread[pos].ScratchPtr )
    PKV_FREE ( pkvthread[pos].ScratchPtr );
  pthread_mutex_lock ( &q->complete_mutex );
  _free_thr ( pos );
  q->nthr --;
  if ( q->waiting3 && q->nthr == 0 ) {
    pthread_cond_signal ( &q->complete_cond );
    q->waiting3 = false;
  }
  pthread_mutex_unlock ( &q->complete_mutex );
  pthread_exit ( NULL );
} /*_pkv_BlackPill*/

boolean pkv_SetPThreadsToWork ( int jobdim, int4 *jobsize, int npthreads,
                                size_t stacksize, size_t scratchmemsize,
                                void *usrdata,
                                PKVThreadWorkToDo jobproc,
                                void *extradata,
                                PKVThreadWorkToDo extrajob,
                                boolean *success )
{
  PThreadsQueue *q;
  short int     *tpos;
  pthread_t     thread;
  int           rc, i, pos, njobs;
  int4          jobnum;
  boolean       _success;

  if ( jobdim < 1 || jobdim > 4 )
    exit ( 1 );
  njobs = jobsize->x;
  for ( i = 1; i < jobdim; i++ )
    njobs *= ((int*)jobsize)[i];
  if ( njobs < npthreads )
    npthreads = njobs;
        /* create and initialise the queue */
  PKV_MALLOC ( q, sizeof(PThreadsQueue)+(2*npthreads+1)*sizeof(short int) );
  if ( !q )
    return false;
  tpos = &q->qtab[npthreads+1];
  q->npos = q->nthr = npthreads;
  q->fr = q->en = 0;
  q->totaljobs = njobs;
  q->jobsdone = 0;
  q->waiting1 = q->waiting2 = false;
  rc = pthread_mutex_init ( &q->queue_mutex, NULL );
  if ( rc )
    goto failure1;
  rc = pthread_cond_init ( &q->queue_cond, NULL );
  if ( rc )
    goto failure2;
  rc = pthread_mutex_init ( &q->complete_mutex, NULL );
  if ( rc )
    goto failure3;
  rc = pthread_cond_init ( &q->complete_cond, NULL );
  if ( rc )
    goto failure4;
        /* create the working threads - each new thread goes to the queue */
  for ( i = 0; i < npthreads; i++ ) {
    tpos[i] = _pkv_NewThread ( PTHREAD_CREATE_DETACHED, _pkv_DetachedSR,
                              stacksize, scratchmemsize, NULL, NULL, NULL,
                              (void*)q, &thread );
    if ( tpos[i] < 0 ) {
      npthreads = q->nthr = i;
      _success = false;
      goto dismiss;
    }
  }
        /* assign and supervise the work */
  _success = true;
  memset ( &jobnum, 0, jobdim*sizeof(int) );
  do {
      /* wait until there is a thread in the queue */
    pos = _pkv_DequeuePThread ( q );
      /* find out, if so far the work went o.k. */
    if ( !pkvthread[pos].success ) {
        /* if not, make the thread terminate and go out */
      pthread_mutex_lock ( &pkvthread[pos].mutex );
      pkvthread[pos].jobproc = _pkv_BlackPill;
      pkvthread[pos].jobdata = (void*)q;
      if ( pkvthread[pos].waiting ) {
        pthread_cond_signal ( &pkvthread[pos].cond );
        pkvthread[pos].waiting = false;
      }
      pthread_mutex_unlock ( &pkvthread[pos].mutex );
      npthreads --;
      _success = false;
      goto dismiss;
    }
      /* set the thread to work */
    pthread_mutex_lock ( &pkvthread[pos].mutex );
    pkvthread[pos].jobdata = usrdata;
    pkvthread[pos].jobnum = jobnum;
    pkvthread[pos].jobproc = jobproc;
    if ( pkvthread[pos].waiting ) {
      pthread_cond_signal ( &pkvthread[pos].cond );
      pkvthread[pos].waiting = false;
    }
    pthread_mutex_unlock ( &pkvthread[pos].mutex );
  } while ( pkv_IncMultiCounter ( jobdim, (int*)jobsize, (int*)&jobnum ) );
        /* do the extra job, if present */
  if ( extrajob )
    _success = extrajob ( extradata, NULL );

        /* wait until all jobs are complete */
  pthread_mutex_lock ( &q->complete_mutex );
  if ( q->jobsdone < q->totaljobs ) {
    q->waiting2 = true;
    pthread_cond_wait ( &q->complete_cond, &q->complete_mutex );
  }
  pthread_mutex_unlock ( &q->complete_mutex );
        /* check if everything went well */
  for ( i = 0; i < npthreads; i++ ) {
    if ( !pkvthread[tpos[i]].success ) {
      _success = false;
      break;
    }
  }

dismiss:
        /* destroy the threads, or rather let them destroy themselves */
        /* using the black pill */
  q->waiting3 = false;
  for ( i = 0; i < npthreads; i++ ) {
    pos = _pkv_DequeuePThread ( q );
    pthread_mutex_lock ( &pkvthread[pos].mutex );
    pkvthread[pos].jobproc = _pkv_BlackPill;
    pkvthread[pos].jobdata = (void*)q;
    if ( pkvthread[pos].waiting ) {
      pthread_cond_signal ( &pkvthread[pos].cond );
      pkvthread[pos].waiting = false;
    }
    pthread_mutex_unlock ( &pkvthread[pos].mutex );
  }
        /* wait until all working threads free their entries */
        /* in the pkvthread array */
  pthread_mutex_lock ( &q->complete_mutex );
  if ( q->nthr > 0 ) {
    q->waiting3 = true;
    pthread_cond_wait ( &q->complete_cond, &q->complete_mutex );
  }
  pthread_mutex_unlock ( &q->complete_mutex );

        /* destroy the queue */
  pthread_cond_destroy ( &q->complete_cond );
  pthread_mutex_destroy ( &q->complete_mutex );
  pthread_cond_destroy ( &q->queue_cond );
  pthread_mutex_destroy ( &q->queue_mutex );
  PKV_FREE ( q );
  if ( success )
    *success = _success;
#ifdef DEBUG
printf ( "\n" );
#endif
  return true;

failure4:
  pthread_mutex_destroy ( &q->complete_mutex );
failure3:
  pthread_cond_destroy ( &q->queue_cond );
failure2:
  pthread_mutex_destroy ( &q->queue_mutex );
failure1:
  PKV_FREE ( q );
  return false;
} /*pkv_SetPThreadsToWork*/

/* ////////////////////////////////////////////////////////////////////////// */
/* the start routine for joinable threads; after thread creation it will do   */
/* the job, then release the thread resources and terminate the thread */
static void* _pkv_JoinableSR ( void *usrdata )
{
  pkv_thread        *thr;
  pthread_t         myself;
  short int         pos;
  PKVThreadWorkToDo jobproc;
  void              *jobdata;
  int4              jobnum;

  thr = (pkv_thread*)usrdata;
  thr->valid = true;
  pos = _pkv_PThreadMyPos ( &myself );
  RegisterThreadStack ( pos, (char*)&thr );
  jobproc = thr->jobproc;
  jobdata = thr->jobdata;
  jobnum = thr->jobnum;
        /* do the job */
  jobproc ( jobdata, &jobnum );
        /* cleanup after the job is done and terminate */
  WithdrawThreadStack ( pos );
  pthread_mutex_destroy ( &pkvthread[pos].mutex );
  pthread_cond_destroy ( &pkvthread[pos].cond );
  if ( thr->ScratchPtr )
    PKV_FREE ( thr->ScratchPtr );
  _free_thr ( pos );
  pthread_exit ( NULL );
} /*_pkv_JoinableSR*/

short int pkv_NewJoinablePThread ( size_t stacksize, size_t scratchmemsize,
                                   PKVThreadWorkToDo jobproc, void *jobdata,
                                   int4 *jobnum, void *auxdata,
                                   pthread_t *thread )
{
  return _pkv_NewThread ( PTHREAD_CREATE_JOINABLE, _pkv_JoinableSR,
                          stacksize, scratchmemsize, jobproc, jobdata, jobnum,
                          auxdata, thread );
} /*pkv_NewJoinablePThread*/

/* ////////////////////////////////////////////////////////////////////////// */
void pkv_CancelPThread ( short int pos )
{
  if ( pos >= 0 ) {
    pthread_cancel ( pkvthread[pos].thread );
    pthread_mutex_destroy ( &pkvthread[pos].mutex );
    pthread_cond_destroy ( &pkvthread[pos].cond );
    if ( pkvthread[pos].ScratchPtr )
      PKV_FREE ( pkvthread[pos].ScratchPtr );
    WithdrawThreadStack ( pos );
    _free_thr ( pos );
  }
} /*pkv_CancelPThread*/

void pkv_CancelPThreads ( void )
{
  int pos;

  for ( pos = pkvt_busylist; pos >= 0; pos = pkvt_busylist )
    pkv_CancelPThread ( pos );
} /*pkv_CancelPThreads*/

