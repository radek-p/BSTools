
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "g2blendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "bsfile.h"

#include "pomnijipc.h"
#include "pomnij_proc.h"
#include "proc_regmem.h"

jmp_buf escape_place;
boolean jump_ready = false;
boolean got_interrupted = false;

ipc_options options =
  { -1.0, 20, 3, 8, {3, 1000, 3, 1000}, 0, true, false, false };

int degu, degv, lknu, lknv, dim;
double *knotsu = NULL, *knotsv = NULL, *cpoints = NULL;
boolean clu, clv;  /* clamped boundary indicators */

/* constraints data: curves of constant "u" parameter */
int nucknots;
double *ucknots = NULL, *ucurvcp = NULL;
/* constraint equations matrix and right hand side */
double *cmat = NULL, *crhs = NULL;

/* ////////////////////////////////////////////////////////////////////////// */
void pomnij_proc_callback ( int msg, int size )
{
  void    *sp;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  if ( got_interrupted )
    goto handle_interrupt;
  if ( setjmp ( escape_place ) == 0 )  {
    jump_ready = true;
/*
printf ( "msg = %d, size = %d\n", msg, size );
*/
    switch ( msg ) {
  case cmdGET_OPTIONS:
      if ( InputOptions ( size ) )
        xge_CallTheParent ( cmdSUCCESS, 0 );
      else
        xge_CallTheParent ( cmdERROR, 0 );
      break;

  case cmdGET_BSPATCH:
      if ( InputBSPatch ( size ) )
        xge_CallTheParent ( cmdSUCCESS, 0 );
      else {
          /* the input pipe ought to be emptied - TODO */
        xge_CallTheParent ( cmdERROR, 0 );
      }
      break;

  case cmdGET_CONSTRAINTS:
      if ( InputConstraints ( size ) )
        xge_CallTheParent ( cmdSUCCESS, 0 );
      else
        xge_CallTheParent ( cmdERROR, 0 );
      break;

  case cmdG2BL_OPTIMIZE_LMT:
/*
printf ( "child: init\n" );
*/
      switch ( options.gcont ) {
    case 1:
        result = G1OptimizeLMTInit ();
        break;
    case 2:
        result = G2OptimizeLMTInit ();
        break;
    default:
        result = false;
        break;
      }
      if ( result )
        xge_ChildCallYourself ( cmdCONTINUE_OPT );
      else
        xge_CallTheParent ( cmdERROR, 0 );
      break;

  case cmdCONTINUE_OPT:
/*
printf ( "child: continue\n" );
*/
        /* make one iteration */
      switch ( options.gcont ) {
    case 1:
        result = G1OptimizeLMTIter ();
        break;
    case 2:
        result = G2OptimizeLMTIter ();
        break;
    default:
        result = false;
        break;
      }
      if ( result ) {
        if ( finished )
          xge_CallTheParent ( cmdSUCCESS, 0 );
        else if ( options.send_partial )
          xge_CallTheParent ( cmdG2BL_PARTIAL_RESULT, 0 );
        else
          xge_ChildCallYourself ( cmdCONTINUE_OPT );
      }
      else
        xge_CallTheParent ( cmdERROR, 0 );
      break;

  case cmdSEND_BSPATCH:
/*
printf ( "child: send\n" );
*/
      xge_CallTheParent ( cmdG2BL_FINAL_RESULT, PatchDataSize () );
      xge_ChildCallYourself ( cmdSEND_THE_BSPATCH );
      break;

  case cmdSEND_CONSTRAINTS:
      xge_CallTheParent ( cmdG2BL_CONSTRAINTS, ConstraintsSize () );
      xge_ChildCallYourself ( cmdSEND_THE_CONSTRAINTS );
      break;

  case cmdSEND_THE_BSPATCH:
/*
printf ( "Begun outputting\n" );
*/
      OutputBSPatch ();
/*
printf ( "Finished outputting\n" );
*/
      break;

  case cmdSEND_THE_CONSTRAINTS:
      OutputConstraints ();
      break;

  case cmdINTERRUPT:  /* this command is sent by the parent in order to */
                      /* force return from XNextEvent in case the child */
                      /* is waiting for an event, and in any case */
                      /* it should be ignored */
      break;

  case cmdTERMINATE:
/*
printf ( "Child terminated\n" );
*/
      exit ( 0 );

  default:
      xge_CallTheParent ( cmdERROR, 0 );
      xge_ChildFlushPipe ();
      break;
    }
    jump_ready = false;
  }
  else {
handle_interrupt:
    xge_ChildFlushPipe ();
    xge_CallTheParent ( cmdGOT_INTERRUPT, 0 );
    got_interrupted = false;
  }
  pkv_SetScratchMemTop ( sp );
} /*pomnij_proc_callback*/

void my_signal_handler ( void )
{
  if ( pkv_critical ) {
    pkv_signal = true;
  }
  else {
    got_interrupted = true;
    pkv_critical = true;
    MemBlockFreeAll ();
      /* the memory blocks pointed by the global pointers have been  */
      /* deallocated; the pointers must be cleared in order to avoid */
      /* attempts of deallocating them again */
    knotsu = knotsv = cpoints = ucknots = ucurvcp = cmat = crhs = NULL;
    pkv_critical = false;
    if ( jump_ready ) {
      jump_ready = false;
      longjmp ( escape_place, 1 );
    }
  }
} /*my_signal_handler*/

#ifdef _FPU_CONTROL_H
void SetupFPEInterrupt ( void )
{
  unsigned short cw;

  _FPU_GETCW ( cw );
  cw &= ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  _FPU_SETCW ( cw );
} /*SetupFPEInterrupt*/
#endif

void pomnij_proc_signal_handler ( int sig )
{
  switch ( sig ) {
case SIGUSR1:
    signal ( SIGUSR1, pomnij_proc_signal_handler );
printf ( "Child interrupt\n" );
    my_signal_handler ();
    break;
case SIGFPE:
    signal ( SIGFPE, pomnij_proc_signal_handler );
    SetupFPEInterrupt ();
printf ( "Floating point exception\n" );
    my_signal_handler ();
    break;
default:
    break;
  }
} /*pomnij_proc_signal_handler*/

void lib_error_handler ( int module, const char *file, int line,
                         int errcode, const char *errstr )
{
  fprintf ( stderr, "Error in module %d, file %s, line %d: %s\n",
                    module, file, line, errstr );
  switch ( module ) {
case LIB_G2BLENDING:
    break;
default:
    exit ( 1 );
  }
} /*lib_error_handler*/

/* ////////////////////////////////////////////////////////////////////////// */
int main ( int argc, char **argv )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  if ( !xge_ChildInit ( argc, argv, POMNIJ_IPC_MAGIC,
                        pomnij_proc_callback ) ) {
    exit ( 1 );
  }
  if ( !pkv_InitScratchMem ( SCRATCH_MEM_SIZE ) ) {
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }
  pkv_SetErrorHandler ( lib_error_handler );
        /* from now on all memory blocks allocated by malloc */
        /* have to be deallocated if an interrupt comes */
  MemBlockListInit ();
  pkv_register_memblock = MemBlockRegister;
  pkv_signal_handler = my_signal_handler;
  signal ( SIGUSR1, pomnij_proc_signal_handler );
  signal ( SIGFPE, pomnij_proc_signal_handler );
#ifdef _FPU_CONTROL_H
  SetupFPEInterrupt ();
#endif
printf ( "New child process launched\n" );
  xge_CallTheParent ( cmdCHILD_READY, 0 );
  xge_ChildMessageLoop ();
  exit ( 0 );
} /*main*/

