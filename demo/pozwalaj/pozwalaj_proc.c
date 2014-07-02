
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "eg2holed.h"
#include "bsmesh.h"
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"
#include "mengerc.h"
#include "xgedit.h"

#define CHILD_SIDE
#include "pozwalajipc.h"

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void pozwalaj_proc_callback ( int msg, int size )
{
#ifdef DEBUG_IPC
printf ( "message for CHILD %d, size %d\n", msg, size );
#endif
  switch ( msg ) {
case ipccmd_GET_DATA:
    GetData ( size );
    break;

case ipccmd_BEGIN_BLP:
    BeginBLPOptimization ();
    break;

case ipccmd_CONTINUE_BLP1:
    ContinueBLPOptimization1 ();
    break;

case ipccmd_CONTINUE_BLP2:
    ContinueBLPOptimization2 ();
    break;

case ipccmd_BEGIN_BSM:
    BeginBSMOptimization ();
    break;

case ipccmd_CONTINUE_BSM:
    ContinueBSMOptimization ();
    break;

case ipccmd_BEGIN_BSCMC:
    BeginBSCMCOptimization ();
    break;

case ipccmd_CONTINUE_BSCMC:
    ContinueBSCMCOptimization ();
    break;

case ipccmd_SEND_RESULT:
    IPCSendData ();
    break;

default:
    break;
  }
} /*pozwalaj_proc_callback*/

void my_signal_handler ( void )
{
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

void pozwalaj_proc_signal_handler ( int sig )
{
  switch ( sig ) {
case SIGUSR1:
    signal ( SIGUSR1, pozwalaj_proc_signal_handler );
    my_signal_handler ();
    break;
case SIGFPE:
    signal ( SIGFPE, pozwalaj_proc_signal_handler );
#ifdef _FPU_CONTROL_H
    SetupFPEInterrupt ();
#endif
    my_signal_handler ();
printf ( "SIGFPU\n" );
exit ( 1 );
    break;
default:
    break;
  }
} /*pozwalaj_proc_signal_handler*/

int main ( int argc, char **argv )
{
  int ncpu;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  if ( !xge_ChildInit ( argc, argv, POZWALAJ_IPC_MAGIC, pozwalaj_proc_callback ) )
    exit ( 1 );
  if ( !pkv_InitScratchMem ( SCRATCH_MEM_SIZE ) ) {
    printf ( "Error: cannot allocate scratch memory stack\n" );
    exit ( 1 );
  }
  signal ( SIGUSR1, pozwalaj_proc_signal_handler );
  signal ( SIGFPE, pozwalaj_proc_signal_handler );
#ifdef _FPU_CONTROL_H
  SetupFPEInterrupt ();
#endif
  printf ( "New child process launched, %d\n", xge_child_pid );
  ncpu = pkv_FindNCPU ();
  if ( ncpu > 1 ) {
    printf ( "Found %d CPUs\n", ncpu );
    if ( pkv_InitPThreads ( 0x02 << G2MBL_MAX_LEVELS ) ) {
      printf ( "pthreads will be used\n" );
      _g2mbl_npthreads = ncpu;
    }
    else {
      printf ( "failed to prepare pthreads\n" );
    }
  }
  xge_CallTheParent ( ipccmd_CHILD_READY, 0 );
  xge_ChildMessageLoop ();
  exit ( 0 );
} /*main*/

