
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <setjmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "xgedit.h"
#include "xgeipc.h"


boolean got_a_signal = false;
jmp_buf loop_env;
boolean handler_ready = false;


void ReceiveSignal ( int sig )
{
  got_a_signal = true;
  signal ( SIGUSR1, ReceiveSignal );
  if ( handler_ready ) {
    handler_ready = false;
    longjmp ( loop_env, 1 );
  }
  /* else ignore it */
} /*ReceiveSignal*/

int count = 0;

void callback ( int msg, int size )
{
  if ( setjmp ( loop_env ) == 0 ) {
    handler_ready = true;
    switch ( msg ) {
  case 0:
      sleep ( 1 );
      xge_CallTheParent ( count++, 0 );
      break;
  default:
      break;
    }
  }
  else {
    xge_CallTheParent ( -1, 0 );
  }
} /*callback*/

int main ( int argc, char **argv )
{
  if ( !xge_ChildInit ( argc, argv, 1, callback ) ) {
        /* the parent process should be informed */
    exit ( 1 );
  }
  signal ( SIGUSR1, ReceiveSignal );
  xge_ChildMessageLoop ();
  exit ( 0 );
} /*main*/

