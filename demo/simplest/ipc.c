
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "ipc.h"

boolean child_busy = false;

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xgeSetForeground ( xgec_Blue4 );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  xgeSetForeground ( xgec_CornflowerBlue );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*RysujOkno*/

boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y )
{
/*  printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
  return true;
} /*OknoMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
/*    printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case 0:
        xge_done = 1;
        break;

  case 1:
        if ( !child_busy ) {
          xge_CallTheChild ( 0, 0 );
          child_busy = true;
        }
        else
          xge_DisplayErrorMessage ( "Error: The child process is busy.", -1 );
        break;

  case 2:
        if ( child_busy )
          xge_SignalTheChild ();
        else
          xge_DisplayErrorMessage ( "Error: The child process is not busy.", -1 );
        break;

  default:
        break;
      }
      break;

default:
      break;
    }
  }
  else {
    switch ( msg ) {
case xgemsg_KEY:
      switch ( key ) {
  case 'Q': case 'q':
        xge_done = 1;
        break;
  default:
        break;
      }
      break;

case xgemsg_RESIZE:
      break;

case xgemsg_CHILD_MESSAGE:
      ProcessChildMessage ( xgeevent.xclient.data.l[0],
                            xgeevent.xclient.data.l[1] );
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
static char b0[] = "Quit";
static char b1[] = "Call";
static char b2[] = "Interrupt";

void init_edwin ( void )
{
  xge_widget *rp, *edr;

  rp = xge_NewButton ( 0, NULL, 0, 60, 18, 0, 0, b0 );
  rp = xge_NewButton ( 0, rp, 1, 60, 18, 0, 20, b1 );
  rp = xge_NewButton ( 0, rp, 2, 60, 18, 0, 40, b2 );
  edr = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, rp );
  edr = xge_NewWidget ( 0, edr, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                        NULL, NULL, OknoMsg, RysujOkno );

  xge_SetWinEdRect ( edr );
  xge_Redraw ();
} /*init_edwin*/

/* ////////////////////////////////////////////////////////////////////////// */
void ProcessChildMessage ( int msg, int size )
{
  printf ( "child message %d, size = %d\n", msg, size );
  child_busy = false;
} /*ProcessChildMessage*/

/* ////////////////////////////////////////////////////////////////////////// */
void destroy_edwin ( void )
{
void kill ( pid_t pid, int sig );
  kill ( xge_child_pid, SIGKILL );
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  if ( !xge_MakeTheChild ( argv[0], "_proc", 1 ) ) {
    xge_Cleanup ();  
    exit ( 1 );
  }
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xge_Cleanup ();  
  exit ( 0 );
} /*main*/   

