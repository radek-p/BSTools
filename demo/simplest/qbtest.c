
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkgeom.h"
#include "xgedit.h"

xge_quatrotballf qballf[9];
quaternionf      qf[9];
trans3f          qtrf[9];

xge_quatrotballd qballd[9];
quaterniond      qd[9];
trans3d          qtrd[9];

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      xge_done = 1;
      return 1;

case xgemsg_QUATROTBALL_COMMAND:
      return 1;

default:
      return 0;
    }
  }
  else {
 
    switch ( msg ) {
case xgemsg_KEY:
      switch ( key ) {
  case 'Q': case 'q':
        xge_done = 1;
        return 1;
  default: 
        return 0;
      }
default:
      return 0;
    }
  }
} /*CallBack*/

void init_edwin ( void  )
{
  xge_widget *w, *m;

  pkv_InitScratchMem ( 655360 );
  w = xge_NewButton ( 0, NULL, 9, 60, 18,
                      (xge_WIDTH-60)/2, (xge_HEIGHT-18)/2+60, "Quit" );
  w = xge_NewQuatRotBalld ( 0, w, 0, 80, 80, 5, 5, 80,
                            &qballd[0], &qd[0], &qtrd[0], NULL );
  w = xge_NewQuatRotBalld ( 0, w, 1, 80, 80, xge_WIDTH-85, 5, 80,
                            &qballd[1], &qd[1], &qtrd[1], NULL );
  w = xge_NewQuatRotBalld ( 0, w, 2, 80, 80, 5, xge_HEIGHT-85, 80,
                            &qballd[2], &qd[2], &qtrd[2], NULL );
  w = xge_NewQuatRotBalld ( 0, w, 3, 80, 80, xge_WIDTH-85, xge_HEIGHT-85, 80,
                            &qballd[3], &qd[3], &qtrd[3], NULL );
  w = xge_NewQuatRotBalld ( 0, w, 4, 80, 80, (xge_WIDTH-80)/2, (xge_HEIGHT-80)/2, 80,
                            &qballd[4], &qd[4], &qtrd[4], NULL );
  w = xge_NewQuatRotBallf ( 0, w, 5, 80, 80, (xge_WIDTH-80)/2, 5, 80,
                            &qballf[5], &qf[5], &qtrf[5], NULL );
  w = xge_NewQuatRotBallf ( 0, w, 6, 80, 80, xge_WIDTH-85, (xge_HEIGHT-80)/2, 80,
                            &qballf[6], &qf[6], &qtrf[6], NULL );
  w = xge_NewQuatRotBallf ( 0, w, 7, 80, 80, (xge_WIDTH-80)/2, xge_HEIGHT-85, 80,
                            &qballf[7], &qf[7], &qtrf[7], NULL );
  w = xge_NewQuatRotBallf ( 0, w, 8, 80, 80, 5, (xge_HEIGHT-80)/2, 80,
                            &qballf[8], &qf[8], &qtrf[8], NULL );
  m = xge_NewMenu ( 0, NULL, 10, xge_WIDTH, xge_HEIGHT, 0, 0, w );
  xge_SetWinEdRect ( m );
  xge_Redraw ();
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

