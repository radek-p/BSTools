
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "xgledit.h"

#include "glsimplest.h"

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xgleClearColor3fv ( xglec_Blue6 );
  glClear ( GL_COLOR_BUFFER_BIT );
  xgle_SetIdentMapping ( er );
  glColor3fv ( xglec_Yellow );
  xgleDrawString ( "OpenGL-rendered text using an X font", 120, 20 );
  glFinish ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );
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

  default:
        break;
      }
      break;

default:
      break;
    }
  }
  else {
/*    printf ( "%d, %d, %d, %d\n", msg, key, x, y ); */
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

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
static char b0[] = "Quit";

void init_edwin ( void )
{
  xge_widget *rp, *edr;
  int dw, dh, dwmm, dhmm;

  rp = xge_NewButton ( 0, NULL, 0, 60, 18, 0, 0, b0 );
  edr = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, rp );
  edr = xge_NewWidget ( 0, edr, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                        NULL, NULL, OknoMsg, RysujOkno );

  xge_SetWinEdRect ( edr );
  xge_Redraw ();

  dw = DisplayWidth ( xgedisplay, xgescreen );
  dh = DisplayHeight ( xgedisplay, xgescreen );
  dwmm = DisplayWidthMM ( xgedisplay, xgescreen );
  dhmm = DisplayHeightMM ( xgedisplay, xgescreen );
  printf ( "%d, %d, %d, %d\n", dw, dh, dwmm, dhmm );
} /*init_edwin*/

void destroy_edwin ( void )
{
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xgle_Init ( argc, argv, CallBack, NULL, true, false, false );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xgle_Cleanup ();  
  exit ( 0 );
} /*main*/   

