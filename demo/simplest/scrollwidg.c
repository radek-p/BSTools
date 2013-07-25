
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

#include "pkgeom.h"
#include "xgedit.h"

#include "scrollwidg.h"

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *okno, *menu, *scrollmenu;

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
  switch ( msg ) {
case xgemsg_RESIZE:
    er->w = x;
    er->h = y;
    if ( key ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;
default:
    return false;
  }
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
  case 2:
        printf ( "\a" );
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
      okno->msgproc ( okno, xgemsg_RESIZE, 0, x-110, y );
      menu->msgproc ( menu, xgemsg_RESIZE, 0, 110, y );
      scrollmenu->msgproc ( scrollmenu, xgemsg_RESIZE, 0, 110, y-30 );
      xge_Redraw ();
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
static char b0[] = "Quit";
static char b1[] = "Beep";

xge_scroll_widget sw;

void init_edwin ( void )
{
  xge_widget *bt1, *bt2, *cont;
/*  int dw, dh, dwmm, dhmm; */

  setvbuf ( stdout, NULL, _IONBF, 0 );
  bt1 = xge_NewButton ( 0, NULL, 2, 60, 18, 0, 40, b1 );
  xge_SetWidgetPositioning ( bt1, 0, 2, 2 );
  cont = xge_NewMenu ( 0, NULL, 3, 140, 400, 0, 0, bt1 );
  scrollmenu = xge_NewScrollWidget ( 0, NULL, 4, 110, 160, 0, 30, &sw, cont );

  bt2 = xge_NewButton ( 0, scrollmenu, 0, 60, 18, 0, 0, b0 );
  menu = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, bt2 );
  okno = xge_NewWidget ( 0, menu, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                         NULL, NULL, OknoMsg, RysujOkno );

  xge_SetWinEdRect ( okno );
  xge_Redraw ();
/*
  dw = DisplayWidth ( xgedisplay, xgescreen );
  dh = DisplayHeight ( xgedisplay, xgescreen );
  dwmm = DisplayWidthMM ( xgedisplay, xgescreen );
  dhmm = DisplayHeightMM ( xgedisplay, xgescreen );
  printf ( "%d, %d, %d, %d\n", dw, dh, dwmm, dhmm );
*/
} /*init_edwin*/

void destroy_edwin ( void )
{
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xge_Cleanup ();  
  exit ( 0 );
} /*main*/   

