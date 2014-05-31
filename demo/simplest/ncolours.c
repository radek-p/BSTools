
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

#include "ncolours.h"

#define CSW     130
#define CSH      75
#define HDIST   150
#define VDIST   100
#define HMARGIN  10
#define LMARGIN  20

#define STATE_SCROLLING_WIN (xgestate_LAST+1)

xge_widget *cwin, *menu, *scb;

int ycmin;
int wheight;
int nrows;  /* number of rows */
int ncols;  /* number of columns */
int nfr;  /* number of the first row */
short lasty;

xge_vscrollbar vsb;

/* ///////////////////////////////////////////////////////////////////////// */
void DrawColourSample ( int nc, short x, short y )
{
  int w;

  if ( nc >= 0 && nc < XGE_PALETTE_LENGTH ) {
    xgeSetForeground ( xge_palette[nc] );
    xgeFillRectangle ( CSW-1, CSH-1, x, y );
    xgeSetForeground ( xgec_White );
    w = strlen(xge_colour_name[nc]);
    xgeDrawString ( xge_colour_name[nc], x + (CSW-6*w)/2, y + CSH+13 );
  }
} /*DrawColourSample*/

void RysujOkno ( xge_widget *er, boolean onscreen )
{
  int   i, j, nlr, nc;
  short x, y;

  xgeSetForeground ( xgec_Black );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );

  nrows = (XGE_PALETTE_LENGTH+ncols-1) / ncols;
  wheight = 2*HMARGIN+nrows*VDIST;
  nfr = (HMARGIN-ycmin)/VDIST;
  nlr = nfr + er->h/VDIST + 2;
  for ( i = nfr, nc = nfr*ncols, y = (short)(er->y+ycmin+nfr*VDIST+HMARGIN);
        i < nlr;
        i++, y += VDIST )
    for ( j = 0, x = (short)(er->x+LMARGIN);  j < ncols;  j++, nc++, x += HDIST )
      DrawColourSample ( nc, x, y );

  xgeSetForeground ( xgec_CornflowerBlue );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*RysujOkno*/

boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  int lycmin;

  switch ( er->state ) {
case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_DOWN) &&
           (key & xgemouse_LBUTTON_CHANGE) ) {
        lasty = y;
        er->state = STATE_SCROLLING_WIN;
        xge_GrabFocus ( er, true );
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        lycmin = ycmin;
        ycmin += 25;
        ycmin = min ( ycmin, 0 );
        goto redraw_it;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        lycmin = ycmin;
        ycmin -= 25;
        ycmin = max ( ycmin, (er->h-wheight) );
        goto redraw_it;
      }
      return true;

 default:
      return false;
    }

case STATE_SCROLLING_WIN:
    switch ( msg ) {
  case xgemsg_MMOVE:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y != lasty ) {
          lycmin = ycmin;
          ycmin += y-lasty;
          ycmin = min ( ycmin, 0 );
          ycmin = max ( ycmin, (er->h-wheight) );
          lasty = y;
redraw_it:
          if ( ycmin != lycmin ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
            wheight = 2*HMARGIN+nrows*VDIST;
            xge_VScrollBarSetPagePos ( scb, -ycmin, true );
          }
        }
      }
      else {
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
      }
      return true;

  case xgemsg_MCLICK:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
      }
      return true;

  default:
      return false;
    }

default:
    return false;
  }
} /*OknoMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int lycmin;

  if ( er ) {
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id = 0 ) {
  case 0:
        xge_done = 1;
        break;

  default:
        break;
      }
      break;

case xgemsg_VSCROLLBAR_COMMAND:
      if ( er->id == 1 ) {
        wheight = 2*HMARGIN+nrows*VDIST;
        lycmin = ycmin;
        ycmin = -vsb.pagevpos;
        if ( ycmin != lycmin ) {
          xge_SetClipping ( cwin );
          cwin->redraw ( cwin, true );
        }
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
      cwin->w = scb->x = (short)(xge_current_width-10);
      menu->w = xge_current_width;
      cwin->h = scb->h = (short)(xge_current_height-20);
      ncols = cwin->w / HDIST;
      nrows   = (XGE_PALETTE_LENGTH+ncols-1) / ncols;
      wheight = 2*HMARGIN+nrows*VDIST;
      xge_VScrollBarResizeWindow ( scb, cwin->h );
      xge_VScrollBarResizePage ( scb, wheight );
      ycmin = -vsb.pagevpos;
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

void init_edwin ( void )
{
  xge_widget *rp;

  ncols   = (xge_WIDTH-10) / HDIST;
  nrows   = (XGE_PALETTE_LENGTH+ncols-1) / ncols;
  wheight = 2*HMARGIN+nrows*VDIST;
  ycmin   = 0;

  scb = xge_NewVScrollBar ( 0, NULL, 1, 10, xge_HEIGHT-20, xge_WIDTH-10, 20,
                            wheight, xge_HEIGHT-20, &vsb, NULL );
  rp = xge_NewButton ( 0, NULL, 0, 60, 18, 0, 0, b0 );
  menu = xge_NewMenu ( 0, scb, 2, xge_WIDTH, 20, 0, 0, rp );
  cwin = xge_NewWidget ( 0, menu, 0, xge_WIDTH-10, xge_HEIGHT-20, 0, 20,
                        NULL, NULL, OknoMsg, RysujOkno );
  nfr     = 0;

  xge_SetWinEdRect ( cwin );
  xge_Redraw ();
} /*init_edwin*/

void destroy_edwin ()
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

