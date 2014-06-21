
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "xgedit.h"
#include "xgeprivate.h"


void xge_DrawSlidebard ( xge_widget *er, boolean onscreen )
{
  int    x;
  double *slipos;

  slipos = er->data0;
  if ( er->state == xgestate_MOVINGSLIDE )
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SLIDEBAR_A, xgec_SLIDEBAR_B, er->h-2 );
  else
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SLIDEBAR_C, xgec_SLIDEBAR_D, er->h-2 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  x = er->x + 2 + (int)((*slipos)*(double)(er->w - 10));
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  xgeFillRectangle ( 6, 6, x, er->y+2 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawSlidebard*/

boolean xge_SlidebardMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  double  z;
  double  *slipos;
  boolean inval;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  slipos = er->data0;
  inval = false;
  if ( er->state == xgestate_NOTHING ) {
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( x < er->x+5 ) x = (short)(er->x+5);
        else if ( x > er->x+er->w-5 ) x = (short)(er->x+er->w-5);
        er->state = xgestate_MOVINGSLIDE;
        xge_GrabFocus ( er, true );
        inval = true;
        goto update1;
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        z = *slipos - 1.0/(double)(er->w-10);
        z = max ( z, 0.0 );
        goto update2;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        z = *slipos + 1.0/(double)(er->w-10);
        z = min ( z, 1.0 );
        goto update2;
      }

  default:
      return false;
    }
  }
  else if ( er->state == xgestate_MOVINGSLIDE ) {
    switch ( msg ) {
  case xgemsg_MMOVE:
  case  xgemsg_MCLICK:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( x < er->x+5 ) x = (short)(er->x+5);
        else if ( x > er->x+er->w-5 ) x = (short)(er->x+er->w-5);
        if ( x != xge_xx ) {
          inval = false;
update1:
          xge_xx = x;
          z = (double)(x-er->x-5)/(double)(er->w-10);
update2:
          if ( *slipos != z ) {
            inval = true;
            *slipos = z;
                      /* double slidebar position is not passed as a parameter */
            if ( !er->up ||
                 !er->up->msgproc ( er->up, xgemsg_SLIDEBAR_COMMAND, er->id, x, y ) )
              xge_callback ( er, xgemsg_SLIDEBAR_COMMAND, 0, x, y );
          }
          if ( inval ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else {
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      return true;

  default:
      return false;
    }
  }
  return false;
} /*xge_SlidebardMsg*/

double xge_LinSlidebarValued ( double xmin, double xmax, double t )
{
  return xmin + t*(xmax-xmin);
} /*xge_LinSlidebarValued*/

double xge_LinSlidebarPosd ( double xmin, double xmax, double x )
{
  return (x-xmin)/(xmax-xmin);
} /*xge_LinSlidebarPosd*/

double xge_LogSlidebarValued ( double xmin, double xmax, double t )
{
  double a, b;

  b = log ( xmin );
  a = log ( xmax ) - b;
  return exp ( a*t + b );
} /*xge_LogSlidebarValued*/

double xge_LogSlidebarPosd ( double xmin, double xmax, double x )
{
  double a, b;

  b = log ( xmin );
  a = log ( xmax ) - b;
  return (log ( x ) - b)/a;
} /*xge_LogSlidebarPosd*/

xge_widget *xge_NewSlidebard ( char window_num, xge_widget *prev, int id,
                               short w, short h, short x, short y,
                               double *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, NULL,
                         xge_SlidebardMsg, xge_DrawSlidebard );
} /*xge_NewSlidebard*/

