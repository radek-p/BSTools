
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

/* ///////////////////////////////////////////////////////////////////////// */
/* four-window widget, allows resizing four windows, tiled in its rectangle  */

/* minimal widget width and height */
#define MINW  20
#define MINH  20
#define HDIST  2

static int xge_GetFourWW ( xge_widget *er, xge_widget **ww )
{
  int cnt;
  xge_widget *w;

  for ( cnt = 0, w = er->data0;  cnt < 4;  cnt++ ) {
    if ( (ww[cnt] = w) )
      w = w->next;
    else
      break;
  }
  return cnt;
} /*xge_GetFourWW*/

boolean xge_CompSizeFourWW ( xge_widget *er, char cs )
{
  xge_widget *fww[4];
  xge_fourww *fwwdata;
  short      splitx, splity, w, h, x0, x1, y0, y1, w0, w1, h0, h1, z;
  float      xsfr, ysfr;
  boolean    resized;

  fwwdata = er->data1;
  w = er->w;  h = er->h;
  z = fwwdata->zoomwin;
  if ( z >= 0 && z < 4 ) {
    fwwdata->win[z]->w = er->w;
    fwwdata->win[z]->h = er->h;
    fwwdata->win[z]->x = er->x;
    fwwdata->win[z]->y = er->y;
    fwwdata->resized = true;
    return true;
  }

  x0 = er->x;  y0 = er->y;
  switch ( cs ) {
case 0:
    splitx = fwwdata->splitx;  splity = fwwdata->splity;
    w0 = (short)(splitx-x0-HDIST);
    h0 = (short)(splity-y0-HDIST);
    break;
case 1:
    xsfr = fwwdata->xsfr;  w0 = (short)(xsfr*w-HDIST);
    ysfr = fwwdata->ysfr;  h0 = (short)(ysfr*h-HDIST);
    break;
default:
    return false;
  }
  w0 = (short)max ( w0, MINW );  w0 = (short)min ( w0, w-(MINW+2*HDIST) );
  h0 = (short)max ( h0, MINH );  h0 = (short)min ( h0, h-(MINH+2*HDIST) );
  w1 = (short)(w-w0-(2*HDIST));
  h1 = (short)(h-h0-(2*HDIST));
  x1 = (short)(x0+w0+(2*HDIST));
  y1 = (short)(y0+h0+(2*HDIST));
  fwwdata->splitx = (short)(x0+w0+HDIST);
  fwwdata->splity = (short)(y0+h0+HDIST);
  switch ( cs ) {
case 0:
    fwwdata->xsfr = (float)(w0+HDIST)/(float)w;
    fwwdata->ysfr = (float)(h0+HDIST)/(float)h;
    break;
case 1:
    break;
  }
  xge_GetFourWW ( er, fww );
  resized = false;
  if ( fww[0] ) {
    if ( fww[0]->w != w1 || fww[0]->h != h1 ||
         fww[0]->x != x1 || fww[0]->y != y1 ) resized = true;
    fww[0]->w = w1;  fww[0]->h = h1;  fww[0]->x = x1;  fww[0]->y = y1;
  }
  if ( fww[1] ) {
    if ( fww[1]->w != w0 || fww[1]->h != h1 ||
         fww[1]->x != x0 || fww[1]->y != y1 ) resized = true;
    fww[1]->w = w0;  fww[1]->h = h1;  fww[1]->x = x0;  fww[1]->y = y1;
  }
  if ( fww[2] ) {
    if ( fww[2]->w != w1 || fww[2]->h != h0 ||
         fww[2]->x != x1 || fww[2]->y != y0 ) resized = true;
    fww[2]->w = w1;  fww[2]->h = h0;  fww[2]->x = x1;  fww[2]->y = y0;
  }
  if ( fww[3] ) {
    if ( fww[3]->w != w0 || fww[3]->h != h0 ||
         fww[3]->x != x0 || fww[3]->y != y0 ) resized = true;
    fww[3]->w = w0;  fww[3]->h = h0;  fww[3]->x = x0;  fww[3]->y = y0;
  }
  return (fwwdata->resized = resized);
} /*xge_CompSizeFourWW*/

void xge_DrawFourWW ( xge_widget *er, boolean onscreen )
{
  xge_widget *fww[4];
  xge_fourww *fwwdata;
  int i, nww;

  nww = xge_GetFourWW ( er, fww );
  fwwdata = er->data1;
  i = fwwdata->zoomwin;
  if ( i >= 0 && i < 4 ) {
    if ( fwwdata->resized ) {
      xge_SetClipping ( fwwdata->win[i] );
      xge_CallMsgProc ( (xge_widget**)&er->data2, fwwdata->win[i],
                        xgemsg_RESIZE, 0, fwwdata->win[i]->w, fwwdata->win[i]->h );
    }
    xge_SetClipping ( fwwdata->win[i] );
    fwwdata->win[i]->redraw ( fwwdata->win[i], false );
    goto finish_it;
  }
  if ( fwwdata->resized ) {
    for ( i = 0; i < nww; i++ ) {
      xge_SetClipping ( fww[i] );
      xge_CallMsgProc ( (xge_widget**)&er->data2, fww[i],
                        xgemsg_RESIZE, 0, fww[i]->w, fww[i]->h );
    }
    fwwdata->resized = false;
  }
  for ( i = 0; i < nww; i++ ) {
    xge_SetClipping ( fww[i] );
    fww[i]->redraw ( fww[i], false );
  }
  xge_SetClipping ( er );
  xgeSetForeground ( xgec_SPLIT_WIN_BK );
  xgeFillRectangle ( er->w, 2*HDIST, er->x, fwwdata->splity-HDIST );
  xgeFillRectangle ( 2*HDIST, er->h, fwwdata->splitx-HDIST, er->y );
finish_it:
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawFourWW*/

boolean xge_FourWWMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_fourww *fww;
  xge_widget *ww;
  boolean    cx, cy;
  int        i;

  fww = er->data1;
  switch ( er->state ) {
case xgestate_RESIZING_X:
    switch ( msg ) {
  case xgemsg_MCLICK:
      goto process_click;
  case xgemsg_MMOVE:
      fww->splitx = x;
      if ( xge_CompSizeFourWW ( er, 0 ) )
        goto process_resizing;
  case xgemsg_KEY:
      goto process_key;
  default:
      break;
    }
    break;

case xgestate_RESIZING_Y:
    switch ( msg ) {
  case xgemsg_MCLICK:
      goto process_click;
  case xgemsg_MMOVE:
      fww->splity = y;
      if ( xge_CompSizeFourWW ( er, 0 ) )
        goto process_resizing;
  case xgemsg_KEY:
      goto process_key;
  default:
      break;
    }
    break;

case xgestate_RESIZING_XY:
    switch ( msg ) {
  case xgemsg_MCLICK:
process_click:
      if ( !(key & xgemouse_LBUTTON_DOWN) ) {
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
      }
      break;
  case xgemsg_MMOVE:
      fww->splitx = x;
      fww->splity = y;
      if ( xge_CompSizeFourWW ( er, 0 ) )
        goto process_resizing;
  case xgemsg_KEY:
      goto process_key;
  default:
      break;
    }
    break;

default:
    switch ( msg ) {
  case xgemsg_ENTERING:
      er->data2 = NULL;
      xge_SetCurrentWindowCursor ( xgeCURSOR_CROSSHAIR );
      return true;

  case xgemsg_EXITING:
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
      xge_CallMsgProc ( (xge_widget**)&er->data2, NULL, xgemsg_NULL, 0, 0, 0 );
      return true;

  case xgemsg_RESIZE:
      er->w = x;  er->h = y;
      if ( xge_CompSizeFourWW ( er, 1 ) ) {
process_resizing:
        ww = er->data0;
        while ( ww ) {
          xge_CallMsgProc ( (xge_widget**)&er->data2, ww,
                            xgemsg_RESIZE, 0, ww->w, ww->h );
          ww = ww->next;
        }
        fww->resized = false;
        if ( key ) {
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
      }
      return true;

  case xgemsg_MOVE:
      return true;

  case xgemsg_KEY:
      switch ( key ) {
    case 'S':  case 's':
        if ( fww->zoomwin == -1 ) {
          for ( i = 0; i < 4; i++ )
            if ( xge_PointInRect ( fww->win[i], x, y ) ) {
              fww->zoomwin = i;
              break;
            }
        }
        else
          fww->zoomwin = -1;
        xge_CompSizeFourWW ( er, 1 );
        goto process_resizing;
    default:
        break;
      }
      break;

  default:
      break;
    }
                    /* try to send it to one of your widgets */
    if ( fww->zoomwin >= 0 && fww->zoomwin < 4 ) {
      ww = fww->win[(int)fww->zoomwin];
      return xge_CallMsgProc ( (xge_widget**)&er->data2, ww, msg, key, x, y );
    }
    else {
      ww = er->data0;
      while ( ww ) {
        if ( xge_PointInRect ( ww, x, y ) )
          return xge_CallMsgProc ( (xge_widget**)&er->data2, ww, msg, key, x, y );
        ww = ww->next;
      }
    }
                    /* process the message sent to you */
    if ( (xge_widget**)&er->data2 ) {
      xge_CallMsgProc ( (xge_widget**)&er->data2, NULL, xgemsg_NULL, 0, 0, 0 );
      xge_SetCurrentWindowCursor ( xgeCURSOR_CROSSHAIR );
    }
    switch ( msg ) {
  case xgemsg_MMOVE:
      break;

  case xgemsg_MCLICK:
      if ( (er->state == xgestate_NOTHING) &&
           (key & xgemouse_LBUTTON_DOWN) ) {
        cx = (boolean)(abs(x-fww->splitx) <= 10);
        cy = (boolean)(abs(y-fww->splity) <= 10);
        if ( cx ) {
          if ( cy )
            er->state = xgestate_RESIZING_XY;
          else
            er->state = xgestate_RESIZING_X;
        }
        else if ( cy )
          er->state = xgestate_RESIZING_Y;
      }
      if ( er->state != xgestate_NOTHING )
        xge_GrabFocus ( er, true );
      break;

  case xgemsg_KEY:
process_key:
      switch ( key ) {
    case 'R':  case 'r':
        fww->xsfr = fww->ysfr = 0.5;
        fww->zoomwin = -1;
        if ( xge_CompSizeFourWW ( er, 1 ) ) {
          xge_dispatch_message ( xgemsg_MMOVE, 0, xge_prevx, xge_prevy );
          goto process_resizing;
        }
        break;
    default:
        break;
      }
      break;

  default:
      break;
    }
    break;
  }
  return true;
} /*xge_FourWWMsg*/

xge_widget *xge_NewFourWW ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            xge_widget *ww, xge_fourww *fwwdata )
{
  xge_widget *fww, *first;
  int        i;

  first = ww;
  if ( first )
    while ( first->prev )
      first = first->prev;
  if ( (fww = xge_NewWidget ( window_num, prev, id, w, h, x, y, first, fwwdata,
                              xge_FourWWMsg, xge_DrawFourWW )) ) {
    fwwdata->er = fww;
    fww->data2 = NULL;
    fwwdata->xsfr = fwwdata->ysfr = 0.5;
    fwwdata->zoomwin = -1;
    xge_CompSizeFourWW ( fww, 1 );
    for ( i = 0;  i < 4 && first;  i++, first = first->next )
      fwwdata->win[3-i] = first;
  }
  return fww;
} /*xge_NewFourWW*/

