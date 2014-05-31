
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
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


void xge_SetupScrollWidgetPos ( xge_widget *er )
{
  xge_scroll_widget *sw;
  xge_widget *contents, *clipw, *xsb, *ysb;
  short      ws, hs, wc, hc;

  sw = er->data0;
  contents = sw->contents;
  clipw = sw->clipw;
  xsb = sw->xsb;
  ysb = sw->ysb;
  clipw->x = er->x;
  clipw->y = er->y;
  ws = er->w;
  hs = er->h;
  wc = contents->w;  hc = contents->h;
  if ( wc <= ws ) {
    if ( hc <= hs ) {
      xge_HScrollBarSetPagePos ( xsb, 0, false );
      xge_VScrollBarSetPagePos ( ysb, 0, false );
    }
    else {
      if ( wc > ws-10 )
        goto both;
      ws -= 10;
      xge_HScrollBarSetPagePos ( xsb, 0, false );
    }
  }
  else {
    if ( hc <= hs-10 ) {
      hs -= 10;
      xge_VScrollBarSetPagePos ( ysb, 0, false );
    }
    else {
both:
      ws -= 10;
      hs -= 10;
    }
  }
  sw->xslon = ws < contents->w;
  sw->yslon = hs < contents->h;
  contents->x = er->x;
  contents->y = er->y;
  if ( sw->yslon )
    contents->y -= sw->ys.pagevpos;
  if ( sw->xslon )
    contents->x -= sw->xs.pagehpos;
  if ( sw->yslon ) {
    clipw->x++;  clipw->w--;  clipw->y++;
    if ( sw->xslon )
      hs--;
    else
      hs -= 2;
  }
  else {
    if ( sw->xslon ) {
      clipw->y++;  hs--;  clipw->x++;  ws -= 2;
    }
  }
  clipw->w = ws;
  clipw->h = hs;
  if ( sw->xslon || sw->yslon ) {
    contents->x += 2;
    contents->y += 2;
    if ( sw->xslon ) {
      xsb->w = ws;
      xge_HScrollBarResizeWindow ( xsb, ws );
    }
    if ( sw->yslon ) {
      ysb->h = hs;
      xge_VScrollBarResizeWindow ( ysb, hs );
    }
  }
  contents->msgproc ( contents, xgemsg_RESIZE, 0, contents->w, contents->h );
} /*xge_SetupScrollWidgetPos*/

void xge_DrawScrollWidget ( xge_widget *er, boolean onscreen )
{
  xge_scroll_widget *sw;

  sw = er->data0;
  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );
  xge_SetClipping ( sw->contents );
  sw->contents->redraw ( sw->contents, false );
  if ( sw->xslon || sw->yslon ) {
    xge_SetClipping ( er );
    xgeSetForeground ( xgec_Grey2 );
    xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
    if ( sw->xslon ) {
      xge_SetClipping ( sw->xsb );
      sw->xsb->redraw ( sw->xsb, false );
    }
    if ( sw->yslon ) {
      xge_SetClipping ( sw->ysb );
      sw->ysb->redraw ( sw->ysb, false );
    }
  }
  if ( onscreen ) {
    xge_SetClipping ( er );
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
  }
} /*xge_DrawScrollWidget*/

boolean xge_ScrollWidgetMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_scroll_widget *sw;
  xge_widget *clipw;

  sw = er->data0;
  clipw = sw->clipw;
  switch ( msg ) {
case xgemsg_HSCROLLBAR_COMMAND:
case xgemsg_VSCROLLBAR_COMMAND:
    xge_SetupScrollWidgetPos ( er );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return true;

case xgemsg_RESIZE:
    clipw->x = er->x;
    clipw->y = er->y;
    er->w = clipw->w = x;
    if ( sw->yslon ) clipw->w -= 10;
    er->h = clipw->h = y;
    if ( sw->xslon ) clipw->h -= 10;
    sw->xsb->y = er->y+er->h-10;
    sw->ysb->x = er->x+er->w-10;
    xge_SetupScrollWidgetPos ( er );
    if ( key ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;

case xgemsg_KEY:
case xgemsg_SPECIAL_KEY:
case xgemsg_MMOVE:
case xgemsg_MCLICK:
case xgemsg_OTHEREVENT:
    if ( sw->xslon ) {
      if ( xge_PointInRect ( sw->xsb, x, y ) )
        return sw->xsb->msgproc ( sw->xsb, msg, key, x, y );
    }
    if ( sw->yslon ) {
      if ( xge_PointInRect ( sw->ysb, x, y ) )
        return sw->ysb->msgproc ( sw->ysb, msg, key, x, y );
    }
    if ( xge_PointInRect ( clipw, x, y ) &&
         xge_PointInRect ( sw->contents, x, y ) )
      return sw->contents->msgproc ( sw->contents, msg, key, x, y );
    return false;

default:
    return false;
  }
} /*xge_ScrollWidgetMsg*/

xge_widget *xge_NewScrollWidget ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y,
                                  xge_scroll_widget *sw, xge_widget *contents )
{
  xge_widget *nsw, *clipw, *xsb, *ysb;

  nsw = xge_NewWidget ( window_num, prev, id, w, h, x, y, sw, contents,
                        xge_ScrollWidgetMsg, xge_DrawScrollWidget );
  sw->er = nsw;
  if ( nsw ) {
    clipw = sw->clipw = xge_NewEmptyWidget ( window_num, NULL, 0, w, h, x, y );
    xsb = sw->xsb = xge_NewHScrollBar ( window_num, clipw, 1, w, 10, x, y+h-10,
                                        contents->w, clipw->w, &sw->xs, nsw );
    ysb = sw->ysb = xge_NewVScrollBar ( window_num, xsb, 2, 10, h, x+w-10, y,
                                        contents->h, clipw->h, &sw->ys, nsw );
    if ( clipw && xsb && ysb ) {
      xsb->up = ysb->up = clipw->up = nsw;
      contents->up = clipw;
      sw->xs.slpos = sw->ys.slpos = 0;
      sw->xslon = sw->yslon = false;
      sw->contents = contents;
      xge_SetupScrollWidgetPos ( nsw );
    }
    else
      return NULL;
  }
  return nsw;
} /*xge_NewScrollWidget*/

