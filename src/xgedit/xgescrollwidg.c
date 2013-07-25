
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
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
  xge_widget *contents, *clipw, *xsl, *ysl;
  short      ws, hs, wc, hc, w, h;

  sw = er->data0;
  contents = sw->contents;
  clipw = sw->clipw;
  xsl = sw->xsl;
  ysl = sw->ysl;
  clipw->x = er->x;
  clipw->y = er->y;
  ws = clipw->w = er->w;
  hs = clipw->h = er->h;
  wc = contents->w;  hc = contents->h;
  if ( wc <= ws ) {
    if ( hc <= hs ) {
      sw->x = sw->y = 0.0;
    }
    else {
      if ( wc > ws-10 )
        goto both;
      clipw->w = ws -= 10;
      sw->x = 0.0;
      ysl->h = hs;
    }
  }
  else {
    if ( hc <= hs-10 ) {
      clipw->h = hs -= 10;
      xsl->w = ws;
      sw->y = 0.0;
    }
    else {
both:
      clipw->w = xsl->w = ws -= 10;
      clipw->h = ysl->h = hs -= 10;
    }
  }
  sw->xslon = ws < contents->w;
  sw->yslon = hs < contents->h;
  contents->x = er->x;
  contents->y = er->y;
  if ( sw->yslon ) {
    h = contents->h-hs;
    contents->y -= (short)((float)h*sw->y);
  }
  if ( sw->xslon ) {
    w = contents->w-ws;
    contents->x -= (short)((float)w*sw->x);
  }
  if ( sw->yslon ) {
    clipw->x++;  clipw->w--;  clipw->y++;
    if ( sw->xslon )
      clipw->h--;
    else
      clipw->h -= 2;
  }
  else {
    if ( sw->xslon ) {
      clipw->y++;  clipw->h--;  clipw->x++;  clipw->w -= 2;
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
    xgeSetForeground ( xgec_White );
    xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
    if ( sw->xslon ) {
      xge_SetClipping ( sw->xsl );
      sw->xsl->redraw ( sw->xsl, false );
    }
    if ( sw->yslon ) {
      xge_SetClipping ( sw->ysl );
      sw->ysl->redraw ( sw->ysl, false );
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
case xgemsg_SLIDEBAR_COMMAND:
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
    sw->xsl->y = er->y+er->h-10;
    sw->ysl->x = er->x+er->w-10;
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
      if ( xge_PointInRect ( sw->xsl, x, y ) )
        return sw->xsl->msgproc ( sw->xsl, msg, key, x, y );
    }
    if ( sw->yslon ) {
      if ( xge_PointInRect ( sw->ysl, x, y ) )
        return sw->ysl->msgproc ( sw->ysl, msg, key, x, y );
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
  xge_widget *nsw, *clipw, *xsl, *ysl;

  nsw = xge_NewWidget ( window_num, prev, id, w, h, x, y, sw, contents,
                        xge_ScrollWidgetMsg, xge_DrawScrollWidget );
  sw->er = nsw;
  if ( nsw ) {
    clipw = sw->clipw = xge_NewEmptyWidget ( window_num, NULL, 0, w, h, x, y );
    xsl = sw->xsl = xge_NewSlidebarf ( window_num, clipw, 0, w, 10, x, y+h-10,
                                       &sw->x );
    ysl = sw->ysl = xge_NewVSlidebarf ( window_num, xsl, 1, 10, h, x+w-10, y,
                                        &sw->y );
    if ( clipw && xsl && ysl ) {
      xsl->up = ysl->up = clipw->up = nsw;
      contents->up = clipw;
      sw->x = sw->y = 0.0;
      sw->xslon = sw->yslon = false;
      sw->contents = contents;
      xge_SetupScrollWidgetPos ( nsw );
    }
    else
      return NULL;
  }
  return nsw;
} /*xge_NewScrollWidget*/

