
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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
static boolean _xge_HScrollbarAuxResize ( xge_hscrollbar *hs, short newwp, short newww )
{
  short oldwp, oldww, oldppos, newwsl, bw, hpos;

  bw = hs->er->w-2;
  oldwp = hs->wpage;
  oldww = hs->wwin;
  if ( newwp != oldwp || newww != oldww ) {
    hs->wpage = newwp;
    hs->wwin = newww;
    oldppos = hs->pagehpos;
    hs->wsl = newwsl = ((int)bw*newww)/(int)newwp;
    hpos = oldwp > oldww ? ((int)oldppos*(newwp-newww))/(oldwp-oldww) : 0;
    if ( hpos+newww > newwp )
      hpos = newwp-newww;
    hs->pagehpos = hpos = max ( hpos, 0 );
    if ( newwp > newww )
      hs->slpos = ((int)(bw-newwsl)*hpos)/(int)(newwp-newww);
    else
      hs->slpos = 0;
  }
  return hs->wpage > hs->wwin;
} /*_xge_HScrollbarAuxResize*/

boolean xge_HScrollBarResizeWindow ( xge_widget *er, short newwidth )
{
  xge_hscrollbar *hs;

  hs = er->data0;
  return _xge_HScrollbarAuxResize ( hs, hs->wpage, newwidth );
} /*xge_HScrollBarResizeWindow*/

boolean xge_HScrollBarResizePage ( xge_widget *er, short newwidth )
{
  xge_hscrollbar *hs;

  hs = er->data0;
  return _xge_HScrollbarAuxResize ( hs, newwidth, hs->wwin );
} /*xge_HScrollBarResizePage*/

short xge_HScrollBarSetPagePos ( xge_widget *er, short newpos, boolean redraw )
{
  xge_hscrollbar *hs;

  hs = er->data0;
  if ( hs->wwin + newpos > hs->wpage )
    newpos = hs->wpage - hs->wwin;
  if ( newpos < 0 )
    newpos = 0;
  hs->pagehpos = newpos;
  hs->slpos = ((int)(er->w-hs->wsl)*newpos)/(hs->wpage-hs->wwin);
  if ( redraw ) {
    xge_SetClipping ( er );
    er->redraw ( er, true );
  }
  return newpos;
} /*xge_HScrollBarSetPagePos*/

void xge_DrawHScrollBar ( xge_widget *er, boolean onscreen )
{
  xge_hscrollbar *hs;
  int            x;

  hs = er->data0;
  if ( er->state == xgestate_SCROLLING ) {
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SCROLLBAR_A, xgec_SCROLLBAR_B, er->h-2 );
    xgeSetForeground ( xgec_MENU_DARKBKG );
  }
  else {
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SCROLLBAR_C, xgec_SCROLLBAR_D, er->h-2 );
    xgeSetForeground ( xgec_INFOMSG_BACKGROUND );
  }
  x = er->x+hs->slpos+1;
  xgeFillRectangle ( hs->wsl-2, er->h-3, x+1, er->y+1 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  xgeDrawRectangle ( hs->wsl-1, er->h-3, x, er->y+1 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawHScrollBar*/

boolean xge_HScrollBarMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_hscrollbar *hs;
  xge_widget     *window;
  boolean        inval;
  int            slpos;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  hs = er->data0;
  window = er->data1;
  inval = false;
  switch ( er->state ) {
case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        hs->slpoint = x - (er->x+hs->slpos);
        if ( hs->slpoint < 3 )
          hs->slpoint = 3;
        else if ( hs->slpoint >= hs->wsl+1 )
          hs->slpoint = hs->wsl;
        slpos = x - er->x - hs->slpoint;
        er->state = xgestate_SCROLLING;
        xge_GrabFocus ( er, true );
        inval = true;
        goto update1;
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        slpos = hs->slpos - 1;
        goto update2;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        slpos = hs->slpos + 1;
        goto update2;
      }

  default:
      return false;
    }

case xgestate_SCROLLING:
    switch ( msg ) {
  case xgemsg_MMOVE:
  case xgemsg_MCLICK:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        slpos = x - (er->x+hs->slpoint);
update1:
        xge_xx = x;
update2:
        if ( slpos < 0 )
          slpos = 0;
        else if ( slpos > er->w-2-hs->wsl )
          slpos = er->w-2-hs->wsl;
        if ( slpos != hs->slpos ) {
          inval = true;
          hs->slpos = slpos;
          hs->pagehpos = ((int)slpos*(hs->wpage-hs->wwin))/(er->w-2-hs->wsl);
/*printf ( "slpos = %3d, pagepos = %3d\n", hs->slpos, hs->pagehpos );*/
          if ( window )
            window->msgproc ( window, xgemsg_HSCROLLBAR_COMMAND, er->id, x, y );
          else
            xge_callback ( er, xgemsg_HSCROLLBAR_COMMAND, er->id, x, y );
        }
        if ( inval ) {
          xge_SetClipping ( er );
          er->redraw ( er, true );
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

default:
    return false;
  }
} /*xge_HScrollBarMsg*/

xge_widget *xge_NewHScrollBar ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                short wpage, short wwin,
                                xge_hscrollbar *hs, xge_widget *window )
{
  xge_widget *hsb;

  if ( (hsb = xge_NewWidget ( window_num, prev, id, w, h, x, y, hs, window,
                              xge_HScrollBarMsg, xge_DrawHScrollBar )) ) {
    hs->wpage = wpage;
    hs->wwin = wwin;
    hs->wsl = wpage > wwin ? (int)(w-2)*wwin/wpage : w-2;
    hs->slpos = 0;
    hs->slpoint = 0;
    hs->pagehpos = 0;
  }
  hs->er = hsb;
  return hsb;
} /*xge_NewHScrollBar*/

