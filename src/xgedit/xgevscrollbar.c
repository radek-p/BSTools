
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
static boolean _xge_VScrollbarAuxResize ( xge_vscrollbar *vs, short newhp, short newhw )
{
  short oldhp, oldhw, oldppos, newhsl, bh, vpos;

  bh = vs->er->h-2;
  oldhp = vs->hpage;    vs->hpage = newhp;
  oldhw = vs->hwin;     vs->hwin = newhw;
  oldppos = vs->pagevpos;
  vs->hsl = newhsl = ((int)bh*newhw)/(int)newhp;
  vpos = oldhp > oldhw ? ((int)oldppos*(newhp-newhw))/(oldhp-oldhw) : 0;
  if ( vpos+newhw > newhp )
    vpos = newhp-newhw;
  vs->pagevpos = vpos = max ( vpos, 0 );
  if ( newhp > newhw )
    vs->slpos = ((int)(bh-newhsl)*vpos)/(int)(newhp-newhw);
  else
    vs->slpos = 0;
  return vs->hpage > vs->hwin;
} /*_xge_VScrollbarAuxResize*/

boolean xge_VScrollBarResizeWindow ( xge_widget *er, short newheight )
{
  xge_vscrollbar *vs;

  vs = er->data0;
  return _xge_VScrollbarAuxResize ( vs, vs->hpage, newheight );
} /*xge_VScrollBarResizeWindow*/

boolean xge_VScrollBarResizePage ( xge_widget *er, short newheight )
{
  xge_vscrollbar *vs;

  vs = er->data0;
  return _xge_VScrollbarAuxResize ( vs, newheight, vs->hwin );
} /*xge_VScrollBarResizePage*/

short xge_VScrollBarSetPagePos ( xge_widget *er, short newpos, boolean redraw )
{
  xge_vscrollbar *vs;

  vs = er->data0;
  if ( vs->hwin + newpos > vs->hpage )
    newpos = vs->hpage - vs->hwin;
  if ( newpos < 0 )
    newpos = 0;
  vs->pagevpos = newpos;
  vs->slpos = ((int)(er->h-vs->hsl)*newpos)/(vs->hpage-vs->hwin);
  if ( redraw ) {
    xge_SetClipping ( er );
    er->redraw ( er, true );
  }
  return newpos;
} /*xge_VScrollBarSetPagePos*/

void xge_DrawVScrollBar ( xge_widget *er, boolean onscreen )
{
  xge_vscrollbar *vs;
  int            y;

  vs = er->data0;
  if ( er->state == xgestate_SCROLLING )
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue7, xgec_Blue4, er->h-2 );
  else
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue3, xgec_Blue5, er->h-2 );
  y = er->y+vs->slpos+1;
  xgeSetForeground ( xgec_INFOMSG_BACKGROUND );
  xgeFillRectangle ( er->w-3, vs->hsl-2, er->x+1, y+1 );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  xgeDrawRectangle ( er->w-3, vs->hsl-1, er->x+1, y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawVScrollBar*/

boolean xge_VScrollBarMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_vscrollbar *vs;
  xge_widget     *window;
  boolean        inval;
  int            slpos;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  vs = er->data0;
  window = er->data1;
  inval = false;
  switch ( er->state ) {
case xgestate_NOTHING:
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( key & xgemouse_LBUTTON_DOWN ) {
        vs->slpoint = y - (er->y+vs->slpos);
        if ( vs->slpoint < 3 )
          vs->slpoint = 3;
        else if ( vs->slpoint >= vs->hsl+1 )
          vs->slpoint = vs->hsl;
        slpos = y - er->y - vs->slpoint;
        er->state = xgestate_SCROLLING;
        xge_GrabFocus ( er, true );
        inval = true;
        goto update1;
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        slpos = vs->slpos - 1;
        goto update2;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        slpos = vs->slpos + 1;
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
        slpos = y - (er->y+vs->slpoint);
update1:
        xge_yy = y;
update2:
        if ( slpos < 0 )
          slpos = 0;
        else if ( slpos > er->h-2-vs->hsl )
          slpos = er->h-2-vs->hsl;
        if ( slpos != vs->slpos ) {
          inval = true;
          vs->slpos = slpos;
          vs->pagevpos = ((int)slpos*(vs->hpage-vs->hwin))/(er->h-2-vs->hsl);
/*printf ( "slpos = %3d, pagepos = %3d\n", vs->slpos, vs->pagevpos );*/
          if ( window )
            window->msgproc ( window, xgemsg_VSCROLLBAR_COMMAND, er->id, x, y );
          else
            xge_callback ( er, xgemsg_VSCROLLBAR_COMMAND, er->id, x, y );
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
} /*xge_VScrollBarMsg*/

xge_widget *xge_NewVScrollBar ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                short hpage, short hwin,
                                xge_vscrollbar *vs, xge_widget *window )
{
  xge_widget *vsb;

  if ( (vsb = xge_NewWidget ( window_num, prev, id, w, h, x, y, vs, window,
                              xge_VScrollBarMsg, xge_DrawVScrollBar )) ) {
    vs->hpage = hpage;
    vs->hwin = hwin;
    vs->hsl = hpage > hwin ? (int)(h-2)*hwin/hpage : h-2;
    vs->slpos = 0;
    vs->slpoint = 0;
    vs->pagevpos = 0;
  }
  vs->er = vsb;
  return vsb;
} /*xge_NewVScrollBar*/

