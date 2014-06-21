
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


void xge_DrawVSlidebarf ( xge_widget *er, boolean onscreen )
{
  int   y;
  float *slipos;

  slipos = er->data0;
  if ( er->state == xgestate_MOVINGSLIDE )
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SLIDEBAR_A, xgec_SLIDEBAR_B, er->h-2 );
  else
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_SLIDEBAR_A, xgec_SLIDEBAR_D, er->h-2 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  y = er->y + 2 + (int)((*slipos)*(float)(er->h - 10));
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  xgeFillRectangle ( 6, 6, er->x+2, y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawVSlidebarf*/

boolean xge_VSlidebarfMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  float   z;
  float   *slipos;
  boolean inval;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  slipos = er->data0;
  inval = false;
  if ( er->state == xgestate_NOTHING ) {
    if ( msg == xgemsg_MCLICK ) {
      if ( (key & xgemouse_LBUTTON_DOWN) &&
           (key & xgemouse_LBUTTON_CHANGE) ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        er->state = xgestate_MOVINGSLIDE;
        xge_GrabFocus ( er, true );
        inval = true;
        goto update1;
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        z = *slipos - 1.0/(float)(er->h-10);
        z = max ( z, 0.0 );
        goto update2;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        z = *slipos + 1.0/(float)(er->h-10);
        z = min ( z, 1.0 );
        goto update2;
      }
    }
  }
  else if ( er->state == xgestate_MOVINGSLIDE ) {
    if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        if ( y != xge_yy ) {
          inval = false;
update1:
          xge_yy = y;
          z = (float)(y-er->y-5)/(float)(er->h-10);
update2:
          if ( *slipos != z ) {
            inval = true;
            *slipos = z;
                      /* float slidebar position is not passed as a parameter */
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
    }
  }
  return true;
} /*xge_VSlidebarfMsg*/

xge_widget *xge_NewVSlidebarf ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                float *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, NULL,
                         xge_VSlidebarfMsg, xge_DrawVSlidebarf );
} /*xge_NewVSlidebarf*/

