
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


/* angle increments for a wheel impulse */
#define STEP0 (10.0*PI/180.0)
#define STEP1 (PI/180.0)
#define STEP2 (0.1*PI/180.0)

void xge_DrawDiald ( xge_widget *er, boolean onscreen )
{
  short  xc, yc, x, y, d, r;
  double *dialpos;
  char   *title;

  dialpos = er->data0;
  title   = er->data1;
  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->x+er->w-1, er->y+er->h-1, er->x, er->y );
  d = (short)min ( er->w, er->h );
  if ( !(d & 0x0001) )
    d --;
  r = (short)(d/2);
  xc = (short)(er->x + r - 1);
  yc = (short)(er->y + r - 1);
  if ( er->state == xgestate_TURNINGDIAL )
    xgeSetForeground ( xgec_Blue6 );
  else
    xgeSetForeground ( xgec_Blue3 );
  xgeFillArc ( d-1, d-1, er->x, er->y, 0, 360*64 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawArc ( d-1, d-1, er->x, er->y, 0, 360*64 );
  x = (short)(xc + (r-5)*cos(*dialpos));
  y = (short)(yc - (r-5)*sin(*dialpos));
  xgeSetForeground ( xgec_White );
  xgeFillArc ( 7, 7, x-2, y-2, 0, 360*64 );
  if ( title ) {
    if ( er->w > er->h )  /* title aside */
      { x = (short)(er->x+d+2);  y = (short)(er->y+r+4); }
    else                  /* title below */
      { x = (short)(er->x+r-6*strlen(title)/2);  y = (short)(er->y+d+11); }
    xgeDrawString ( title, x, y );
  }

  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawDiald*/

boolean xge_DialdMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  double  z;
  double  *dialpos;
  boolean inval;
  short   xc, yc, d, r;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  dialpos = er->data0;
  inval = false;
  if ( er->state == xgestate_NOTHING ) {
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( key & xgemouse_WHEELFW_CHANGE ) {
        if ( key & xgemouse_LBUTTON_DOWN )      z = *dialpos + STEP0;
        else if ( key & xgemouse_RBUTTON_DOWN ) z = *dialpos + STEP2;
        else                                    z = *dialpos + STEP1;
        while ( z > PI )
          z -= 2.0*PI;
        goto update2;
      }
      else if ( key & xgemouse_WHEELBK_CHANGE) {
        if ( key & xgemouse_LBUTTON_DOWN )      z = *dialpos - STEP0;
        else if ( key & xgemouse_RBUTTON_DOWN ) z = *dialpos - STEP2;
        else                                    z = *dialpos - STEP1;
        while ( z <= -PI )
          z += 2.0*PI;
        goto update2;
      }
      else if ( (key & xgemouse_LBUTTON_DOWN) &&
                (key & xgemouse_LBUTTON_CHANGE) ) {
        er->state = xgestate_TURNINGDIAL;
        xge_GrabFocus ( er, true );
        inval = true;
        goto update1;
      }
      break;
  default:
      return false;
    }
  }
  else if ( er->state == xgestate_TURNINGDIAL ) {
    if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( key & xgemouse_WHEELFW_CHANGE ) {
          z = *dialpos + STEP0;
          while ( z > PI )
            z -= 2.0*PI;
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          goto update2;
        }
        else if ( key & xgemouse_WHEELBK_CHANGE ) {
          z = *dialpos - STEP0;
          while ( z <= -PI )
            z += 2.0*PI;
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          goto update2;
        }
        else if ( x != xge_xx || y != xge_yy ) {
update1:
          xge_xx = x;
          xge_yy = y;
          d = (short)min ( er->w, er->h );
          if ( !(d & 0x0001) )
            d --;
          r = (short)(d/2);
          xc = (short)(er->x + r - 1);
          yc = (short)(er->y + r - 1);
          if ( xc != x || yc != y ) {
            z = atan2 ( yc-y, x-xc );
update2:
            if ( *dialpos != z ) {
              inval = true;
              *dialpos = z;
                      /* double dial position is not passed as a parameter */
              xge_callback ( er, xgemsg_DIAL_COMMAND, 0, x, y );
            }
            if ( inval ) {
              xge_SetClipping ( er );
              er->redraw ( er, true );
            }
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
} /*xge_DialdMsg*/

xge_widget *xge_NewDiald ( char window_num, xge_widget *prev, int id,
                           short w, short h, short x, short y,
                           char *title, double *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, title,
                         xge_DialdMsg, xge_DrawDiald );
} /*xge_NewDiald*/

