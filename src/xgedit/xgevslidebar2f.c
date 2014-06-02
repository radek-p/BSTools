
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


void xge_DrawVSlidebar2f ( xge_widget *er, boolean onscreen )
{
  int    y0, y1;
  float *slipos0, *slipos1;

  slipos0 = (float*)er->data0;
  slipos1 = &((float*)er->data0)[1];
  if ((er->state == xgestate_MOVINGSLIDE2A ||
       er->state == xgestate_MOVINGSLIDE2B ||
       er->state == xgestate_MOVINGSLIDE ) )
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue7, xgec_Blue4, er->h-2 );
  else
    xge_DrawHShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue3, xgec_Blue5, er->h-2 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  y0 = er->y + 2 + (int)((*slipos0)*(float)(er->h - 10));
  y1 = er->y + 8 + (int)((*slipos1)*(float)(er->h - 10));
  xgeSetForeground ( xgec_White );
  xgeFillRectangle ( 6, y1-y0, er->x+2, y0 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawVSlidebar2f*/

boolean xge_VSlidebar2fMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  float  z, z0, z1, dz, mz;
  float  *slipos0, *slipos1;
  boolean inval;
  int     y0, y1;

  if ( msg == xgemsg_SPECIAL_KEY )
    return false;

  slipos0 = (float*)er->data0;
  slipos1 = &((float*)er->data0)[1];
  dz = *slipos1-*slipos0;
  mz = 0.5*(*slipos0+*slipos1);
  y0 = er->y + 5 + (int)((*slipos0)*(float)(er->h - 10));
  y1 = er->y + 5 + (int)((*slipos1)*(float)(er->h - 10));
  inval = false;
  switch ( er->state ) {
case xgestate_NOTHING:
    if ( msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        xge_GrabFocus ( er, true );
        inval = true;
        if ( 3*y <= 2*y0+y1 ) {
          er->state = xgestate_MOVINGSLIDE2A;
          goto update3;
        }
        else if ( 3*y >= y0+2*y1 ) {
          er->state = xgestate_MOVINGSLIDE2B;
          goto update4;
        }
        else {
          er->state = xgestate_MOVINGSLIDE;
          goto update1;
        }
      }
      else if ( (key & xgemouse_WHEELFW_DOWN) &&
                (key & xgemouse_WHEELFW_CHANGE) ) {
        z0 = mz - 1.0/(float)(er->h-10) - 0.5*dz;
        z0 = max ( z0, 0.0 );
        z1 = z0 + dz;
        z = 0.5*(z0+z1);
        goto update2;
      }
      else if ( (key & xgemouse_WHEELBK_DOWN) &&
                (key & xgemouse_WHEELBK_CHANGE) ) {
        z1 = mz + 1.0/(float)(er->h-10) + 0.5*dz;
        z1 = min ( z1, 1.0 );
        z0 = z1 - dz;
        z = 0.5*(z0+z1);
        goto update2;
      }
    }
    break;

case xgestate_MOVINGSLIDE:
    if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        if ( y != xge_yy ) {
          inval = false;
update1:
          xge_yy = y;
          z = (float)(y-er->y-5)/(float)(er->h-10);
          z0 = z - 0.5*dz;
          z1 = z + 0.5*dz;
update2:
          if ( mz != z ) {
            inval = true;
            *slipos0 = max ( z0, 0.0 );
            *slipos1 = min ( z1, 1.0 );
                    /* float slidebar position is not passed as a parameter */
            if ( !er->up ||
                 !er->up->msgproc ( er->up, xgemsg_SLIDEBAR2_COMMAND, er->id, x, y ) )
              xge_callback ( er, xgemsg_SLIDEBAR2_COMMAND, 0, x, y );
          }
          if ( inval ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else
        goto exitmode;
    }
    break;

case xgestate_MOVINGSLIDE2A:
    if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        if ( y != xge_yy ) {
          inval = false;
update3:
          xge_yy = y;
          z = (float)(y-er->y-5)/(float)(er->h-10);
          if ( *slipos0 != z ) {
            inval = true;
            *slipos0 = z;
            *slipos1 = max ( z, *slipos1 );
                    /* float slidebar position is not passed as a parameter */
            if ( !er->up ||
                 !er->up->msgproc ( er->up, xgemsg_SLIDEBAR2_COMMAND, er->id, x, y ) )
              xge_callback ( er, xgemsg_SLIDEBAR2_COMMAND, 0, x, y );
          }
          if ( inval ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else
        goto exitmode;
    }
    break;

case xgestate_MOVINGSLIDE2B:
    if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( y < er->y+5 ) y = (short)(er->y+5);
        else if ( y > er->y+er->h-5 ) y = (short)(er->y+er->h-5);
        if ( y != xge_yy ) {
          inval = false;
update4:
          xge_yy = y;
          z = (float)(y-er->y-5)/(float)(er->h-10);
          if ( *slipos1 != z ) {
            inval = true;
            *slipos1 = z;
            *slipos0 = min ( z, *slipos0 );
                    /* float slidebar position is not passed as a parameter */
            if ( !er->up ||
                 !er->up->msgproc ( er->up, xgemsg_SLIDEBAR2_COMMAND, er->id, x, y ) )
              xge_callback ( er, xgemsg_SLIDEBAR2_COMMAND, 0, x, y );
          }
          if ( inval ) {
            xge_SetClipping ( er );
            er->redraw ( er, true );
          }
        }
      }
      else {
exitmode:
        er->state = xgestate_NOTHING;
        xge_ReleaseFocus ( er );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
    }
    break;

default:
    return false;
  }
  return true;
} /*xge_VSlidebar2fMsg*/

xge_widget *xge_NewVSlidebar2f ( char window_num, xge_widget *prev, int id,
                                 short w, short h, short x, short y,
                                 float *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, NULL,
                         xge_VSlidebar2fMsg, xge_DrawVSlidebar2f );
} /*xge_NewVSlidebar2f*/

