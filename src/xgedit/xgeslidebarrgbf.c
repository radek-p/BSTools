
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2014                            */
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


#define STEPS 32
void xge_DrawSlidebarfRGB ( xge_widget *er, boolean onscreen )
{
  int           x, x1, i;
  float         *slipos, a0, a1;
  xgecolour_int c0, c1;

  slipos = er->data0;
  x = er->x+1;
  if ( er->state == xgestate_MOVINGSLIDE ) {
    a0 = 0.4;  a1 = 0.8;
  }
  else {
    a0 = 1.0;  a1 = 0.6;
  }
  switch ( er->id & 0x03 ) {
case 0:
case 3:
    for ( i = 1; i <= STEPS; i++ ) {
      c0 = xge_PixelColourf ( a0*(float)i/(float)STEPS, 0.0, 0.0 );
      c1 = xge_PixelColourf ( a1*(float)i/(float)STEPS, 0.0, 0.0 );
      x1 = er->x + 1 + ((er->w-2)*i)/STEPS;
      xge_DrawVShadedRect ( x1-x, er->h-2, x, er->y+1, c0, c1, er->h-2 );
      x = x1;
    }
    break;
case 1:
    for ( i = 1; i <= STEPS; i++ ) {
      c0 = xge_PixelColourf ( 0.0, a0*(float)i/(float)STEPS, 0.0 );
      c1 = xge_PixelColourf ( 0.0, a1*(float)i/(float)STEPS, 0.0 );
      x1 = er->x + 1 + ((er->w-2)*i)/STEPS;
      xge_DrawVShadedRect ( x1-x, er->h-2, x, er->y+1, c0, c1, er->h-2 );
      x = x1;
    }
    break;
case 2:
    for ( i = 1; i <= STEPS; i++ ) {
      c0 = xge_PixelColourf ( 0.0, 0.0, a0*(float)i/(float)STEPS );
      c1 = xge_PixelColourf ( 0.0, 0.0, a1*(float)i/(float)STEPS );
      x1 = er->x + 1 + ((er->w-2)*i)/STEPS;
      xge_DrawVShadedRect ( x1-x, er->h-2, x, er->y+1, c0, c1, er->h-2 );
      x = x1;
    }
    break;
  }
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  x = er->x + 2 + (int)((*slipos)*(float)(er->w - 10));
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  xgeFillRectangle ( 6, 6, x, er->y+2 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawSlidebarfRGB*/

xge_widget *xge_NewSlidebarfRGB ( char window_num, xge_widget *prev, int id,
                                  short w, short h, short x, short y,
                                  float *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, NULL,
                         xge_SlidebarfMsg, xge_DrawSlidebarfRGB );
} /*xge_NewSlidebarfRGB*/

