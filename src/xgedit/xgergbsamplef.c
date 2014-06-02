
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


void xge_DrawRGBSamplef ( xge_widget *er, boolean onscreen )
{
  float *colour;

  colour = (float*)er->data0;
  xgeSetForeground ( xge_PixelColourf ( colour[0], colour[1], colour[2] ) );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawRGBSamplef*/

xge_widget *xge_NewRGBSamplef ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                float *data )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, data, NULL,
                         xge_EmptyMsg, xge_DrawRGBSamplef );
} /*xge_NewRGBSamplef*/

