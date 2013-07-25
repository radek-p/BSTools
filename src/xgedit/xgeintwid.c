
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
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


void xge_DrawIntWidget ( xge_widget *er, boolean onscreen )
{
  int            nd, mm, w, val, x, l;
  xge_int_widget *iw;
  char           s[13];

  val = *(int*)(er->data1);
  iw   = er->data0;
  mm = abs(iw->minvalue);  nd = abs(iw->maxvalue);
  mm = max ( mm, nd );
  nd   = 1;
  while ( mm > 9 ) { nd++;  mm /= 10; }

  w = 6*nd+4;
  x = er->x+er->w-w;
  xge_DrawVShadedRect ( w-1, er->h-2, x, er->y+1,
                        xgec_Green3, xgec_Green5, er->h-2 );
  if ( er->w-w-2 > 0 ) {
    xge_DrawVShadedRect ( er->w-w-2, er->h-2, er->x+1, er->y+1,
                          xgec_Green5, xgec_Green7, er->h-2 );
  }
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  xgeDrawLine ( x-1, er->y+er->h-1, x-1, er->y+1 );
  sprintf ( s, "%d", val );
  l = strlen(s);
  xgeDrawString ( s, er->x+er->w-2-6*l, er->y+er->h-5 );
  if ( iw->title )
    xgeDrawString ( iw->title, er->x+2, er->y+er->h-5 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawIntWidget*/

boolean xge_IntWidgetMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_int_widget *iw;
  int            value;

  iw = (xge_int_widget*)er->data0;
  value = *(int*)(er->data1);
  switch ( msg ) {
case xgemsg_MCLICK:
    if ( (key & xgemouse_LBUTTON_DOWN) &&
         (key & xgemouse_LBUTTON_CHANGE) )
      goto increase;
    else if ( (key & xgemouse_WHEELBK_DOWN) &&
              (key & xgemouse_WHEELBK_CHANGE) )
      goto increase;
    else if ( (key & xgemouse_RBUTTON_DOWN) &&
              (key & xgemouse_RBUTTON_CHANGE) )
      goto decrease;
    else if ( (key & xgemouse_WHEELFW_DOWN) &&
              (key & xgemouse_WHEELFW_CHANGE) )
      goto decrease;
    break;

case xgemsg_KEY:
    switch ( key ) {
  case '+':
increase:
      if ( value < iw->maxvalue ) {
        value ++;
        goto processed;
      }
      break;
  case '-':
decrease:
      if ( value > iw->minvalue ) {
        value --;
        goto processed;
      }
      break;
  default:
      break;
    }
    break;

default:
    break;
  }
  return false;

processed:
  xge_callback ( er, xgemsg_INT_WIDGET_COMMAND, value, x, y );
  if ( *(int*)(er->data1) == value )
  xge_SetClipping ( er );
  er->redraw ( er, true );
  return true;
} /*xge_IntWidgetMsg*/

xge_widget *xge_NewIntWidget ( char window_num, xge_widget *prev, int id,
                               short w, short h, short x, short y,
                               int minvalue, int maxvalue,
                               xge_int_widget *iw, char *title, int *valptr )
{
  iw->minvalue = minvalue;
  iw->maxvalue = maxvalue;
  iw->title = title;
  return (iw->er = xge_NewWidget (
                         window_num, prev, id, w, h, x, y, iw, valptr,
                         xge_IntWidgetMsg, xge_DrawIntWidget ) );
} /*xge_NewIntWidget*/

