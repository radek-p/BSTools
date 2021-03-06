
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


void xge_DrawText ( xge_widget *er, boolean onscreen )
{
  char    *title;

  title = er->data0;
  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  xgeDrawString ( title, er->x, er->y+er->h-4 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawText*/

boolean xge_TextWidgetMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( msg == xgemsg_MCLICK &&
       (key & xgemouse_LBUTTON_DOWN) && (key & xgemouse_LBUTTON_CHANGE) ) {
    x = (x-er->x+xge_CHAR_WIDTH-1) / xge_CHAR_WIDTH;
    if ( x <= strlen ( er->data0 ) ) {
      if ( xge_callback ( er, xgemsg_TEXT_WIDGET_CLICK, key, x, y ) ) {
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      return true;
    }
  }
  return false;
} /*xge_TextWidgetMsg*/

xge_widget *xge_NewTextWidget ( char window_num, xge_widget *prev, int id,
                                short w, short h, short x, short y,
                                char *text )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, text, NULL,
                         xge_TextWidgetMsg, xge_DrawText );
} /*xge_NewTextWidget*/

