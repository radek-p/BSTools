
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


void xge_DrawButton ( xge_widget *er, boolean onscreen )
{
  char *title;

  title = er->data0;
  switch ( er->state ) {
case xgestate_BUTTON_DEFAULT:
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_BUTTON_DEFAULT_A, xgec_BUTTON_DEFAULT_B, er->h-2 );
    break;
case xgestate_BUTTON_COMBO_0:
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_BUTTON_COMBO_0A, xgec_BUTTON_COMBO_0B, er->h-2 );
    break;
case xgestate_BUTTON_COMBO_1:
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_BUTTON_COMBO_1A, xgec_BUTTON_COMBO_1B, er->h-2 );
    break;
case xgestate_BUTTON_INACTIVE:
    xge_DrawVShadedRect ( er->w-2, er->h-2, er->x+1, er->y+1,
                          xgec_BUTTON_INACTIVE_A, xgec_BUTTON_INACTIVE_B, er->h-2 );
    break;
default:
    break;
  }
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  xgeSetForeground ( xgec_WIDGET_FOREGROUND );
  if ( title )
    xgeDrawString ( title, er->x+2, er->y+er->h-5 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawButton*/

boolean xge_ButtonMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er->state != xgestate_BUTTON_INACTIVE &&
       ((msg == xgemsg_MCLICK && key & xgemouse_LBUTTON_DOWN) ||
        (msg == xgemsg_KEY && key == 13)) ) {
    xge_callback ( er, xgemsg_BUTTON_COMMAND, 0, x, y );
    return true;
  }
  else
    return false;
} /*xge_ButtonMsg*/

xge_widget *xge_NewButton ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char *title )
{
  xge_widget *wdg;

  wdg = xge_NewWidget ( window_num, prev, id, w, h, x, y, title, NULL,
                        xge_ButtonMsg, xge_DrawButton );
  if ( wdg ) wdg->state = xgestate_BUTTON_DEFAULT;
  return wdg;
} /*xge_NewButton*/

