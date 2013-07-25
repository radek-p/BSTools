
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


void xge_DrawSwitch ( xge_widget *er, boolean onscreen )
{
  char    *title;
  boolean state, *switchvar;

  title     = er->data0;
  switchvar = er->data1;
  state = *switchvar;
  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );
  if ( state )
    xge_DrawVShadedRect ( er->h-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue5, xgec_Blue3, er->h-2 );
  else
    xge_DrawVShadedRect ( er->h-2, er->h-2, er->x+1, er->y+1,
                          xgec_Blue3, xgec_Blue5, er->h-2 );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->h-1, er->h-1, er->x, er->y );
  if ( state )
    xgeFillRectangle ( er->h-10, er->h-10, er->x+5, er->y+5 );
  if ( title )
    xgeDrawString ( title, er->x+er->h+2, er->y+er->h-4 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawSwitch*/

boolean xge_SwitchMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  boolean *switchvar, initval;

  if ( (msg == xgemsg_MCLICK && key & xgemouse_LBUTTON_DOWN) ||
       (msg == xgemsg_KEY && key == 13) ) {
    switchvar = er->data1;
    initval = *switchvar;
    *switchvar = (boolean)(!initval);
    xge_callback ( er, xgemsg_SWITCH_COMMAND, *switchvar, x, y );
    if ( *switchvar != initval ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;
  }
  else
    return false;
} /*xge_SwitchMsg*/

xge_widget *xge_NewSwitch ( char window_num, xge_widget *prev, int id,
                            short w, short h, short x, short y,
                            char *title, boolean *sw )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, title, sw,
                         xge_SwitchMsg, xge_DrawSwitch );
} /*xge_NewSwitch*/

