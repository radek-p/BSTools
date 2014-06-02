
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2010                            */
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


void xge_DrawMenu ( xge_widget *er, boolean onscreen )
{
  xge_widget *widget;

  xgeSetForeground ( xgec_MENU_BACKGROUND );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );

  for ( widget = er->data0;  widget;  widget = widget->next ) {
    xge_SetClipping ( widget );
    widget->redraw ( widget, false );
  }

  xge_SetClipping ( er );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawMenu*/

boolean xge_MenuMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_widget *widget;

  widget = er->data1;

  switch ( msg ) {
case xgemsg_ENTERING:
    xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    er->data2 = NULL;
    return true;

case xgemsg_EXITING:
    xge_CallMsgProc ( (xge_widget**)&er->data2, NULL, xgemsg_NULL, 0, 0, 0 );
    return true;

case xgemsg_RESIZE:
    er->w = x;
    er->h = y;
    if ( widget )
      xge_RepositionWidgets ( x, y, er->x, er->y, widget );
    if ( key ) {
      xge_SetClipping ( er );
      er->redraw ( er, true );
    }
    return true;

case xgemsg_SLIDEBAR_COMMAND:
    return false;

default:
    for ( ;  widget;  widget = widget->prev )
      if ( xge_PointInRect ( widget, x, y ) ) {
        if ( xge_CallMsgProc ( (xge_widget**)&er->data2, widget, msg, key, x, y ) )
          return true;
      }
        /* none found */
    xge_CallMsgProc ( (xge_widget**)&er->data2, NULL, xgemsg_NULL, 0, 0, 0 );
    return false;
  }
} /*xge_MenuMsg*/

boolean xge_PopupMenuMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( msg == xgemsg_MCLICK &&
       key & xgemouse_LBUTTON_DOWN && key & xgemouse_LBUTTON_CHANGE &&
       !xge_PointInRect ( er, x, y ) ) {
    xge_RemovePopup ( true );
    return true; 
  }
  else return xge_MenuMsg ( er, msg, key, x, y );
} /*xge_PopupMenuMsg*/

static void xge_MenuSetUp ( xge_widget *er )
{
  xge_widget *w;

  for ( w = er->data0; w; w = w->next )
    w->up = er;
} /*xge_MenuSetUp*/

xge_widget *xge_NewMenu ( char window_num, xge_widget *prev, int id,
                          short w, short h, short x, short y,
                          xge_widget *widgetlist )
{
  xge_widget *first, *er;

  first = widgetlist;
  if ( first )
    while ( first->prev )
      first = first->prev;
  er = xge_NewWidget ( window_num, prev, id, w, h, x, y, first, widgetlist,
                       xge_MenuMsg, xge_DrawMenu );
  xge_MenuSetUp ( er );
  xge_RepositionWidgets ( er->w, er->h, er->x, er->y, er->data1 );
  return er;
} /*xge_NewMenu*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_DrawFMenu ( xge_widget *er, boolean onscreen )
{
  xge_DrawMenu ( er, false );
  xgeSetForeground ( xgec_WIDGET_FRAME );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawFMenu*/

xge_widget *xge_NewFMenu ( char window_num, xge_widget *prev, int id,
                           short w, short h, short x, short y,
                           xge_widget *widgetlist )
{
  xge_widget *first, *er;

  first = widgetlist;
  if ( first )
    while ( first->prev )
      first = first->prev;
  er = xge_NewWidget ( window_num, prev, id, w, h, x, y, first, widgetlist,
                       xge_MenuMsg, xge_DrawFMenu );
  xge_MenuSetUp ( er );
  return er;
} /*xge_NewFMenu*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_SetMenuWidgets ( xge_widget *menu, xge_widget *widgetlist,
                          boolean redraw )
{
  xge_widget *first;

  first = widgetlist;
  if ( first )
    while ( first->prev )
      first = first->prev;
  menu->data0 = first;
  menu->data1 = widgetlist;
  while ( first ) {
    first->up = menu;
    first = first->next;
  }
  xge_RepositionWidgets ( menu->w, menu->h, menu->x, menu->y, menu->data1 );
  if ( redraw ) {
    xge_SetClipping ( menu );
    menu->redraw ( menu, true );
  }
} /*xge_SetMenuWidgets*/

