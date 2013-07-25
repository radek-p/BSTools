
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

void xge_AddPopup ( xge_widget *er )
{
    /* add a popup widget to the end of the list */
  xge_SetWindow ( er->window_num );
  er->next = NULL;
  er->prev = xge_windesc[xge_current_win].popup1;
  xge_windesc[xge_current_win].popup1 = er;
  if ( !xge_windesc[xge_current_win].popup0 )
    xge_windesc[xge_current_win].popup0 = er;
  xge_SetClipping ( er );
  er->redraw ( er, true );
} /*xge_AddPopup*/

/* the application must keep its own pointers to the widgets */
/* displayed as popups. */

void xge_RemovePopup ( boolean redraw )
{
  xge_widget *er;
  int        win;

      /* remove the last popup widget from the list */
  win = xge_current_win;
  if ( (er = xge_windesc[win].popup1) ) {
    if ( xge_lastwin ) {
      xge_lastwin->msgproc ( xge_lastwin, xgemsg_EXITING, 0, 0, 0 );
      xge_lastwin = NULL;
    }
    xge_callback ( NULL, xgemsg_POPUP_REMOVED, er->id, 0, 0 );
    if ( win != xge_current_win )
      xge_SetWindow ( win );
    xge_ReleaseFocus ( er );
    er = er->prev;
    xge_windesc[xge_current_win].popup1 = er;
    if ( er )
      er->next = NULL;
    else
      xge_windesc[xge_current_win].popup0 = NULL;
      /* the following procedure call sends the ENTERING message to   */
      /* whatever widget is pointed by the cursor; it is assumed that */
      /* the application took care of the focus */
    xge_dispatch_message ( xgemsg_MMOVE, xge_mouse_buttons,
                           (short)xge_mouse_x, (short)xge_mouse_y );
    if ( redraw )
      xge_Redraw ();
  }
} /*xge_RemovePopup*/

void xge_RemovePopups ( boolean redraw )
{
  int win;

      /* remove all popups of the window */
  win = xge_current_win;
  if ( xge_windesc[win].popup1 ) {
    xge_callback ( NULL, xgemsg_POPUPS_REMOVED, win, 0, 0 );
    if ( win != xge_current_win )
      xge_SetWindow ( win );
    while ( xge_windesc[win].popup0 ) {
      xge_RemovePopup ( false );
      if ( xge_current_win != win )
        xge_SetWindow ( win );
    }
    if ( redraw )
      xge_Redraw ();
  }
} /*xge_RemovePopups*/

boolean xge_IsPopupOn ( xge_widget *er )
{
  int win;
  xge_widget *pp;

  win = er->window_num;
  for ( pp = xge_windesc[win].popup1; pp; pp = pp->prev )
    if ( pp == er )
      return true;
  return false;
} /*xge_IsPopupOn*/

