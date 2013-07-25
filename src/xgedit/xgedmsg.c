
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

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget        *xge_errmsg_edr;
int               xge_usr_msg_key = -1;
static int        xge_InfoNLines, xge_InfoMaxLength;

/* ///////////////////////////////////////////////////////////////////////// */
static void xge_RedrawErrorMessage ( xge_widget *er, boolean onscreen )
{
  int sl;

  xge_SetClipping ( xge_errmsg_edr );
  xgeSetForeground ( xge_msgbkcolour );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  sl = strlen(xge_errmsg_msgtext);
  xgeDrawString (  xge_errmsg_msgtext, er->x+((er->w-6*sl)/2), er->y+16 );
  xgeDrawString ( "OK", er->x+(er->w/2)-7, er->y+er->h-6 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_RedrawErrorMessage*/

void xge_RemoveErrorMessage ( void )
{
  if ( xge_errmsg_edr->state == xgestate_MESSAGE ) {
    xge_SetWindow ( xge_errmsg_win );
    xge_errmsg_win = -1;
    xge_errmsg_edr->state = xgestate_NOTHING;
    xge_RemovePopup ( true );
    xge_callback ( NULL, xgemsg_USER_MESSAGE_DISMISSED, xge_usr_msg_key, 0, 0 );
    xge_usr_msg_key = -1;
  }
} /*xge_RemoveErrorMessage*/

static boolean xge_ErrorMessageMsg ( xge_widget *er, int msg, int key,
                                     short x, short y )
{
  Cursor c;

  if ( xge_current_win == xge_errmsg_win &&
       xge_errmsg_edr->state == xgestate_MESSAGE && !xge_notinfocus ) {
    switch ( msg ) {
  case xgemsg_MCLICK:
      if ( xge_PointInRect ( er, x, y ) &&
           key & xgemouse_LBUTTON_DOWN &&
           key & xgemouse_LBUTTON_CHANGE )
        xge_RemoveErrorMessage ();
      return true;
  case xgemsg_MMOVE:
      if ( er->window_num == xge_current_win ) {
        if ( xge_PointInRect ( er, x, y ) )
          c = xgeCURSOR_DEFAULT;
        else
          c = xgeCURSOR_CIRCLE;
        if ( c != xge_windesc[xge_current_win].cursor )
          xge_SetCurrentWindowCursor ( c );
      }
      return true;
  case xgemsg_KEY:
      if ( key == 0x0d ) {
        xge_RemoveErrorMessage ();
        return true;
      }
      break;
  default:
      break;
    }
  }
  return false;
} /*xge_ErrorMessageMsg*/

void _xge_DisplayErrorMessage ( char *message, xgecolour_int bk, int key )
{
  int win;

  win = xge_current_win;
  xge_RemoveErrorMessage ();
  xge_SetWindow ( win );

  xge_usr_msg_key = key;
  xge_errmsg_msgtext = message;
  xge_info_msgtext   = NULL;
  xge_msgbkcolour    = bk;
  xge_SetupEdRect ( xge_current_win, xge_errmsg_edr, 0, 0, xge_WIDTH-20, 60,
                    (short)((xge_current_width-xge_WIDTH+20)/2),
                    (short)((2*xge_current_height)/3-30),
                    &xge_ErrorMessageMsg, &xge_RedrawErrorMessage );
  xge_errmsg_edr->state = xgestate_MESSAGE;
  xge_AddPopup ( xge_errmsg_edr );
  xge_GrabFocus ( xge_errmsg_edr, true );
  if ( xge_PointInRect ( xge_errmsg_edr, xge_mouse_x, xge_mouse_y ) )
    xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
  else
    xge_SetCurrentWindowCursor ( xgeCURSOR_CIRCLE );
  xge_SetOtherWindowsCursor ( xgeCURSOR_CIRCLE );
  xge_errmsg_win = xge_current_win;
  if ( xgeevent.type != ConfigureNotify )
    xgeRaiseWindow ();
} /*_xge_DisplayErrorMessage*/

void xge_DisplayErrorMessage ( char *message, int key )
{
  _xge_DisplayErrorMessage ( message, xgec_ERRORMSG_BACKGROUND, key );
} /*xge_DisplayErrorMessage*/

void xge_DisplayWarningMessage ( char *message, int key )
{
  _xge_DisplayErrorMessage ( message, xgec_WARNINGMSG_BACKGROUND, key );
} /*xge_DisplayWarningMessage*/

static void xge_RedrawInfoMessage ( xge_widget *er, boolean onscreen )
{
  int  i;
  char *msgstr;

  xge_SetClipping ( xge_errmsg_edr );
  xgeSetForeground ( xgec_INFOMSG_BACKGROUND );
  xgeFillRectangle ( er->w, er->h, er->x, er->y );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  for ( i = 0; i < xge_InfoNLines; i++ ) {
    msgstr = xge_info_msgtext[i];
    if ( msgstr[0] )
      xgeDrawString ( msgstr,
        er->x+((er->w-6*xge_InfoMaxLength)/2), er->y+16*(i+1) );
  }

  xgeDrawString ( "OK", er->x+(er->w/2)-7, er->y+er->h-6 );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_RedrawInfoMessage*/

void xge_DisplayInfoMessage ( char **msglines, int key )
{
  short height;
  int   win;

  win = xge_current_win;
  xge_RemoveErrorMessage ();
  xge_SetWindow ( win );

  xge_usr_msg_key = key;
  xge_errmsg_win = xge_current_win;
  xge_errmsg_msgtext = NULL;
  xge_info_msgtext   = msglines;
  for ( xge_InfoNLines = xge_InfoMaxLength = 0;
        *msglines;
         msglines++, xge_InfoNLines++ )
    xge_InfoMaxLength = max ( xge_InfoMaxLength, strlen(*msglines) );

  height = (short)(44 + 16*xge_InfoNLines);
  xge_SetupEdRect ( xge_current_win, xge_errmsg_edr, 0, 0, xge_WIDTH-20, height,
                (short)((xge_current_width-xge_WIDTH+20)/2),
                (short)((xge_current_height-height)/2),
                &xge_ErrorMessageMsg, &xge_RedrawInfoMessage );
  xge_errmsg_edr->state = xgestate_MESSAGE;
  xge_AddPopup ( xge_errmsg_edr );
  xge_GrabFocus ( xge_errmsg_edr, true );
  xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
  xge_SetOtherWindowsCursor ( xgeCURSOR_CIRCLE );
  if ( xgeevent.type != ConfigureNotify )
    xgeRaiseWindow ( );
} /*xge_DisplayInfoMessage*/

