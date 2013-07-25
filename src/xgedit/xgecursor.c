
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2011                            */
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

int xge_NewCursor ( int shape )
{
  Pixmap source, mask;
  XColor black = {0,0,0,0,DoRed|DoGreen|DoBlue,0};
  char data[32] = {0};

  if ( xge_cursnum < xge_MAX_CURSORS ) {
    if ( shape >= 0 )
      xgecursor[xge_cursnum] = XCreateFontCursor ( xgedisplay, shape );
    else {  /* create an invisible cursor */
      source = XCreateBitmapFromData ( xgedisplay, xgewindow, data, 16, 16 );
      mask = XCreateBitmapFromData ( xgedisplay, xgewindow, data, 16, 16 );
      xgecursor[xge_cursnum] = XCreatePixmapCursor ( xgedisplay, source, mask,
                                                     &black, &black, 0, 0 );
      XFreePixmap ( xgedisplay, mask );
      XFreePixmap ( xgedisplay, source );
    }
    xge_cursnum ++;
    return xge_cursnum-1;
  }
  else return -1;
} /*xge_NewCursor*/

void xge_SetWindowCursor ( int win, Cursor cursor )
{
  xge_windesc[win].cursor = cursor;
  XDefineCursor ( xgedisplay, xge_windesc[win].thewindow, cursor );
} /*xge_SetWindowCursor*/

void xge_SetCurrentWindowCursor ( Cursor cursor )
{
  xge_SetWindowCursor ( xge_current_win, cursor );
} /*xge_SetCurrentWindowCursor*/

void xge_SetOtherWindowsCursor ( Cursor cursor )
{
  int i;

  for ( i = 0; i < xge_winnum; i++ )
    if ( i != xge_current_win )
      xge_SetWindowCursor ( i, cursor );
} /*xge_SetOtherWindowsCursor*/

