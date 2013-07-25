
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
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

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

void xge_GrabFocus ( xge_widget *er, boolean all )
{
  int i, fsp;

  if ( er ) {
    if ( all ) {
      for ( i = 0; i < xge_winnum; i++ ) {
        fsp = xge_windesc[i].fsp;
        if ( fsp < xge_FOCUS_DEPTH ) {
          xge_windesc[i].cursorstack[fsp] = xge_windesc[i].cursor;
          xge_windesc[i].focusstack[fsp] = er;
          xge_windesc[i].fsp = fsp + 1;
        }
      }
    }
    else {
      i = er->window_num;
      fsp = xge_windesc[i].fsp;
      if ( fsp < xge_FOCUS_DEPTH ) {
        xge_windesc[i].cursorstack[fsp] = xge_windesc[i].cursor;
        xge_windesc[i].focusstack[fsp] = er;
        xge_windesc[i].fsp = fsp + 1;
      }
    }
  }
} /*xge_GrabFocus*/

void xge_ReleaseFocus ( xge_widget *er )
{
  int i, j, fsp;

        /* search the grabbed focus stacks of all windows */
  for ( i = 0; i < xge_winnum; i++ ) {
    fsp = xge_windesc[i].fsp;
    for ( j = fsp-1; j >= 0; j-- ) {
      if ( xge_windesc[i].focusstack[j] == er ) {
        xge_windesc[i].fsp = j;
        if ( xge_windesc[i].cursor != xge_windesc[i].cursorstack[j] )
          xge_SetWindowCursor ( i, xge_windesc[i].cursorstack[j] );
        xge_windesc[i].focusstack[j] = NULL;
        break;
      }
    }
  }
} /*xge_ReleaseFocus*/

xge_widget *xge_GetFocusWidget ( char win )
{
  int fsp;

  if ( win < 0 || win >= xge_winnum )
    return NULL;
  fsp = xge_windesc[(int)win].fsp;
  if ( fsp > 0 )
    return xge_windesc[(int)win].focusstack[fsp-1];
  else
    return NULL;
} /*xge_GetFocusWidget*/

