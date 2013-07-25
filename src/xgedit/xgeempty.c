
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
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


void xge_DrawEmpty ( xge_widget *er, boolean onscreen )
{
} /*xge_DrawEmpty*/

boolean xge_EmptyMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  return true;
} /*xge_EmptyMsg*/

xge_widget *xge_NewEmptyWidget ( char window_num, xge_widget *prev, int id,
                                 short w, short h, short x, short y )
{
  return xge_NewWidget ( window_num, prev, id, w, h, x, y, NULL, NULL,
                         xge_EmptyMsg, xge_DrawEmpty );
} /*xge_NewEmptyWidget*/

