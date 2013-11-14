
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgedit.h"
#include "xgledit.h"
#include "xgeprivate.h"


void xgle_T2KnotWinfDrawCursorPos ( xge_T2KnotWinf *T2win, short x, short y )
{
  xge_widget *er;
  char       s[10];
  point2f    p, q;
  int        xx;

  if ( T2win->display_coord ) {
    er = T2win->er;
    xgle_SetIdentMapping ( er );
    SetPoint2f ( &p, x, y );
    CameraUnProjectPoint2f ( &T2win->CPos, &p, &q );
    glColor3fv ( xglec_Green5 );
    xgleDrawLine ( x, er->y, x, er->y+er->h );
    xgleDrawLine ( er->x, y, er->x+er->w, y );
    glColor3fv ( xglec_Green );
    sprintf ( s, "%5.3f", q.x );
    xx = er->x+er->w-strlen(s)*6;
    x = (short)(min ( x+2, xx ));
    xgleDrawString ( s, x, er->y+er->h-4 );
    sprintf ( s, "%5.3f", q.y );
    y = (short)(max ( y, er->y+10 ));
    xgleDrawString ( s, er->x, y );
  }
} /*xgle_T2KnotWinfDrawCursorPos*/

