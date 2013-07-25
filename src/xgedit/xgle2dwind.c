
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


void xgle_2DwindDrawCursorPos ( xge_2Dwind *_2Dwin, short x, short y )
{
  xge_widget *er;
  char       s[10];
  point2d    p, q;
  int        xx;

  if ( _2Dwin->display_coord ) {
    er = _2Dwin->er;
    xgle_SetIdentMapping ( er );
    SetPoint2d ( &p, x, y );
    CameraUnProjectPoint2d ( &_2Dwin->CPos, &p, &q );
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
} /*xgle_2DwindDrawCursorPos*/

void xgle_2DwindDrawAxes ( xge_2Dwind *_2Dwin )
{
  xge_widget *er;
  point2d    p, q;

  er = _2Dwin->er;
  SetPoint2d ( &p, 0.0, 0.0 );
  CameraProjectPoint2d ( &_2Dwin->CPos, &p, &q );
  xgle_SetIdentMapping ( er );
  glColor3fv ( xglec_Green5 );
  if ( q.x >= er->x && q.x < er->x+er->w )
    xgleDrawLine ( q.x, er->y, q.x, er->y+er->h );
  if ( q.y >= er->y && q.y < er->y+er->h )
    xgleDrawLine ( er->x, q.y, er->x+er->w, q.y );
} /*xgle_2DwindDrawAxes*/

