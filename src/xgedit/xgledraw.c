
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgledit.h"
#include "xgeprivate.h"

void xgleDrawPoint ( int x, int y )
{
  glBegin ( GL_POINTS );
    glVertex2i ( x, y );
  glEnd ();
} /*xgleDrawPoint*/

void xgleDrawPoints ( int n, XPoint *p )
{
  int i;

  glBegin ( GL_POINTS );
    for ( i = 0; i < n; i++ )
      glVertex2i ( p[i].x, p[i].y );
  glEnd ();
} /*xgleDrawPoints*/

void xgleDrawLine ( int x0, int y0, int x1, int y1 )
{
  glBegin ( GL_LINES );
    glVertex2i ( x0, y0 );
    glVertex2i ( x1, y1 );
  glEnd ();
} /*xgleDrawLine*/

void xgleDrawRectangle ( int w, int h, int x, int y )
{
  glBegin ( GL_LINE_LOOP );
    glVertex2i ( x, y );
    glVertex2i ( x+h, y );
    glVertex2i ( x+h, y+h );
    glVertex2i ( x, y+h );
  glEnd ();
} /*xgleDrawRectangle*/

