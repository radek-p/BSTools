
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

void xge_DrawVShadedRect ( short w, short h, short x, short y,
                           xgecolour_int c0, xgecolour_int c1, short nb )
{
  short i, y0, y1;
  byte  r0, r1, r, g0, g1, g, b0, b1, b;
  float a0, a1;

  if ( nb < 2 ) nb = 2; else if ( nb > h ) nb = h;
  xge_GetPixelColour ( c0, &r0, &g0, &b0 );
  xge_GetPixelColour ( c1, &r1, &g1, &b1 );
  for ( i = 0, y0 = y;  i < nb;  i++, y0 = y1 ) {
    y1 = y + (short)(((float)i+0.5)/(float)(nb-1)*(float)h + 0.5);
    a0 = (float)i/(float)(nb-1);
    a1 = 1.0-a0;
    r = (byte)(a1*(float)r0 + a0*(float)r1);
    g = (byte)(a1*(float)g0 + a0*(float)g1);
    b = (byte)(a1*(float)b0 + a0*(float)b1);
    xgeSetForeground ( xge_PixelColour ( r, g, b ) );
    xgeFillRectangle ( w, y1-y0, x, y0 );
  }
} /*xge_DrawVShadedRect*/

void xge_DrawHShadedRect ( short w, short h, short x, short y,
                           xgecolour_int c0, xgecolour_int c1, short nb )
{
  short i, x0, x1;
  byte  r0, r1, r, g0, g1, g, b0, b1, b;
  float a0, a1;

  if ( nb < 2 ) nb = 2; else if ( nb > w ) nb = w;
  xge_GetPixelColour ( c0, &r0, &g0, &b0 );
  xge_GetPixelColour ( c1, &r1, &g1, &b1 );
  for ( i = 0, x0 = x;  i < nb;  i++, x0 = x1 ) {
    x1 = x + (short)(((float)i+0.5)/(float)(nb-1)*(float)w + 0.5);
    a0 = (float)i/(float)(nb-1);
    a1 = 1.0-a0;
    r = (byte)(a1*(float)r0 + a0*(float)r1);
    g = (byte)(a1*(float)g0 + a0*(float)g1);
    b = (byte)(a1*(float)b0 + a0*(float)b1);
    xgeSetForeground ( xge_PixelColour ( r, g, b ) );
    xgeFillRectangle ( x1-x0, h, x0, y );
  }
} /*xge_DrawHShadedRect*/

