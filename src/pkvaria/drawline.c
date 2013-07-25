
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */
#include <stdlib.h>

#include "pkvaria.h"
#include "pkgeom.h"


/* Bresenham algorithm of drawing straight lines */

void   (*_pkv_OutputPixels)(const xpoint *buf, int n) = NULL;
xpoint *_pkv_pixbuf = NULL;
int    _pkv_npix = 0;

void _pkv_InitPixelBuffer ( void )
{
  _pkv_pixbuf = pkv_GetScratchMem ( PKV_BUFSIZE*sizeof(xpoint) );
  _pkv_npix = 0;
} /*_pkv_InitPixelBuffer*/

void _pkv_DestroyPixelBuffer ( void )
{
  if ( _pkv_pixbuf ) {
    pkv_FreeScratchMem ( PKV_BUFSIZE*sizeof(xpoint) );
    _pkv_pixbuf = NULL;
  }
} /*_pkv_DestroyPixelBuffer*/

void _pkv_DrawLine ( int x1, int y1, int x2, int y2 )
{
  int dx, dy, dx2, dy2, g, h, c;

  dx = x2-x1;
  if ( dx > 0 ) g = 1; else g = -1;
  dx = abs(dx);  dx2 = dx+dx;
  dy = y2-y1;
  if ( dy > 0 ) h = 1; else h = -1;
  dy = abs(dy);  dy2 = dy+dy;

  if ( dx > dy ) {
    c = -dx;
    while ( x1 != x2 ) {
      PKV_SETPIXEL ( x1, y1 );
      c += dy2;
      if ( c > 0 ) { y1 += h;  c -= dx2; }
      x1 += g;
    }
  }
  else {
    c = -dy;
    while ( y1 != y2 ) {
      PKV_SETPIXEL ( x1, y1 );
      c += dx2;
      if ( c > 0 ) { x1 += g;  c -= dy2; }
      y1 += h;
    }
  }
} /*_pkv_DrawLine*/

void pkv_DrawLine ( int x1, int y1, int x2, int y2,
                    void (*output)(const xpoint *buf, int n) )
{
  _pkv_InitPixelBuffer ();
  if ( !_pkv_pixbuf )
    return;
  _pkv_OutputPixels = output;
  _pkv_DrawLine ( x1, y1, x2, y2 );
  PKV_FLUSH
  _pkv_DestroyPixelBuffer ();
} /*pkv_DrawLine*/

