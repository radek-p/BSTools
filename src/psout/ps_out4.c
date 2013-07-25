
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <memory.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"

#include "psprivate.h"

void ps_Draw_BezierCf ( const point2f *p, int n )
{
#define strokes  50
  void    *sp;
  point2f *q, *r;
  int     i, j, k;
  float   s, t;

  if (n == 0)
    return;
  if (n == 1) {
    ps_Draw_Line(p[0].x, p[0].y, p[1].x, p[1].y);
    return;
  }
  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMem ( (strokes + 1)*sizeof(point2f) );
  r = pkv_GetScratchMem ( (n + 1)*sizeof(point2f) );
  if ( q && r ) {
    q[0] = p[0];
    for (i = 1; i <= strokes; i++) {
      t = (float)i/(float)strokes;   /* De Casteljau algorithm */
      s = (float)1.0 - t;
      memmove((char *)r, (char *)p, (n + 1L) * sizeof(point2f));
      k = n;
      while (k > 0) {
	for (j = 0; j <= k - 1; j++) {
	  r[j].x = s * r[j].x + t * r[j+1].x;
	  r[j].y = s * r[j].y + t * r[j+1].y;
	}
	k--;
      }
      q[i] = r[0];
    }
    ps_Draw_Polyline2f ( strokes+1, q );
  }
  pkv_SetScratchMemTop ( sp );
#undef strokes
} /*ps_Draw_BezierCf*/

void ps_Draw_BezierCd ( const point2d *p, int n )
{
#define strokes  50
  void     *sp;
  point2d  *q, *r;
  int      i, j, k;
  double   s, t;

  if (n == 0)
    return;
  if (n == 1) {
    ps_Draw_Line ( (float)p[0].x, (float)p[0].y, (float)p[1].x, (float)p[1].y );
    return;
  }
  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMem ( (strokes + 1)*sizeof(point2d) );
  r = pkv_GetScratchMem ( (n + 1)*sizeof(point2d) );
  if ( q && r ) {
    q[0] = p[0];
    for (i = 1; i <= strokes; i++) {
      t = (double)i / strokes;   /* De Casteljau algorithm */
      s = 1.0 - t;
      memmove((char *)r, (char *)p, (n + 1L) * sizeof(point2d));
      k = n;
      while (k > 0) {
	for (j = 0; j <= k - 1; j++) {
	  r[j].x = s * r[j].x + t * r[j+1].x;
	  r[j].y = s * r[j].y + t * r[j+1].y;
	}
	k--;
      }
      q[i] = r[0];
    }
    ps_Draw_Polyline2d ( strokes+1, q );
  }
  pkv_SetScratchMemTop ( sp );
#undef strokes
} /*ps_Draw_BezierCd*/

