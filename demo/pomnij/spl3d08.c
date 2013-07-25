
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"

/* ///////////////////////////////////////////////////////////////////////// */
void ClearPointMarking ( int npoints, boolean *mkpoints )
{
  memset ( mkpoints, 0, npoints*sizeof(boolean) );
} /*ClearPointMarking*/

boolean FindNearestPoint ( int npoints, const point2d *rpoints,
                           short x, short y, short *dist, int *k )
{
  int   i, j;
  short d, e;

  e = (short)(*dist+1);
  j = -1;
  for ( i = 0; i < npoints; i++ ) {
    d = (short)(fabs(rpoints[i].x-(double)x)+fabs(rpoints[i].y-(double)y));
    if ( d < e ) {
      e = d;
      j = i;
    }
  }
  if ( j >= 0 ) {
    *k = j;
    *dist = e;
    return true;
  }
  else
    return false;
} /*FindNearestPoint*/

void MarkPoints ( int npoints, const point2d *rpoints, boolean *mkpoints,
                  const Box2s *sbox, boolean mk )
{
  int   i;
  short dist;

  if ( sbox->x1-sbox->x0 >= 2 || sbox->y1-sbox->y0 >= 2 ) {
    for ( i = 0; i < npoints; i++ )
      if ( rpoints[i].x >= sbox->x0 && rpoints[i].x <= sbox->x1 &&
           rpoints[i].y >= sbox->y0 && rpoints[i].y <= sbox->y1 )
        mkpoints[i] = mk;
  }
  else {
    dist = xge_MINDIST;
    if ( FindNearestPoint ( npoints, rpoints, sbox->x0, sbox->y0,
                            &dist, &i ) )
      mkpoints[i] = mk;
  }
} /*MarkPoints*/

