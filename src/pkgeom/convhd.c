
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "convh.h"

boolean FindConvexHull2d ( int *n, point2d *p )
{
  double       *q;
  unsigned int *permut;
  int          i, j, nn;
  point2d      a;
  double       b, c;
  void         *top;

  top = pkv_GetScratchMemTop ();
  nn = *n;
  if (nn <= 2)
    goto failure;
  q = pkv_GetScratchMem ( nn*sizeof(double) );
  permut = pkv_GetScratchMem ( nn*sizeof(int) );

  if ( !q || !permut ) {
    PKV_SIGNALERROR ( LIB_GEOM, 2, ERRMSG_2 );
    goto failure;
  }
  for (i = 1; i < nn; i++) {
    if (p[i].y < p[0].y || (p[i].y == p[0].y && p[i].x > p[0].x))
    { a = p[i];  p[i] = p[0];  p[0] = a; }
  }
  q[0] = 0.0;
  for (i = 1; i < nn; i++)
    q[i] = pkv_SqAngle(p[i].x - p[0].x, p[i].y - p[0].y);

/* */
  for ( i = 0; i < nn; i++ )
    permut[i] = i;
  if ( pkv_SortKernel ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                        0, nn, q, permut ) )
    pkv_SortPermute ( sizeof(point2d), nn, p, permut );
  else goto failure;
/* */

/*
  for (i = 2; i < nn; i++) {
    j = i;
    while (q[j] < q[j-1]) {
      a = p[j];  p[j] = p[j-1];  p[j-1] = a;
      b = q[j];  q[j] = q[j-1];  q[j-1] = b;
      j--;
    }
  }
*/
  b = pkv_SqAngle(p[1].x - p[0].x, p[1].y - p[0].y);
  j = 2;
  while ( j < nn ) {
    c = pkv_SqAngle ( p[j].x - p[j-1].x, p[j].y - p[j-1].y );
    if ( c > b ) {
      b = c;
      j++;
    }
    else
    {
      memmove( &(p[j-1]), &(p[j]), (nn-j)*sizeof(point2d) );
      nn--;
      if (j > 2)
        j--;
      b = pkv_SqAngle ( p[j-1].x - p[j-2].x, p[j-1].y - p[j-2].y );
    }
  }
  *n = nn;
  pkv_SetScratchMemTop ( top );
  return true;

failure:
  pkv_SetScratchMemTop ( top );
  return false;
} /*FindConvexHull2d*/

