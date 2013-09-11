
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* constructing closed cubic B-spline curves of interpolation */

boolean mbs_multiBSCubicClosedInterpd ( int lastinterpknot, double *interpknots,
                               int ncurves, int spdimen, int xpitch,
                               const double *x,
                               int *lastbsknot, double *bsknots,
                               int bspitch, double *ctlpoints )
{
  void  *sp;
  int   i, j, l, n, nnz;
  int   lastknot, clcK;
  double bfv[4], m, clcT;
  double *a, *b, *c, *d;
  
  sp = pkv_GetScratchMemTop ();
                                       /* setup curve knot sequence */
  *lastbsknot = lastknot = lastinterpknot+6;
  n = lastinterpknot;
  clcK = lastinterpknot;
  clcT = interpknots[clcK]-interpknots[0];

  memcpy ( &bsknots[3], interpknots, (n+1)*sizeof(double) );
  for ( i = 0; i < 3; i++ )
    bsknots[i] = bsknots[i+clcK]-clcT;
  for ( i = 4; i <= 6; i++ )
    bsknots[i+clcK] = bsknots[i]+clcT;
                                       /* setup the system of equations */
  a = pkv_GetScratchMemd ( 3*n );
  if ( !a ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  b = &a[n];
  c = &b[n];

  for ( i = 0; i < lastinterpknot; i++ ) {
    mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[i], &j, &nnz, bfv );
    a[i] = bfv[0];  b[i] = bfv[1];  c[i] = bfv[2];
  }
  pkv_Selectd ( ncurves, n*spdimen, xpitch, bspitch,
                x, &ctlpoints[spdimen] );

              /* solve the system with the cyclic three-diagonal matrix; */
              /* Gaussian elimination without pivoting */
  for ( i = 0; i < n-2; i++ ) {
    m = a[i+1]/b[i];
    b[i+1] -= m*c[i];
    a[i+1] = -m*a[i];
    for ( l = 0, d = &ctlpoints[spdimen];  l < ncurves;  l++, d += bspitch )
      for ( j = 0; j < spdimen; j++ )
        d[(i+1)*spdimen+j] -= m*d[i*spdimen+j];
    m = c[n-1]/b[i];
    c[n-1] = -m*c[i];
    b[n-1] -= m*a[i];
    for ( l = 0, d = &ctlpoints[spdimen]; l < ncurves; l++, d += bspitch )
      for ( j = 0; j < spdimen; j++ )
        d[(n-1)*spdimen+j] -= m*d[i*spdimen+j];
  }
  c[n-2] += a[n-2];
  a[n-1] += c[n-1];
  m = a[n-1]/b[n-2];
  b[n-1] -= m*c[n-2];
  for ( l = 0, d = &ctlpoints[spdimen];  l < ncurves;  l++, d+= bspitch )
    for ( j = 0; j < spdimen; j++ )
      d[(n-1)*spdimen+j] = (d[(n-1)*spdimen+j]-m*d[(n-2)*spdimen+j])/b[n-1];
  for ( l = 0, d = &ctlpoints[spdimen];  l < ncurves;  l++, d += bspitch )
    for ( j = 0; j < spdimen; j++ )
      d[(n-2)*spdimen+j] = (d[(n-2)*spdimen+j]-c[n-2]*d[(n-1)*spdimen+j])/b[n-2];
  for ( i = n-3; i >= 0; i-- )
    for ( l = 0, d = &ctlpoints[spdimen];  l < ncurves;  l++, d += bspitch )
      for ( j = 0; j < spdimen; j++ )
        d[i*spdimen+j] = (d[i*spdimen+j]-c[i]*d[(i+1)*spdimen+j]
                          -a[i]*d[(n-1)*spdimen+j])/b[i];
  pkv_Selectd ( ncurves, spdimen, bspitch, bspitch,
                &ctlpoints[n*spdimen], ctlpoints );
  pkv_Selectd ( ncurves, 2*spdimen, bspitch, bspitch,
                &ctlpoints[spdimen], &ctlpoints[(n+1)*spdimen] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSCubicClosedInterpd*/

