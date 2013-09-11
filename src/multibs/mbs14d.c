
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Modified by */
/* Pawel Szklarz, 23.03.2004 - added processing of Bessel end conditions    */
/* Pawel Szklarz, 07.04.2004 - added processing of given second derivative  */
/*                             end conditions                               */
/* Przemyslaw Kiciak, 9.04.2004 - added zero first and second order         */
/*                                end conditions                            */
/* Pawel Szklarz, 06.05.2004 - added processing of not-a-knot end condition */
/* Pawel Szklarz, 25.05.2004 - added processing of the third order          */
/*                             derivative end conditions                    */
/* Przemyslaw Kiciak, 30.06.2004 - error correction and some editing        */
/* Przemyslaw Kiciak, 13.07.2007 - removing unused parameters               */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* constructing cubic B-spline curves of interpolation */

/* first come the procedures setting up various boundary condition */
/* equations */

static void _mbs_BSCubicFirstDerLeftd ( double *bsknots,
                                        double *a, double *b, double *c,
                                        int ncurves, int spdimen, int ypitch,
                                        const double *y, int bspitch, double *d,
                                        boolean bc )
{
  if ( bc == BS3_BC_FIRST_DER )       /* given derivative at the end point */
    pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[spdimen] );
  else                                /* zero derivative at the end point */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[spdimen] );
  a[1] = -(b[1] = 3.0/(bsknots[4]-bsknots[1]));  c[1] = 0.0;
} /*_mbs_BSCubicFirstDerLeftd*/

static void _mbs_BSCubicFirstDerRightd ( int n, double *bsknots,
                                         double *a, double *b, double *c,
                                         int ncurves, int spdimen, int ypitch,
                                         const double *y, int bspitch, double *d,
                                         boolean bc )
{
  if ( bc == BS3_BC_FIRST_DER )       /* given derivative at the end point */
    pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[(n-1)*spdimen] );
  else                                /* zero derivative at the end point  */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[(n-1)*spdimen] );
  a[n-1] = 0.0;  b[n-1] = -(c[n-1] = 3.0/(bsknots[n+3]-bsknots[n]));
} /*_mbs_BSCubicFirstDerRightd*/

static void _mbs_BSCubicSecondDerLeftd ( double *bsknots,
                                         double *a, double *b, double *c,
                                         int ncurves, int spdimen, int ypitch,
                                         const double *y, int bspitch, double *d,
                                         boolean bc )
{
  double aa, bb, h;

  if ( bc == BS3_BC_SECOND_DER )      /* given second order derivative */
    pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[spdimen] );
  else                                /* zero second order derivative  */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[spdimen] );
  h = 6.0/(bsknots[4]-bsknots[2]);
  aa = h/(bsknots[5]-bsknots[2]);
  bb = h/(bsknots[4]-bsknots[1]);
  a[1] = bb;  b[1] = -(aa+bb);  c[1] = aa;
} /*_mbs_BSCubicSecondDerLeftd*/

static void _mbs_BSCubicSecondDerRightd ( int n,
                                          double *bsknots,
                                          double *a, double *b, double *c,
                                          int ncurves, int spdimen, int ypitch,
                                          const double *y, int bspitch, double *d,
                                          boolean bc )
{
  double aa, bb, h;

  if ( bc == BS3_BC_SECOND_DER )      /* given second order derivative */
    pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[(n-1)*spdimen] );
  else                                /* zero second order derivative */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[(n-1)*spdimen] );
  h = 6.0/(bsknots[n+2]-bsknots[n]);
  aa = h/(bsknots[n+3]-bsknots[n]);
  bb = h/(bsknots[n+2]-bsknots[n-1]);
  a[n-1] = bb;  b[n-1] = -(aa+bb);  c[n-1] = aa;
} /*_mbs_BSCubicSecondDerRightd*/

static void _mbs_BSCubicThirdDerLeftd ( double *bsknots,
                                        double *a, double *b, double *c,
                                        int ncurves, int spdimen, int ypitch,
                                        const double *y, int bspitch, double *d,
                                        boolean bc )
{
  double aa, bb, cc, dd;
  int   i;

  if ( bc == BS3_BC_THIRD_DER )       /* given third order derivative */
     pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[spdimen] );
  else                                /* zero third order derivative  */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[spdimen] );
  aa = bsknots[4]-bsknots[3];
  bb = bsknots[5]-bsknots[3];
  cc = bsknots[6]-bsknots[3];
  a[1] = -6.0/(aa*aa*aa);
  b[1] = ((6.0/aa+6.0/bb)/aa+6.0/(bb*bb))/aa;
  c[1] = -(6.0/aa+6.0/bb+6.0/cc)/(bb*aa);
  dd = 6.0/(aa*bb*cc);
  dd /= c[2];
  b[1] -= dd*a[2];
  c[1] -= dd*b[2];
  for ( i = 0; i < spdimen; i++ )
    d[spdimen+i] -= dd*d[2*spdimen+i];
} /*_mbs_BSCubicThirdDerLeftd*/

static void _mbs_BSCubicThirdDerRightd ( int n,
                                         double *bsknots,
                                         double *a, double *b, double *c,
                                         int ncurves, int spdimen, int ypitch,
                                         const double *y, int bspitch, double *d,
                                         boolean bc )
{
  double aa, bb, cc, dd;
  int   i;

  if ( bc == BS3_BC_THIRD_DER )       /* given third order derivative */
    pkv_Selectd ( ncurves, spdimen, ypitch, bspitch, y, &d[(n-1)*spdimen] );
  else                                /* zero third order derivative */
    pkv_ZeroMatd ( ncurves, spdimen, bspitch, &d[(n-1)*spdimen] );
  aa = bsknots[n+1]-bsknots[n];
  bb = bsknots[n+1]-bsknots[n-1];
  cc = bsknots[n+1]-bsknots[n-2];
  dd = -6.0/(aa*bb*cc);
  a[n-1] = (6.0/aa+6.0/bb+6.0/cc)/(bb*aa);
  b[n-1] = -((6.0/aa+6.0/bb)/aa+6.0/(bb*bb))/aa;
  c[n-1] = 6.0/(aa*aa*aa);
  dd /= a[n-2];
  a[n-1] -= dd*b[n-2];
  b[n-1] -= dd*c[n-2];
  for ( i = 0; i < spdimen; i++ )
    d[(n-1)*spdimen+i] -= dd*d[(n-2)*spdimen+i];
} /*_mbs_BSCubicThirdDerRightd*/

static void _mbs_BSCubicBesselLeftd ( double *bsknots,
                                      double *a, double *b, double *c,
                                      int ncurves, int spdimen, int xpitch,
                                      const double *x, int bspitch, double *d )
{
  void  *sp;
  double *ybctemp, *ybctemp1;
  double h0, h1;

  sp = pkv_GetScratchMemTop ();
  a[1] = -(b[1] = 3.0/(bsknots[4]-bsknots[1]));  c[1] = 0.0;
  h0 = bsknots[4]-bsknots[3];
  h1 = bsknots[5]-bsknots[4];
  ybctemp = pkv_GetScratchMemd ( ncurves*spdimen );
  ybctemp1 = pkv_GetScratchMemd ( ncurves*spdimen );
  pkn_MatrixMDifferenced ( ncurves, spdimen, xpitch, &x[spdimen], xpitch, x,
                           2.0 + h1/h0, spdimen, ybctemp );
  pkn_MatrixMDifferenced ( ncurves, spdimen, xpitch, &x[2*spdimen],
                           xpitch, &x[spdimen], h0/h1, spdimen, ybctemp1 );
  pkn_MatrixMDifferenced ( ncurves, spdimen, spdimen, ybctemp,
                           spdimen, ybctemp1, 1.0/(h0+h1),
                           bspitch, &d[spdimen] );
  pkv_SetScratchMemTop ( sp );
} /*_mbs_BSCubicBesselLeftd*/

static void _mbs_BSCubicBesselRightd ( int n,
                                       double *bsknots,
                                       double *a, double *b, double *c,
                                       int ncurves, int spdimen, int xpitch,
                                       const double *x, int bspitch, double *d )
{
  void  *sp;
  double *ybctemp, *ybctemp1;
  double hn1, hn2;

  sp = pkv_GetScratchMemTop ();
  b[n-1] = -(c[n-1] = 3.0/(bsknots[n+3]-bsknots[n]));
  a[n-1] = 0.0;
  hn1 = bsknots[n-2]-bsknots[n-3];
  hn2 = bsknots[n-3]-bsknots[n-4];
  ybctemp = pkv_GetScratchMemd ( ncurves*spdimen );
  ybctemp1 = pkv_GetScratchMemd ( ncurves*spdimen );
  pkn_MatrixMDifferenced ( ncurves, spdimen, xpitch, &x[spdimen*(n-2)],
                           xpitch, &x[spdimen*(n-3)],
                           2.0 + hn2/hn1, spdimen, ybctemp );
  pkn_MatrixMDifferenced ( ncurves, spdimen, xpitch, &x[spdimen*(n-3)],
                           xpitch, &x[spdimen*(n-4)], hn1/hn2,
                           spdimen, ybctemp1 );
  pkn_MatrixMDifferenced ( ncurves, spdimen, spdimen, ybctemp,
                           spdimen, ybctemp1, 1.0/(hn1+hn2),
                           bspitch, &d[(n-1)*spdimen] );
  pkv_SetScratchMemTop ( sp );
} /*_mbs_BSCubicBesselRightd*/

static void _mbs_BSCubicNotAKnotLeftd ( double *a, double *b, double *c,
                                        int ncurves, int spdimen,
                                        int bspitch, double *d, double d1 )
{
  int   i, j;
  double w;

  w = d1/c[2];
  b[1] = b[1]-a[2]*w;
  c[1] = c[1]-b[2]*w;
  for ( j = 0;  j < ncurves;  j++, d += bspitch )
    for ( i = 0; i < spdimen; i++ )
      d[spdimen+i] -= w*d[2*spdimen+i];
} /*_mbs_BSCubicNotAKnotLeftd*/

static void _mbs_BSCubicNotAKnotRightd ( int n,
                                         double *a, double *b, double *c,
                                         int ncurves, int spdimen,
                                         int bspitch, double *d, double dn1 )
{
  int   i, j;
  double w;

  w = dn1/a[n-2];
  a[n-1] -= w*b[n-2];
  b[n-1] -= w*c[n-2];
  for ( j = 0;  j < ncurves;  j++, d += bspitch )
    for ( i = 0; i < spdimen; i++ )
      d[(n-1)*spdimen+i] -= w*d[(n-2)*spdimen+i];
} /*_mbs_BSCubicNotAKnotRightd*/

boolean mbs_multiBSCubicInterpd ( int lastinterpknot, double *interpknots,
                                  int ncurves, int spdimen, int xpitch,
                                  const double *x,
                                  int ypitch,
                                  char bcl, const double *ybcl,
                                  char bcr, const double *ybcr,
                                  int *lastbsknot, double *bsknots,
                                  int bspitch, double *ctlpoints )
{
  void   *sp;
  int    i, j, l, n, nnz;
  int    lastknot;
  double bfv[4], m;
  double *a, *b, *c, *d;
  int    size_a;
  int    notaknot;
  double d1 = 0.0, dn1 = 0.0;  /* to suppress a warning */

  sp = pkv_GetScratchMemTop ();
  notaknot = 0;
  if ( bcl != BS3_BC_NOT_A_KNOT && bcr != BS3_BC_NOT_A_KNOT ) {
                                       /* setup curve knot sequence */
    *lastbsknot = lastknot = lastinterpknot+6-notaknot;
    n = lastknot-4;

    memcpy ( &bsknots[3], interpknots, (lastinterpknot+1)*sizeof(double) );
    bsknots[1] = bsknots[2] = bsknots[3];
    bsknots[0] = bsknots[1] - 0.5;
    bsknots[lastknot-1] = bsknots[lastknot-2] = bsknots[lastknot-3];
    bsknots[lastknot] = bsknots[lastknot-1] + 0.5;

                                       /* setup the system of equations */
    a = pkv_GetScratchMemd ( size_a = lastinterpknot+3 );
    b = pkv_GetScratchMemd ( size_a );
    c = pkv_GetScratchMemd ( size_a );
    if ( !a || !b || !c ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }

    b[0] = 1.0;  c[0] = 0.0;
    for ( i = 1; i < lastinterpknot; i++ ) {
      mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[i], &j, &nnz, bfv );
      a[i+1] = bfv[0];  b[i+1] = bfv[1];  c[i+1] = bfv[2];
    }
    a[n] = 0.0;  b[n] = 1.0;

    pkv_Selectd ( ncurves, spdimen, xpitch, bspitch, x, ctlpoints );
    pkv_Selectd ( ncurves, (lastinterpknot-1)*spdimen, xpitch, bspitch,
                  &x[spdimen], &ctlpoints[2*spdimen] );
    pkv_Selectd ( ncurves, spdimen, xpitch, bspitch,
                  &x[lastinterpknot*spdimen],
                  &ctlpoints[(lastknot-4)*spdimen] );
  }
  else if ( bcl == BS3_BC_NOT_A_KNOT && bcr == BS3_BC_NOT_A_KNOT ) {
    notaknot = 2;
                                /* setup curve knot sequence */
    *lastbsknot = lastknot = lastinterpknot+6-notaknot;
    n = lastknot-4;

    memcpy ( &bsknots[3], interpknots,sizeof(double) );
    memcpy ( &bsknots[4], &interpknots[2], (lastinterpknot-2)*sizeof(double) );
    memcpy ( &bsknots[lastinterpknot+1], &interpknots[lastinterpknot],
             sizeof(double) );

    bsknots[1] = bsknots[2] = bsknots[3];
    bsknots[0] = bsknots[1] - 0.5;
    bsknots[lastknot-1] = bsknots[lastknot-2] = bsknots[lastknot-3];
    bsknots[lastknot] = bsknots[lastknot-1] + 0.5;

                                       /* setup the system of equations */
    a = pkv_GetScratchMemd ( size_a = lastinterpknot+3-notaknot );
    b = pkv_GetScratchMemd ( size_a );
    c = pkv_GetScratchMemd ( size_a );
    if ( !a || !b || !c ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }

    b[0] = 1.0;  c[0] = 0.0;
    mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[1], &j, &nnz, bfv );
    a[1] = bfv[0];  b[1] = bfv[1];  c[1] = bfv[2];  d1 = bfv[3];

    for ( i = 2; i < lastinterpknot-1; i++ ) {
      mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[i], &j, &nnz, bfv );
      a[i] = bfv[0];  b[i] = bfv[1];  c[i] = bfv[2];
    }
    mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[lastinterpknot-1],
                       &j, &nnz, bfv );
    a[lastinterpknot-1] = bfv[1];
    b[lastinterpknot-1] = bfv[2];
    c[lastinterpknot-1] = bfv[3];
    dn1 = bfv[0];
    a[n] = 0.0;  b[n] = 1.0;

    pkv_Selectd ( ncurves, (lastinterpknot+1)*spdimen, xpitch, bspitch,
                  x, ctlpoints );
  }
  else if ( bcl == BS3_BC_NOT_A_KNOT ) {
    notaknot = 1;
                                  /* setup curve knot sequence */
    *lastbsknot = lastknot = lastinterpknot+6-notaknot;
    n = lastknot-4;
  
    memcpy ( &bsknots[3], interpknots,sizeof(double) );
    memcpy ( &bsknots[4], &interpknots[2],(lastinterpknot-1)*sizeof(double));

    bsknots[1] = bsknots[2] = bsknots[3];
    bsknots[0] = bsknots[1] - 0.5;
    bsknots[lastknot-1] = bsknots[lastknot-2] = bsknots[lastknot-3];
    bsknots[lastknot] = bsknots[lastknot-1] + 0.5;

                                       /* setup the system of equations */
    a = pkv_GetScratchMemd ( size_a = lastinterpknot+3-notaknot );
    b = pkv_GetScratchMemd ( size_a );
    c = pkv_GetScratchMemd ( size_a );
    if ( !a || !b || !c ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }

    b[0] = 1.0;  c[0] = 0.0;
    mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[1], &j, &nnz, bfv );
      a[1] = bfv[0];  b[1] = bfv[1];  c[1] = bfv[2];  d1 = bfv[3];

    for ( i = 2; i < lastinterpknot; i++ ) {
      mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[i], &j, &nnz, bfv );
      a[i] = bfv[0];  b[i] = bfv[1];  c[i] = bfv[2];
    }
    a[n] = 0.0;  b[n] = 1.0;

    pkv_Selectd ( ncurves, (lastinterpknot)*spdimen, xpitch, bspitch,
                  x, ctlpoints );
    pkv_Selectd ( ncurves, spdimen, xpitch, bspitch,
                  &x[lastinterpknot*spdimen],
                  &ctlpoints[(lastknot-4)*spdimen] );
  }
  else {
    notaknot = 1;
                                /* setup curve knot sequence */
    *lastbsknot = lastknot = lastinterpknot+6-notaknot;
    n = lastknot-4;

    memcpy ( &bsknots[3], interpknots, (lastinterpknot-1)*sizeof(double) );
    memcpy ( &bsknots[lastinterpknot+2], &interpknots[lastinterpknot],
             sizeof(double) );

    bsknots[1] = bsknots[2] = bsknots[3];
    bsknots[0] = bsknots[1] - 0.5;
    bsknots[lastknot-1] = bsknots[lastknot-2] = bsknots[lastknot-3];
    bsknots[lastknot] = bsknots[lastknot-1] + 0.5;

                                       /* setup the system of equations */
    a = pkv_GetScratchMemd ( size_a = lastinterpknot+3-notaknot );
    b = pkv_GetScratchMemd ( size_a );
    c = pkv_GetScratchMemd ( size_a );
    if ( !a || !b || !c ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }

    b[0] = 1.0;  c[0] = 0.0;
    for ( i = 1; i < lastinterpknot-1; i++ ) {
      mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[i], &j, &nnz, bfv );
      a[i+1] = bfv[0];  b[i+1] = bfv[1];  c[i+1] = bfv[2];
    }
    mbs_deBoorBasisd ( 3, lastknot, bsknots, interpknots[lastinterpknot-1], &j,
                       &nnz, bfv );
    a[lastinterpknot] = bfv[1];
    b[lastinterpknot] = bfv[2];
    c[lastinterpknot] = bfv[3];
    dn1 = bfv[0];

    a[n] = 0.0;  b[n] = 1.0;

    pkv_Selectd ( ncurves, spdimen, xpitch, bspitch, x, ctlpoints );
    pkv_Selectd ( ncurves, (lastinterpknot)*spdimen, xpitch, bspitch,
                  &x[spdimen], &ctlpoints[2*spdimen] );
  }

 /* Setup boundary conditions */

  switch ( bcl ) {

case BS3_BC_FIRST_DER:           /* given derivative at the end point */
case BS3_BC_FIRST_DER0:          /* zero derivative at the end point  */
    _mbs_BSCubicFirstDerLeftd ( bsknots, a, b, c, ncurves, spdimen,
                                ypitch, ybcl, bspitch, ctlpoints, bcl );
    break;

case BS3_BC_SECOND_DER:          /* given second order derivative */
case BS3_BC_SECOND_DER0:         /* zero second order derivative  */
    _mbs_BSCubicSecondDerLeftd ( bsknots, a, b, c, ncurves, spdimen,
                                 ypitch, ybcl, bspitch, ctlpoints, bcl );
    break;

case BS3_BC_THIRD_DER:           /* given third order derivative */
case BS3_BC_THIRD_DER0:          /* zero third order derivative */
    _mbs_BSCubicThirdDerLeftd ( bsknots, a, b, c, ncurves, spdimen,
                                ypitch, ybcl, bspitch, ctlpoints, bcl );
    break;

case BS3_BC_BESSEL:              /* Bessel end condition */
    _mbs_BSCubicBesselLeftd ( bsknots, a, b, c, ncurves, spdimen,
                              xpitch, x, bspitch, ctlpoints );
    break;
    
case BS3_BC_NOT_A_KNOT:          /* not a knot condition */
    _mbs_BSCubicNotAKnotLeftd ( a, b, c, ncurves, spdimen,
                                bspitch, ctlpoints, d1 );
    break;

default:
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    goto failure;   /* unknown boundary condition */
  }

  switch ( bcr ) {

case BS3_BC_FIRST_DER:           /* given derivative at the end point */
case BS3_BC_FIRST_DER0:          /* zero derivative at the end point  */
    _mbs_BSCubicFirstDerRightd ( n, bsknots, a, b, c, ncurves, spdimen,
                                 ypitch, ybcr, bspitch, ctlpoints, bcr );
    break;

case BS3_BC_SECOND_DER:          /* given second order derivative */
case BS3_BC_SECOND_DER0:         /* zero second order derivative */
    _mbs_BSCubicSecondDerRightd ( n, bsknots, a, b, c, ncurves, spdimen,
                                  ypitch, ybcr, bspitch, ctlpoints, bcr );
    break;

case BS3_BC_THIRD_DER:          /* given third order derivative */
case BS3_BC_THIRD_DER0:         /* zero third order derivative */
    _mbs_BSCubicThirdDerRightd ( n, bsknots, a, b, c, ncurves, spdimen,
                                 ypitch, ybcr, bspitch, ctlpoints, bcr );
    break;

case BS3_BC_BESSEL:             /* Bessel end condition */
    _mbs_BSCubicBesselRightd ( n, bsknots, a, b, c, ncurves, spdimen,
                               xpitch, x, bspitch, ctlpoints );
    break;

case BS3_BC_NOT_A_KNOT:         /* not a knot condition */
    _mbs_BSCubicNotAKnotRightd ( n, a, b, c, ncurves, spdimen,
                                 bspitch, ctlpoints, dn1 );
    break;

default:
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    goto failure;   /* unknown boundary condition */
  }

              /* solve the system with the three-diagonal matrix; */
              /* Gaussian elimination without pivoting */
  for ( i = 0; i < n; i++ ) {
    m = a[i+1]/b[i];
    b[i+1] -= m*c[i];
    for ( l = 0, d = &ctlpoints[i*spdimen];  l < ncurves;  l++, d += bspitch )
      for ( j = 0;  j < spdimen;  j++ )
        d[j+spdimen] -= m*d[j];
  }

  for ( l = 0, d = &ctlpoints[n*spdimen];  l < ncurves;  l++, d+= bspitch )
    for ( j = 0; j < spdimen; j++ )
      d[j] /= b[n];

  for ( i = n-1; i >= 0; i-- )
    for ( l = 0, d = &ctlpoints[i*spdimen];  l < ncurves;  l++, d += bspitch )
      for ( j = 0; j < spdimen; j++ )
        d[j] = (d[j]-c[i]*d[j+spdimen])/b[i];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBSCubicInterpd*/

