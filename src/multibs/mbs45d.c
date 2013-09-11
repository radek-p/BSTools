
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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

#include "msgpool.h"

/* /////////////////////////////////////////// */
static int nsigns ( int n, const double *a )
{
  int     i, ch;

  ch = 0;
  if ( a[0] < 0.0 ) {
    for ( i = 1; i < n; i++ )
      if ( a[i] >= 0.0 ) ch = 1;
      else if ( ch && a[i] < 0.0 ) return 2;
  }
  else {
    for ( i = 1; i < n; i++ )
      if ( a[i] < 0.0 ) ch = 1;
      else if ( ch && a[i] >= 0.0 ) return 2;
  }
  return ch;
} /*nsigns*/

static double *_coeff;
static int _deg;

static double poly ( double x )
{
  double f;

  mbs_BCHornerC1d ( _deg, _coeff, x, &f );
  return f;
} /*poly*/

boolean mbs_FindPolynomialZerosd ( int degree, const double *coeff,
                                   int *nzeros, double *zeros, double eps )
{
  void    *sp;
  double   *a, *b, x;
  int     nz, stp;
  boolean error;

  sp = pkv_GetScratchMemTop ();

        /* push the interval [0,1] and the related polynomial coefficients */
  if ( !(a = pkv_GetScratchMemd ( degree+3 )) )
    goto failure;
  _deg = degree;
  memcpy ( &a[2], coeff, (degree+1)*sizeof(double) );
  a[0] = 0.0;
  a[1] = 1.0;
  stp = 1;

  nz = 0;
  do {
    stp --;
    switch ( nsigns ( degree+1, &a[2] ) ) {
  case 0:    /* no sign change - reject the interval */
      pkv_FreeScratchMemd ( degree+3 );
      a -= degree+3;
      break;

  case 1:    /* one sign change - try to solve */
      _coeff = &a[2];
      x = pkn_Illinoisd ( poly, 0.0, 1.0, eps, &error );
      if ( error )
        goto divide;
      zeros[nz++] = (1.0-x)*a[0] + x*a[1];
      if ( nz > degree )
        goto failure;
      pkv_FreeScratchMemd ( degree+3 );
      a -= degree+3;
      break;

  default:   /* more than one sign change - bisect and push */
divide:
      if ( !(b = pkv_GetScratchMemd ( degree+3 )) )
        goto failure;
      stp += 2;
      b[0] = a[0];
      a[0] = b[1] = 0.5*(b[0]+a[1]);
      mbs_BisectBC1d ( degree, &a[2], &b[2] );
      a = b;
      break;
    }
  } while ( stp != 0 );

  *nzeros = nz;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_FindPolynomialZerosd*/


/* /////////////////////////////////////////// */
#define EPS 1.0e-6

boolean mbs_ClipBC2d ( int ncplanes, const vector3d *cplanes,
                       int degree, const point2d *cpoints,
                       void (*output) (int degree, const point2d *cpoints) )
{
  void    *sp;
  point2d *p, *q;
  double   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemd ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemd ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = DotProduct2d ( &cpoints[i], (point2d*)&cplanes[j] ) +
               cplanes[j].z;

      if ( !mbs_FindPolynomialZerosd ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                   0, nt+1,  t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point2d) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = 0.5*(t[i]+t[i+1]);
      mbs_BCHornerC2d ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = DotProduct2d ( p, (vector2d*)&cplanes[j] ) + cplanes[j].z;
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point2d) );
        mbs_DivideBC2d ( degree, t[i], p, q );
        mbs_DivideBC2d ( degree, (t[i+1]-t[i])/(1.0-t[i]), p, q );
        output ( degree, q );
      }
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_ClipBC2d*/

boolean mbs_ClipBC2Rd ( int ncplanes, const vector3d *cplanes,
                        int degree, const point3d *cpoints,
                        void (*output) (int degree, const point3d *cpoints) )
{
  void    *sp;
  point3d *p, *q;
  double   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemd ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemd ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = DotProduct3d ( &cpoints[i], &cplanes[j] );

      if ( !mbs_FindPolynomialZerosd ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                   0, nt+1,  t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point3d) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = 0.5*(t[i]+t[i+1]);
      mbs_BCHornerC3d ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = DotProduct3d ( p, &cplanes[j] );
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point3d) );
        mbs_DivideBC3d ( degree, t[i], p, q );
        mbs_DivideBC3d ( degree, (t[i+1]-t[i])/(1.0-t[i]), p, q );
        output ( degree, q );
      }
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_ClipBC2Rd*/

boolean mbs_ClipBC3d ( int ncplanes, const vector4d *cplanes,
                       int degree, const point3d *cpoints,
                       void (*output) (int degree, const point3d *cpoints) )
{
  void    *sp;
  point3d *p, *q;
  double   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemd ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemd ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = DotProduct3d ( &cpoints[i], (point3d*)&cplanes[j] ) +
               cplanes[j].w;

      if ( !mbs_FindPolynomialZerosd ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                   0, nt+1,  t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point3d) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = 0.5*(t[i]+t[i+1]);
      mbs_BCHornerC3d ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = DotProduct3d ( p, (vector3d*)&cplanes[j] ) + cplanes[j].w;
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point3d) );
        mbs_DivideBC3d ( degree, t[i], p, q );
        mbs_DivideBC3d ( degree, (t[i+1]-t[i])/(1.0-t[i]), p, q );
        output ( degree, q );
      }
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_ClipBC3d*/

boolean mbs_ClipBC3Rd ( int ncplanes, const vector4d *cplanes,
                        int degree, const point4d *cpoints,
                        void (*output) (int degree, const point4d *cpoints) )
{
  void     *sp;
  point4d  *p, *q;
  double   *a, *t, x;
  int      i, j, nt, nz;
  boolean  visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemd ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemd ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = DotProduct4d ( &cpoints[i], &cplanes[j] );

      if ( !mbs_FindPolynomialZerosd ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double),
                   0, nt+1,  t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point4d) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = 0.5*(t[i]+t[i+1]);
      mbs_BCHornerC4d ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = DotProduct4d ( p, &cplanes[j] );
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point4d) );
        mbs_DivideBC4d ( degree, t[i], p, q );
        mbs_DivideBC4d ( degree, (t[i+1]-t[i])/(1.0-t[i]), p, q );
        output ( degree, q );
      }
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_ClipBC3Rd*/

