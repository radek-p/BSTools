
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
static int nsigns ( int n, const float *a )
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

static float *_coeff;
static int _deg;

static float poly ( float x )
{
  float f;

  mbs_BCHornerC1f ( _deg, _coeff, x, &f );
  return f;
} /*poly*/

boolean mbs_FindPolynomialZerosf ( int degree, const float *coeff,
                                   int *nzeros, float *zeros, float eps )
{
  void    *sp;
  float   *a, *b, x;
  int     nz, stp;
  boolean error;

  sp = pkv_GetScratchMemTop ();

        /* push the interval [0,1] and the related polynomial coefficients */
  if ( !(a = pkv_GetScratchMemf ( degree+3 )) )
    goto failure;
  _deg = degree;
  memcpy ( &a[2], coeff, (degree+1)*sizeof(float) );
  a[0] = 0.0;
  a[1] = 1.0;
  stp = 1;

  nz = 0;
  do {
    stp --;
    switch ( nsigns ( degree+1, &a[2] ) ) {
  case 0:    /* no sign change - reject the interval */
      pkv_FreeScratchMemf ( degree+3 );
      a -= degree+3;
      break;

  case 1:    /* one sign change - try to solve */
      _coeff = &a[2];
      x = pkn_Illinoisf ( poly, 0.0, 1.0, eps, &error );
      if ( error )
        goto divide;
      zeros[nz++] = (float)((1.0-x)*a[0] + x*a[1]);
      if ( nz > degree )
        goto failure;
      pkv_FreeScratchMemf ( degree+3 );
      a -= degree+3;
      break;

  default:   /* more than one sign change - bisect and push */
divide:
      if ( !(b = pkv_GetScratchMemf ( degree+3 )) )
        goto failure;
      stp += 2;
      b[0] = a[0];
      a[0] = b[1] = (float)(0.5*(b[0]+a[1]));
      mbs_BisectBC1f ( degree, &a[2], &b[2] );
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
} /*mbs_FindPolynomialZerosf*/


/* /////////////////////////////////////////// */
#define EPS 1.0e-6

boolean mbs_ClipBC2f ( int ncplanes, const vector3f *cplanes,
                       int degree, const point2f *cpoints,
                       void (*output) (int degree, const point2f *cpoints) )
{
  void    *sp;
  point2f *p, *q;
  float   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemf ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemf ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = (float)DotProduct2f ( &cpoints[i], (point2f*)&cplanes[j] ) +
                        cplanes[j].z;

      if ( !mbs_FindPolynomialZerosf ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                   0, nt+1, t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point2f) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = (float)(0.5*(t[i]+t[i+1]));
      mbs_BCHornerC2f ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = (float)DotProduct2f ( p, (vector2f*)&cplanes[j] ) + cplanes[j].z;
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point2f) );
        mbs_DivideBC2f ( degree, t[i], p, q );
        mbs_DivideBC2f ( degree, (float)((t[i+1]-t[i])/(1.0-t[i])), p, q );
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
} /*mbs_ClipBC2f*/

boolean mbs_ClipBC2Rf ( int ncplanes, const vector3f *cplanes,
                        int degree, const point3f *cpoints,
                        void (*output) (int degree, const point3f *cpoints) )
{
  void    *sp;
  point3f *p, *q;
  float   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemf ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemf ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = (float)DotProduct3f ( &cpoints[i], &cplanes[j] );

      if ( !mbs_FindPolynomialZerosf ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                   0, nt+1, t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point3f) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = (float)(0.5*(t[i]+t[i+1]));
      mbs_BCHornerC3f ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = (float)DotProduct3f ( p, &cplanes[j] );
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point3f) );
        mbs_DivideBC3f ( degree, t[i], p, q );
        mbs_DivideBC3f ( degree, (float)((t[i+1]-t[i])/(1.0-t[i])), p, q );
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
} /*mbs_ClipBC2Rf*/

boolean mbs_ClipBC3f ( int ncplanes, const vector4f *cplanes,
                       int degree, const point3f *cpoints,
                       void (*output) (int degree, const point3f *cpoints) )
{
  void    *sp;
  point3f *p, *q;
  float   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemf ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemf ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = (float)DotProduct3f ( &cpoints[i], (point3f*)&cplanes[j] ) +
                        cplanes[j].w;

      if ( !mbs_FindPolynomialZerosf ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                   0, nt+1, t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point3f) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = (float)(0.5*(t[i]+t[i+1]));
      mbs_BCHornerC3f ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = (float)DotProduct3f ( p, (vector3f*)&cplanes[j] ) + cplanes[j].w;
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point3f) );
        mbs_DivideBC3f ( degree, t[i], p, q );
        mbs_DivideBC3f ( degree, (float)((t[i+1]-t[i])/(1.0-t[i])), p, q );
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
} /*mbs_ClipBC3f*/

boolean mbs_ClipBC3Rf ( int ncplanes, const vector4f *cplanes,
                        int degree, const point4f *cpoints,
                        void (*output) (int degree, const point4f *cpoints) )
{
  void    *sp;
  point4f *p, *q;
  float   *a, *t, x;
  int     i, j, nt, nz;
  boolean visible;

  sp = pkv_GetScratchMemTop ();

  t = pkv_GetScratchMemf ( ncplanes*degree + 2 );
  if ( t ) {
    t[0] = 0.0;
    nt = 1;

    for ( j = 0; j < ncplanes; j++ ) {
        /* compute the polynomial coefficients */
      if ( !(a = pkv_GetScratchMemf ( degree+1)) )
        goto failure;
      for ( i = 0; i <= degree; i++ )
        a[i] = (float)DotProduct4f ( &cpoints[i], &cplanes[j] );

      if ( !mbs_FindPolynomialZerosf ( degree, a, &nz, &t[nt], EPS ) )
        goto failure;
      nt += nz;
    }
    t[nt] = 1.0;
    pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float),
                   0, nt+1, t );
/*
printf ( "nt = %d\n", nt );
for ( i = 0; i <= nt; i++ )
  printf ( "%f\n", t[i] );
*/
        /* now output the arcs */
    p = pkv_GetScratchMem ( 2*(degree+1)*sizeof(point4f) );
    if ( !p )
      goto failure;
    q = &p[degree+1];
    for ( i = 0; i < nt; i++ ) {
          /* is it visible? */
      x = (float)(0.5*(t[i]+t[i+1]));
      mbs_BCHornerC4f ( degree, cpoints, x, p );
      for ( j = 0, visible = true;  j < ncplanes && visible;  j++ ) {
        x = (float)DotProduct4f ( p, &cplanes[j] );
        visible = (boolean)(x > 0.0);
      }
      if ( visible ) {
        memcpy ( p, cpoints, (degree+1)*sizeof(point4f) );
        mbs_DivideBC4f ( degree, t[i], p, q );
        mbs_DivideBC4f ( degree, (float)((t[i+1]-t[i])/(1.0-t[i])), p, q );
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
} /*mbs_ClipBC3Rf*/

