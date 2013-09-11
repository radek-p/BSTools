
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */
#include <stdlib.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* planar Bezier and rational Bezier curve rasterization */

#define T0 0.5
#define T1 0.53125
#define T2 0.46875

static char xround, yround;

static void _mbs_RoundPixelf ( point3f *p, xpoint *px )
{
  float x, fx, y, fy, t;

  x = p->x/p->z;  fx = (float)floor ( x );
  y = p->y/p->z;  fy = (float)floor ( y );
  switch ( xround ) {
 case 0: t = T0;  break;
 case 1: t = T1;  break;
 case 2: t = T2;  break;
default: t = T0;
  }
  if ( x-fx < t ) {
    px->x = (short)fx;      xround = 1;
  }
  else {
    px->x = (short)(fx+1);  xround = 2;
  }
  switch ( yround ) {
case 0: t = T0;  break;
case 1: t = T1;  break;
case 2: t = T2;  break;
  }
  if ( y-fy < t ) {
    px->y = (short)fy;      yround = 1;
  }
  else {
    px->y = (short)(fy+1);  yround = 2;
  }
} /*_mbs_RoundPixelf*/

static boolean _mbs_IsOnePixel ( int degree, vector3f *cp )
{
  xpoint a, b;
  int    i;

  _mbs_RoundPixelf ( &cp[0], &a );
  for ( i = 1; i <= degree; i++ ) {
    _mbs_RoundPixelf ( &cp[i], &b );
    if ( a.x != b.x || a.y != b.y )
      return false;
  }
  return true;
} /*_mbs_IsOnePixel*/

static boolean _mbs_RasterizeBCf ( int degree, vector3f *acp )
{
#define STKSIZE  32

  typedef struct {
    xpoint px;
    float  vx;
  } stackel;

  void     *sp;
  xpoint   last, next;
  float    v0, v1;
  stackel  *stack;
  int      stp;
  vector3f p;
  int      i;

#define NEIGHBOURS(a,b) \
  abs(a.x-b.x) <= 1 && abs(a.y-b.y) <= 1
#define OUTPUT(px) \
  { if ( _pkv_npix == PKV_BUFSIZE ) { _pkv_npix -= 2;  PKV_FLUSH \
      memcpy ( &_pkv_pixbuf[0], &_pkv_pixbuf[PKV_BUFSIZE-2], 2*sizeof(xpoint) ); \
      _pkv_npix = 2; } \
    if ( _pkv_npix > 1 ) \
      { if ( NEIGHBOURS ( px, _pkv_pixbuf[_pkv_npix-2] ) ) _pkv_npix--; } \
    _pkv_pixbuf[_pkv_npix] = px;  _pkv_npix++; }
#define PUSH(p,v) \
  { if (stp < STKSIZE) { stack[stp].px = p;  stack[stp].vx = v;  stp++; } }
#define READ(p,v) \
  { p = stack[stp-1].px;  v = stack[stp-1].vx; }
#define POP \
  stp--;

        /* allocate stack */
  sp = pkv_GetScratchMemTop ();
  stack = (stackel*)pkv_GetScratchMem ( STKSIZE*sizeof(stackel) );
  if ( !stack ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

        /* convert control points to the scaled basis */
  mbs_multiBezScalef ( degree, 1, 1, 3, 0, (float*)acp );

        /* draw the first half of the arc */
  _mbs_RoundPixelf ( &acp[0], &last );  v0 = 0.0;
  OUTPUT ( last )
  stp = 0;         /* stack becomes empty */
  p = acp[0];
  for ( i = 1; i <= degree; i++ )
    AddVector3f ( &p, &acp[i], &p );
  _mbs_RoundPixelf ( &p, &next );
  PUSH ( next, 1.0 )
  do {
    READ ( next, v1 )
    if ( NEIGHBOURS ( last, next ) ) {
      POP
      last = next;  v0 = v1;
      OUTPUT ( last )
    }
    else {
      v1 = (float)(0.5*(v0+v1));
      p = acp[degree];     /* Horner scheme for the power basis */
      for ( i = degree-1; i >= 0; i-- )
        AddVector3Mf ( &acp[i], &p, v1, &p );
      _mbs_RoundPixelf ( &p, &next );
      PUSH ( next, v1 )
    }
  } while ( stp );

        /* draw the second half of the arc */
  _mbs_RoundPixelf ( &acp[degree], &next );
  PUSH ( next, 0.0 )
  do {
    READ ( next, v1 );
    if ( NEIGHBOURS ( last, next ) ) {
      POP
      last = next;  v0 = v1;
      OUTPUT ( last )
    }
    else {
      v1 = (float)(0.5*(v0+v1));
      p = acp[0];          /* Horner scheme for the power basis */
      for ( i = 1; i <= degree; i++ )
        AddVector3Mf ( &acp[i], &p, v1, &p );
      _mbs_RoundPixelf ( &p, &next );
      PUSH ( next, v1 )
    }
  } while ( stp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef STKSIZE
#undef NEIGHBOURS
#undef OUTPUT
#undef PUSH
#undef READ
#undef POP
} /*_mbs_RasterizeBCf*/

static boolean _mbs_RasterizeBC2Rf ( int degree, const point3f *cpoints )
{
  void     *sp;
  char     *cp, *cq;
  vector3f *acp;
  int      size, stp;
  vector2f v, w;

  if ( degree == 1 ) {
    _pkv_DrawLine ( (int)(cpoints[0].x/cpoints[0].z+0.5),
                    (int)(cpoints[0].y/cpoints[0].z+0.5),
                    (int)(cpoints[1].x/cpoints[1].z+0.5),
                    (int)(cpoints[1].y/cpoints[1].z+0.5) );
    return true;
  }

  sp = pkv_GetScratchMemTop ();
  acp = (vector3f*)pkv_GetScratchMem ( (degree+1)*sizeof(vector3f) );
  cp = pkv_GetScratchMem ( size = (degree+1)*sizeof(point3f) );
  if ( !acp || !cp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

  memcpy ( cp, cpoints, size );  cp += size;
  stp = 1;
  do {
    stp--;  cp -= size;
    if ( _mbs_IsOnePixel ( degree, (point3f*)cp ) )
      goto draw_it;
    if ( ((point3f*)cp)[0].z )
      Point3to2f ( &((point3f*)cp)[0], &v );
    else {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_7, ERRMSG_7 );
      goto failure;
    }
    if ( ((point3f*)cp)[degree].z )
      Point3to2f ( &((point3f*)cp)[degree], &w );
    else {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_7, ERRMSG_7 );
      goto failure;
    }
    SubtractPoints2f ( &w, &v, &v );
    if ( mbs_MonotonicPolylineRf ( 3, degree+1, 3, (float*)cp, (float*)&v ) ) {
draw_it:
      memcpy ( acp, cpoints, (degree+1)*sizeof(vector3f) );
      _mbs_RasterizeBCf ( degree, acp );
      pkv_FreeScratchMem ( size );
    }
    else {
      cq = pkv_GetScratchMem ( size );
      if ( !cq ) {
        PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
        goto failure;
      }
      mbs_BisectBC3f ( degree, cp, cq );
      stp += 2;  cp = cq+size;
    }
  } while ( stp );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_RasterizeBC2Rf*/

static boolean _mbs_RasterizeBC2f ( int degree, const point2f *cpoints )
{
  void    *sp;
  point3f *cp;
  int     i;

  sp = pkv_GetScratchMemTop ();
  cp = pkv_GetScratchMem ( (degree+1)*sizeof(point3f) );
  if ( !cp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pkv_Selectf ( degree+1, 2, 2, 3, (float*)cpoints, (float*)cp );
  for ( i = 0; i <= degree; i++ )
    cp[i].z = 1.0;
  _mbs_RasterizeBC2Rf ( degree, cp );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_RasterizeBC2f*/

boolean mbs_RasterizeBC2f ( int degree, const point2f *cpoints,
                            void (*output)(const xpoint *buf, int n),
                            boolean outlast )
{
  void *sp;

  sp = pkv_GetScratchMemTop ();
  xround = yround = 0;
  _pkv_InitPixelBuffer ();
  if ( _pkv_pixbuf ) {
    _pkv_OutputPixels = output;
    if ( !_mbs_RasterizeBC2f ( degree, cpoints ) )
      goto failure;
    if ( _pkv_npix && !outlast )
      _pkv_npix--;
    PKV_FLUSH
  }
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_RasterizeBC2f*/

boolean mbs_RasterizeBC2Rf ( int degree, const point3f *cpoints,
                          void (*output)(const xpoint *buf, int n),
                          boolean outlast )
{
  void *sp;

  sp = pkv_GetScratchMemTop ();
  xround = yround = 0;
  _pkv_InitPixelBuffer ();
  if ( !_pkv_pixbuf )
    return false;
  _pkv_OutputPixels = output;
  if ( !_mbs_RasterizeBC2Rf ( degree, cpoints ) )
    goto failure;
  if ( _pkv_npix && !outlast )
    _pkv_npix--;
  PKV_FLUSH
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
 _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_RasterizeBC2Rf*/

boolean mbs_RasterizeBS2f ( int degree, int lastknot, const float *knots,
                            const point2f *cpoints,
                            void (*output)(const xpoint *buf, int n),
                            boolean outlast )
{
  void    *sp;
  int     ku, i;
  point2f *cp;

  sp = pkv_GetScratchMemTop ();
  xround = yround = 0;
  _pkv_InitPixelBuffer ();
  if ( !_pkv_pixbuf )
    return false;
  ku = mbs_NumKnotIntervalsf ( degree, lastknot, knots );
  cp = (point2f*)pkv_GetScratchMem ( ku*(degree+1)*sizeof(point2f) );
  if ( !cp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  mbs_BSToBezC2f ( degree, lastknot, knots, cpoints, &ku, NULL, NULL, cp );

  _pkv_OutputPixels = output;
  for ( i = 0; i < ku; i++, cp += degree+1 )
    if ( !_mbs_RasterizeBC2f ( degree, cp ) )
      goto failure;
  if ( _pkv_npix && !outlast )
    _pkv_npix--;
  PKV_FLUSH
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_RasterizeBS2f*/

boolean mbs_RasterizeBS2Rf ( int degree, int lastknot, const float *knots,
                             const point3f *cpoints,
                             void (*output)(const xpoint *buf, int n),
                             boolean outlast )
{
  void    *sp;
  int     ku, i;
  point3f *cp;

  sp = pkv_GetScratchMemTop ();
  xround = yround = 0;
  _pkv_InitPixelBuffer ();
  if ( !_pkv_pixbuf )
    return false;
  ku = mbs_NumKnotIntervalsf ( degree, lastknot, knots );
  cp = (point3f*)pkv_GetScratchMem ( ku*(degree+1)*sizeof(point3f) );
  if ( !cp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  mbs_BSToBezC3f ( degree, lastknot, knots, cpoints, &ku, NULL, NULL, cp );

  _pkv_OutputPixels = output;
  for ( i = 0; i < ku; i++, cp += degree+1 )
    if ( !_mbs_RasterizeBC2Rf ( degree, cp ) )
      goto failure;
  if ( _pkv_npix && !outlast )
    _pkv_npix--;
  PKV_FLUSH
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_RasterizeBS2Rf*/

