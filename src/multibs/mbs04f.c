
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

/* /////////////////////////////////////////// */
/* de Boor algorithm                           */
void _mbs_multideBoorKernelf ( int degree, const float *knots,
                                int ncurves, int spdimen,
                                int pitch, const float *ctlpoints,
                                float t, int k, int r, int lj,
                                int dpitch, float *d )
{
  int    i, j, m, ik, l, ll;
  double alpha, beta;

  /* copy control points for the deBoor algorithm */
  pkv_Selectf ( ncurves, dpitch, pitch, dpitch,
                &ctlpoints[(k-degree)*spdimen], &d[0] );

      /* de Boor algorithm without the last lj steps */
  for ( j = 1; j <= degree-r-lj; j++ )
    for ( i = k-degree+j, ik = 0;  i <= k-r;  i++, ik += spdimen ) {
      alpha = (t-knots[i])/(knots[i+degree+1-j]-knots[i]);
      beta = 1.0-alpha;
      for ( l = ll = 0;  l < ncurves;  l++, ll += dpitch )
        for ( m = 0;  m < spdimen;  m++ )
          d[ll+ik+m] = (float)(beta*d[ll+ik+m] + alpha*d[ll+ik+spdimen+m]);
    }
} /*_mbs_multideBoorKernelf*/

int mbs_multideBoorf ( int degree, int lastknot,
                       const float *knots,
                       int ncurves, int spdimen, int pitch,
                       const float *ctlpoints,
                       float t,
                       float *cpoints )
{
  void  *sp;
  int   k, r;
  int   dpitch;
  float *d;

  /* find the proper interval between the knots and the multiplicity of the
     left-side knot. Knots are numbered from 0 to lastknot ! */

  sp = pkv_GetScratchMemTop ();
  k = mbs_FindKnotIntervalf ( degree, lastknot, knots, t, &r );
  if ( r > degree )
    r = degree;

  /* copy control points for the deBoor algorithm */
  dpitch = (degree-r+1)*spdimen;
  if ( !( d = pkv_GetScratchMemf(ncurves*dpitch)) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, 2, ERRMSG_2 );
    goto failure;
  }
  _mbs_multideBoorKernelf ( degree, knots, ncurves, spdimen,
                             pitch, ctlpoints, t, k, r, 0, dpitch, d );
  /* copy the result and cleanup */
  pkv_Selectf ( ncurves, spdimen, dpitch, spdimen, &d[0], &cpoints[0] );
  pkv_SetScratchMemTop ( sp );
  return degree-r;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*mbs_multideBoorf*/

void mbs_deBoorC2Rf ( int degree, int lastknot,   
                      const float *knots, const point3f *ctlpoints,  
                      float t, point2f *cpoint )
{
  point3f p;

  mbs_multideBoorf ( degree, lastknot, knots, 1, 3, 0, (float*)ctlpoints, t,
                     &p.x );
  Point3to2f ( &p, cpoint );
} /*mbs_deBoorC2Rf*/

void mbs_deBoorC3Rf ( int degree, int lastknot,   
                      const float *knots, const point4f *ctlpoints,  
                      float t, point3f *cpoint )
{
  point4f p;

  mbs_multideBoorf ( degree, lastknot, knots, 1, 4, 0, (float*)ctlpoints, t,
                     &p.x );
  Point4to3f ( &p, cpoint );
} /*mbs_deBoorC3Rf*/

boolean mbs_deBoorPf ( int degreeu, int lastknotu, const float *knotsu,
                       int degreev, int lastknotv, const float *knotsv,
                       int spdimen, int pitch,
                       const float *ctlpoints,
                       float u, float v, float *ppoint )
{
  void   *sp;
  int    k, l;
  float  *q;

  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMemf ( (degreeu+1)*spdimen );
  if ( !q )
    return false;

                                       /* find knot intervals */
  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoorf ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, spdimen, pitch,
      &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen], v, q );

                                       /* find the point of this arc */
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1,
                     spdimen, 0, q, u, (float*)ppoint );

  pkv_SetScratchMemTop ( sp );
  return true;
} /*mbs_deBoorP3f*/

boolean mbs_deBoorP3f ( int degreeu, int lastknotu, const float *knotsu,
                        int degreev, int lastknotv, const float *knotsv,
                        int pitch,
                        const point3f *ctlpoints,
                        float u, float v, point3f *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  float  *q;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point3f) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoorf ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 3, pitch,
      (float*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 3, 0, q,
                     u, (float*)ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  return mbs_deBoorPf ( degreeu, lastknotu, knotsu,
                        degreev, lastknotv, knotsv,
                        3, pitch, (float*) ctlpoints, u, v, (float*)ppoint );
#endif
} /*mbs_deBoorP3f*/

boolean mbs_deBoorP3Rf ( int degreeu, int lastknotu, const float *knotsu,
                         int degreev, int lastknotv, const float *knotsv,
                         int pitch,
                         const point4f *ctlpoints,
                         float u, float v, point3f *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  float  *q;
  point4f r;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point4f) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoorf ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 4, pitch,
      (float*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 4, 0, q,
                     u, (float*)&r );
  Point4to3f ( &r, ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  point4f r;

  if ( mbs_deBoorPf ( degreeu, lastknotu, knotsu,
                      degreev, lastknotv, knotsv,
                      4, pitch, (float*) ctlpoints, u, v, (float*)&r ) {
    Point4To3f ( &r, ppoint );
    return true;
  }
  else
    return false;
#endif
} /*mbs_deBoorP3Rf*/

boolean mbs_deBoorP4f ( int degreeu, int lastknotu, const float *knotsu,
                        int degreev, int lastknotv, const float *knotsv,
                        int pitch,
                        const point4f *ctlpoints,
                        float u, float v, point4f *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  float  *q;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point4f) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoorf ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 4, pitch,
      (float*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 4, 0, q,
                     u, (float*)ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  return mbs_deBoorPf ( degreeu, lastknotu, knotsu,
                        degreev, lastknotv, knotsv,
                        4, pitch, (float*) ctlpoints, u, v, (float*)ppoint );
#endif
} /*mbs_deBoorP4f*/

