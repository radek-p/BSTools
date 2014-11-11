
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
void _mbs_multideBoorKerneld ( int degree, const double *knots,
                                int ncurves, int spdimen,
                                int pitch, const double *ctlpoints,
                                double t, int k, int r, int lj,
                                int dpitch, double *d )
{
  int    i, j, m, ik, l, ll;
  double alpha, beta;

  /* copy control points for the deBoor algorithm */
  pkv_Selectd ( ncurves, dpitch, pitch, dpitch,
                &ctlpoints[(k-degree)*spdimen], &d[0] );

      /* de Boor algorithm without the last lj steps */
  for ( j = 1; j <= degree-r-lj; j++ )
    for ( i = k-degree+j, ik = 0;  i <= k-r;  i++, ik += spdimen ) {
      alpha = (t-knots[i])/(knots[i+degree+1-j]-knots[i]);
      beta = 1.0-alpha;
      for ( l = ll = 0;  l < ncurves;  l++, ll += dpitch )
        for ( m = 0;  m < spdimen;  m++ )
          d[ll+ik+m] = beta*d[ll+ik+m] + alpha*d[ll+ik+spdimen+m];
    }
} /*_mbs_multideBoorKerneld*/

int mbs_multideBoord ( int degree, int lastknot,
                       const double *knots,
                       int ncurves, int spdimen, int pitch,
                       const double *ctlpoints,
                       double t,
                       double *cpoints )
{
  void  *sp;
  int   k, r;
  int   dpitch;
  double *d;

  /* find the proper interval between the knots and the multiplicity of the
     left-side knot. Knots are numbered from 0 to lastknot ! */

  sp = pkv_GetScratchMemTop ();
  k = mbs_FindKnotIntervald ( degree, lastknot, knots, t, &r );
  if ( r > degree )
    r = degree;

  /* copy control points for the deBoor algorithm */
  dpitch = (degree-r+1)*spdimen;
  if ( !( d = pkv_GetScratchMemd(ncurves*dpitch)) ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, 2, ERRMSG_2 );
    goto failure;
  }
  _mbs_multideBoorKerneld ( degree, knots, ncurves, spdimen,
                             pitch, ctlpoints, t, k, r, 0, dpitch, d );
  /* copy the result and cleanup */
  pkv_Selectd ( ncurves, spdimen, dpitch, spdimen, &d[0], &cpoints[0] );
  pkv_SetScratchMemTop ( sp );
  return degree-r;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*mbs_multideBoord*/

void mbs_deBoorC2Rd ( int degree, int lastknot,   
                      const double *knots, const point3d *ctlpoints,  
                      double t, point2d *cpoint )
{
  point3d p;

  mbs_multideBoord ( degree, lastknot, knots, 1, 3, 0, (double*)ctlpoints, t,
                     &p.x );
  Point3to2d ( &p, cpoint );
} /*mbs_deBoorC2Rd*/

void mbs_deBoorC3Rd ( int degree, int lastknot,   
                      const double *knots, const point4d *ctlpoints,  
                      double t, point3d *cpoint )
{
  point4d p;

  mbs_multideBoord ( degree, lastknot, knots, 1, 4, 0, (double*)ctlpoints, t,
                     &p.x );
  Point4to3d ( &p, cpoint );
} /*mbs_deBoorC3Rd*/

boolean mbs_deBoorPd ( int degreeu, int lastknotu, const double *knotsu,
                       int degreev, int lastknotv, const double *knotsv,
                       int spdimen, int pitch,
                       const double *ctlpoints,
                       double u, double v, double *ppoint )
{
  void   *sp;
  int    k, l;
  double *q;

  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMemd ( (degreeu+1)*spdimen );
  if ( !q )
    return false;

                                       /* find knot intervals */
  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoord ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, spdimen, pitch,
      &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen], v, q );

                                       /* find the point of this arc */
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1,
                     spdimen, 0, q, u, (double*)ppoint );

  pkv_SetScratchMemTop ( sp );
  return true;
} /*mbs_deBoorP3d*/

boolean mbs_deBoorP3d ( int degreeu, int lastknotu, const double *knotsu,
                        int degreev, int lastknotv, const double *knotsv,
                        int pitch,
                        const point3d *ctlpoints,
                        double u, double v, point3d *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  double *q;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point3d) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoord ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 3, pitch,
      (double*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 3, 0, q,
                     u, (double*)ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  return mbs_deBoorPd ( degreeu, lastknotu, knotsu,
                        degreev, lastknotv, knotsv,
                        3, pitch, (double*) ctlpoints, u, v, (double*)ppoint );
#endif
} /*mbs_deBoorP3d*/

boolean mbs_deBoorP3Rd ( int degreeu, int lastknotu, const double *knotsu,
                         int degreev, int lastknotv, const double *knotsv,
                         int pitch,
                         const point4d *ctlpoints,
                         double u, double v, point3d *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  double *q;
  point4d r;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point4d) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoord ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 4, pitch,
      (double*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 4, 0, q,
                     u, (double*)&r );
  Point4to3d ( &r, ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  point4d r;

  if ( mbs_deBoorPd ( degreeu, lastknotu, knotsu,
                      degreev, lastknotv, knotsv,
                      4, pitch, (double*) ctlpoints, u, v, (double*)&r ) {
    Point4To3d ( &r, ppoint );
    return true;
  }
  else
    return false;
#endif
} /*mbs_deBoorP3Rd*/

boolean mbs_deBoorP4d ( int degreeu, int lastknotu, const double *knotsu,
                        int degreev, int lastknotv, const double *knotsv,
                        int pitch,
                        const point4d *ctlpoints,
                        double u, double v, point4d *ppoint )
{
#ifndef OLD
  int    k, l;
  int    size_q;
  double *q;

  q = pkv_GetScratchMem ( size_q = (degreeu+1)*sizeof(point4d) );

                                       /* find knot intervals */
  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
  mbs_multideBoord ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, 4, pitch,
      (double*)&ctlpoints[(k-degreeu)*(lastknotv-degreev)+(l-degreev)], v, q );

                                       /* find the point of this arc */
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, 4, 0, q,
                     u, (double*)ppoint );

  pkv_FreeScratchMem ( size_q );
  return true;
#else
  return mbs_deBoorPd ( degreeu, lastknotu, knotsu,
                        degreev, lastknotv, knotsv,
                        4, pitch, (double*) ctlpoints, u, v, (double*)ppoint );
#endif
} /*mbs_deBoorP4d*/

