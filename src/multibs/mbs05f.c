
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
/* de Boor algorithm for computing 1st derivative  */
int mbs_multideBoorDerf ( int degree, int lastknot,
                          const float *knots,
                          int ncurves, int spdimen, int pitch,
                          const float *ctlpoints,
                          float t,
                          float *cpoints, float *dervect )
{
  void *sp;
  int i, j, k, r, l, m, ll, ik;
  int dpitch;
  float *d;
  double alpha, beta, factor;

  if ( degree < 1 ) {
    memset ( dervect, 0, ncurves*spdimen*sizeof(float) );
    return mbs_multideBoorf ( degree, lastknot, knots, ncurves, spdimen,
                              pitch, ctlpoints, t, cpoints );
  }

  /* find the proper interval between the knots and the multiplicity of the
     left-side knot. Knots are numbered from 0 to lastknot ! */

  sp = pkv_GetScratchMemTop ();
  k = mbs_FindKnotIntervalf ( degree, lastknot, knots, t, &r );

  if ( r < degree ) { /* a regular case */

    dpitch = (degree-r+1)*spdimen;
    if ( !( d = pkv_GetScratchMemf(ncurves*dpitch)) ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, 2, ERRMSG_2 );
      exit ( 1 );
    }

    _mbs_multideBoorKernelf ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 1, dpitch, d );

  /* the last step of the de Boor algorithm; compute also the derivative  */
    factor = (float)degree/(knots[k+1]-knots[k-r]);
    alpha = (t-knots[k-r])/(knots[k+1]-knots[k-r]);
    beta = 1.0-alpha;
    for ( l = ll = k = 0;  l < ncurves;  l++, ll+= dpitch )
      for ( m = 0;  m < spdimen;  m++, k++ ) {
        dervect[k] = (float)(factor*(d[ll+spdimen+m]-d[ll+m]));
        d[ll+m] = (float)(beta*d[ll+m] + alpha*d[ll+spdimen+m]);
      }
    pkv_Selectf ( ncurves, spdimen, dpitch, spdimen, &d[0], &cpoints[0] );

                                       /* cleanup and return */
  }
  else { /* a singular case, the knot multiplicity is not less than degree */
           /* the curve point is a control point, just copy data */
    pkv_Selectf ( ncurves, spdimen, pitch, spdimen,
                  &ctlpoints[(k-degree)*spdimen], cpoints );

           /* unless t is at the end of the curve domain, compute the */
           /* right-sided derivative */
    if ( k >= lastknot-degree ) {
      k = lastknot-2;
      factor = degree/(knots[lastknot-degree]-knots[lastknot-degree-1]);
    }
    else
      factor = degree/(knots[k+1]-knots[k]);

    for ( l = 0,  j = 0,  ik = (k-degree)*spdimen;
          l < ncurves;
          l++,  ik += pitch )
      for ( i = 0; i < spdimen; i++, j++ ) 
        dervect[j] = (float)factor*(ctlpoints[ik+spdimen+i]-ctlpoints[ik+i]);

  }
  pkv_SetScratchMemTop ( sp );
  return degree-r;
} /*mbs_multideBoorDerf*/

boolean mbs_deBoorDerPf ( int degreeu, int lastknotu, const float *knotsu,
                          int degreev, int lastknotv, const float *knotsv,
                          int spdimen, int pitch, const float *ctlpoints,
                          float u, float v,
                          float *ppoint, float *uder, float *vder )
{
  void   *sp;
  int    k, l;
  float  *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemf ( 2*(degreeu+1)*spdimen )) )
    goto failure;

                                       /* find knot intervals */
  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
                                       /* and the arc of partial derivative */
  mbs_multideBoorDerf ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, spdimen, pitch,
      &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
      v, q, &q[spdimen*(degreeu+1)] );

                                       /* find the point of this arc */
                                       /* and its derivative - it is */
                                       /* the partial derivative with */
                                       /* respect to u */
  mbs_multideBoorDerf ( degreeu, 2*degreeu+1,
                        &knotsu[k-degreeu], 1, spdimen, 0,
                        q, u, ppoint, uder );
  mbs_multideBoorf ( degreeu, 2*degreeu+1,
                     &knotsu[k-degreeu], 1, spdimen, 0,
                     &q[spdimen*(degreeu+1)], u, vder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDerPf*/

