
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
int mbs_multideBoorDerd ( int degree, int lastknot,
                          const double *knots,
                          int ncurves, int spdimen, int pitch,
                          const double *ctlpoints,
                          double t,
                          double *cpoints, double *dervect )
{
  void *sp;
  int i, j, k, r, l, m, ll, ik;
  int dpitch;
  double *d;
  double alpha, beta, factor;

  if ( degree < 1 ) {
    memset ( dervect, 0, ncurves*spdimen*sizeof(double) );
    return mbs_multideBoord ( degree, lastknot, knots, ncurves, spdimen,
                              pitch, ctlpoints, t, cpoints );
  }

  /* find the proper interval between the knots and the multiplicity of the
     left-side knot. Knots are numbered from 0 to lastknot ! */

  sp = pkv_GetScratchMemTop ();
  k = mbs_FindKnotIntervald ( degree, lastknot, knots, t, &r );

  if ( r < degree ) { /* a regular case */

    dpitch = (degree-r+1)*spdimen;
    if ( !( d = pkv_GetScratchMemd(ncurves*dpitch)) ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, 2, ERRMSG_2 );
      goto failure;
    }

    _mbs_multideBoorKerneld ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 1, dpitch, d );

  /* the last step of the de Boor algorithm; compute also the derivative  */
    factor = (double)degree/(knots[k+1]-knots[k-r]);
    alpha = (t-knots[k-r])/(knots[k+1]-knots[k-r]);
    beta = 1.0-alpha;
    for ( l = ll = k = 0;  l < ncurves;  l++, ll+= dpitch )
      for ( m = 0;  m < spdimen;  m++, k++ ) {
        dervect[k] = factor*(d[ll+spdimen+m]-d[ll+m]);
        d[ll+m] = beta*d[ll+m] + alpha*d[ll+spdimen+m];
      }
    pkv_Selectd ( ncurves, spdimen, dpitch, spdimen, &d[0], &cpoints[0] );

                                       /* cleanup and return */
  }
  else { /* a singular case, the knot multiplicity is not less than degree */
           /* the curve point is a control point, just copy data */
    pkv_Selectd ( ncurves, spdimen, pitch, spdimen,
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
        dervect[j] = factor*(ctlpoints[ik+spdimen+i]-ctlpoints[ik+i]);

  }
  pkv_SetScratchMemTop ( sp );
  return degree-r;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*mbs_multideBoorDerd*/

boolean mbs_deBoorDerPd ( int degreeu, int lastknotu, const double *knotsu,
                          int degreev, int lastknotv, const double *knotsv,
                          int spdimen, int pitch, const double *ctlpoints,
                          double u, double v,
                          double *ppoint, double *uder, double *vder )
{
  void   *sp;
  int    k, l;
  double *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemd ( 2*(degreeu+1)*spdimen )) )
    goto failure;

                                       /* find knot intervals */
  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );

                                       /* find arc of constant parameter v */
                                       /* and the arc of partial derivative */
  mbs_multideBoorDerd ( degreev, 2*degreev+1,
      &knotsv[l-degreev], degreeu+1, spdimen, pitch,
      &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
      v, q, &q[spdimen*(degreeu+1)] );

                                       /* find the point of this arc */
                                       /* and its derivative - it is */
                                       /* the partial derivative with */
                                       /* respect to u */
  mbs_multideBoorDerd ( degreeu, 2*degreeu+1,
                        &knotsu[k-degreeu], 1, spdimen, 0,
                        q, u, ppoint, uder );
  mbs_multideBoord ( degreeu, 2*degreeu+1,
                     &knotsu[k-degreeu], 1, spdimen, 0,
                     &q[spdimen*(degreeu+1)], u, vder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDerPd*/

