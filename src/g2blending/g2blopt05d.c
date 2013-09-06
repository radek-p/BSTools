
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"

#include "g2blprivated.h"

#define _DEBUG

/* ////////////////////////////////////////////////////////////////////////// */
double g2bl_ClosedSurfNetDiameterSqd ( int lastknotu, int lastknotv,
                                       int pitch, const point3d *cp )
{
  void     *sp;
  double   *ccp, *knotsv, *knotsa;
  point3d  *acp;
  double   d, dmax;
  vector3d v;
  int      i, j, k, pitch1;

  sp = pkv_GetScratchMemTop ();
        /* Find a clamped boundary representation */
  knotsv = pkv_GetScratchMemd ( lastknotv+1 );
  knotsa = pkv_GetScratchMemd ( 4 );
  ccp = pkv_GetScratchMemd ( (lastknotu-6)*(lastknotv-3)*3 );
  if ( !knotsv || !knotsa || !ccp ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  for ( i = 0; i <= lastknotv; i++ )
    knotsv[i] = (double)i;
  pitch1 = 3*(lastknotv-3);
  pkv_Selectd ( lastknotu-6, pitch1, pitch, pitch1, cp, ccp );

  knotsa[0] = 0.0;
  for ( i = 1; i < 4; i++ )
    knotsa[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( lastknotu-6, 3, 3, knotsv,
                                pitch1, ccp, knotsa );
  for ( i = 0; i < 3; i++ )
    knotsa[i] = (double)(lastknotv-3);
  knotsa[3] = lastknotv;
  mbs_multiBSChangeRightKnotsd ( lastknotu-6, 3, 3, lastknotv, knotsv,
                                 pitch1, ccp, knotsa );

        /* Now extract the control points, which determine */
        /* the boundary conditions */
  pkv_Moved ( lastknotu-6, 3, pitch1, 6-pitch1, &ccp[pitch1-3] );
  pkv_Rearranged ( lastknotu-6, 6, pitch1, 6, ccp );

        /* number of points */
  k = 2*(lastknotu-6);
  acp = (point3d*)ccp;

        /* find the maximal distance between the k points in the acp array*/
  dmax = 0.0;
  for ( i = 1; i < k; i++ )
    for ( j = 0; j < i; j++ ) {
      SubtractPoints3d ( &acp[i], &acp[j], &v );
      d = DotProduct3d ( &v, &v );
      dmax = max ( d, dmax );
    }
#ifdef _DEBUG
printf ( "dmax = %f\n", dmax );
#endif
  pkv_SetScratchMemTop ( sp );
  return dmax;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g2bl_ClosedSurfNetDiameterSqd*/

/* ////////////////////////////////////////////////////////////////////////// */
double _g2bl_ClosedAuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab )
{
  if ( _g2bl_ShiftDecompHessiand ( neqs, hsize, prof, hessian,
                                   Lhessian, Lhrows, nu ) ) {
    pkn_NRBLowerTrSolved ( neqs, prof, Lhessian, Lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( neqs, prof, Lhessian, Lhrows, 1, 1, dcoeff, 1, dcoeff );
    pkn_SubtractMatrixd ( lastknotu-6, 3*(lastknotv-9), pitch, &cp[3].x,
                          3*(lastknotv-9), dcoeff, 3*(lastknotv-3), &acp[3].x );
    memcpy ( &acp[(lastknotu-6)*(lastknotv-3)], acp,
             9*(lastknotv-3)*sizeof(double) );
#ifdef _DEBUG
printf ( "F" );
#endif
    return g2bl_UFuncd ( nknots, qcoeff, Nitab, lastknotu, lastknotv,
                         3*(lastknotv-3), acp, NULL, tC, ftab );
  }
  else
    return MYINFINITY;
} /*_g2bl_ClosedAuxNuFuncd*/

