
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
double g2bl_SurfNetDiameterSqd ( int lastknotu, int lastknotv,
                                 int pitch, const point3d *cp )
{
  void     *sp;
  double   *ccp, *knotsu, *knotsv, *knotsa;
  point3d  *acp;
  double   d, dmax;
  vector3d v;
  int      i, j, k, pitch1;

  sp = pkv_GetScratchMemTop ();
        /* Find a clamped boundary representation */
  knotsu = pkv_GetScratchMemd ( lastknotu+1 );
  knotsv = pkv_GetScratchMemd ( lastknotv+1 );
  knotsa = pkv_GetScratchMemd ( 4 );
  ccp = pkv_GetScratchMemd ( (lastknotu-3)*(lastknotv-3)*3 );
  if ( !knotsu || !knotsv || !knotsa || !ccp ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  for ( i = 0; i <= lastknotu; i++ )
    knotsu[i] = (double)i;
  for ( i = 0; i <= lastknotv; i++ )
    knotsv[i] = (double)i;
  pitch1 = 3*(lastknotv-3);
  pkv_Selectd ( lastknotu-3, pitch1, pitch, pitch1, cp, ccp );

  knotsa[0] = 0.0;
  for ( i = 1; i < 4; i++ )
    knotsa[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( lastknotu-3, 3, 3, knotsv,
                                pitch1, ccp, knotsa );
  for ( i = 0; i < 3; i++ )
    knotsa[i] = (double)(lastknotv-3);
  knotsa[3] = lastknotv;
  mbs_multiBSChangeRightKnotsd ( lastknotu-3, 3, 3, lastknotv, knotsv,
                                 pitch1, ccp, knotsa );
  knotsa[0] = 0.0;
  for ( i = 1; i < 4; i++ )
    knotsa[i] = 3.0;
  mbs_multiBSChangeLeftKnotsd ( 1, pitch1, 3, knotsu,
                                0, ccp, knotsa );
  for ( i = 0; i < 3; i++ )
    knotsa[i] = (double)(lastknotu-3);
  knotsa[3] = lastknotu;
  mbs_multiBSChangeRightKnotsd ( 1, pitch1, 3, lastknotu, knotsu,
                                 0, ccp, knotsa );

        /* Now extract the control points, which determine */
        /* the boundary conditions */
  pkv_Moved ( lastknotu-9, 9, pitch1, -3*(lastknotv-9), &ccp[4*pitch1-9] );
  pkv_Rearranged ( lastknotu-9, 18, pitch1, 18, &ccp[3*pitch1] );
  memmove ( &ccp[3*pitch1+18*(lastknotu-9)], &ccp[(lastknotu-6)*pitch1],
            3*pitch1*sizeof(double) );
        /* number of points */
  k = 6*((lastknotv-3)+(lastknotu-9));
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
} /*g2bl_SurfNetDiameterSqd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean _g2bl_ComputeDeltaQd ( int n, const int *prof, double **hrows,
                               const double *grad, const double *dcoeff,
                               double *dq )
{
  void   *sp;
  double *aux;

  sp = pkv_GetScratchMemTop ();
  aux = pkv_GetScratchMemd ( 2*n );
  if ( !aux ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pkn_NRBSymMultd ( n, prof, hrows[0], hrows, 1, 1, dcoeff, 1, aux );
  pkn_MatrixLinCombd ( 1, n, 0, aux, 0.5, 0, grad, -1.0, 0, aux );
  *dq = pkn_ScalarProductd ( n, aux, dcoeff );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2bl_ComputeDeltaQd*/

boolean _g2bl_ShiftDecompHessiand ( int neqs, int hsize, int *prof, double *hessian,
                                    double *Lhessian, double **Lhrows, double nu )
{
  int i;

  memcpy ( Lhessian, hessian, hsize*sizeof(double) );
  for ( i = 0; i < neqs; i++ )
    Lhrows[i][i] += nu;
  return pkn_NRBSymCholeskyDecompd ( neqs, prof, Lhessian, Lhrows, NULL );
} /*_g2bl_ShiftDecompHessiand*/

double _g2bl_AuxNuFuncd ( int nknots,
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
    pkn_SubtractMatrixd ( lastknotu-9, 3*(lastknotv-9), pitch, &cp[pitch+3].x,
                          3*(lastknotv-9), dcoeff,
                          3*(lastknotv-3), &acp[3*(lastknotv-3)+3].x );
#ifdef _DEBUG
printf ( "F" );
#endif
    return g2bl_UFuncd ( nknots, qcoeff, Nitab, lastknotu, lastknotv,
                         3*(lastknotv-3), acp, NULL, tC, ftab );
  }
  else
    return MYINFINITY;
} /*_g2bl_AuxNuFuncd*/

