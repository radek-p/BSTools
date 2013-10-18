
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */
/* and modified by Przemyslaw Kiciak                                         */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"

#include "g1blprivated.h"

#define _DEBUG

/* ////////////////////////////////////////////////////////////////////////// */
double g1bl_SurfNetDiameterSqd ( int lastknotu, int lastknotv,
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
  knotsa = pkv_GetScratchMemd ( 3 );
  ccp = pkv_GetScratchMemd ( (lastknotu - DEG)*(lastknotv - DEG)*3 );
  if ( !knotsu || !knotsv || !knotsa || !ccp ) {
    PKV_SIGNALERROR ( LIB_G1BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  for ( i = 0; i <= lastknotu; i++ )
    knotsu[i] = (double)i;
  for ( i = 0; i <= lastknotv; i++ )
    knotsv[i] = (double)i;
  pitch1 = 3*(lastknotv- DEG);
  pkv_Selectd ( lastknotu-DEG, pitch1, pitch, pitch1, cp, ccp );

  /*  przejscie do reprezentacji dziedziny: 0, 2,2, ..... , podwojne wezly*/
  knotsa[0] = 0.0;
  for ( i = 1; i < 3; i++ )
    knotsa[i] = 2.0;
  /* zamiana punktow kontrolnych na brzegu tak by TYLKO pierwsza kolumna reprezentowala brzeg */
  if ( !mbs_multiBSChangeLeftKnotsd ( lastknotu-DEG, 3, 2, knotsv,
                                      pitch1, ccp, knotsa ) )
    goto failure;
  for ( i = 0; i < 2; i++ )
    knotsa[i] = (double)(lastknotv-DEG);
  knotsa[2] = lastknotv;
  if ( !mbs_multiBSChangeRightKnotsd ( lastknotu-DEG, 3, 2, lastknotv, knotsv,
                                       pitch1, ccp, knotsa ) )
    goto failure;
  knotsa[0] = 0.0;
  for ( i = 1; i < 3; i++ )
    knotsa[i] = 2.0;
  if ( !mbs_multiBSChangeLeftKnotsd ( 1, pitch1, 3, knotsu,
                                      0, ccp, knotsa ) )
    goto failure;
  for ( i = 0; i < 2; i++ )
    knotsa[i] = (double)(lastknotu-2);
  knotsa[2] = lastknotu;
  if ( !mbs_multiBSChangeRightKnotsd ( 1, pitch1, 2, lastknotu, knotsu,
                                       0, ccp, knotsa ) )
    goto failure;

        /* Now extract the control points, which determine */
        /* the boundary conditions */
  pkv_Moved ( lastknotu-4, 3, pitch1, -3*(lastknotv-4), &ccp[2*pitch1-3] );/* przesunięcie punktów kontrolnych brzegowych */
  pkv_Rearranged ( lastknotu-4, 6, pitch1, 6, &ccp[pitch1] ); /* zsunięcie punktów brzegowych*/
  memmove ( &ccp[pitch1 + 6*(lastknotu-4)], &ccp[(lastknotu-3)*pitch1], pitch1*sizeof(double) );
        /* number of points */
  k = 2*((lastknotv-2)+2*(lastknotu-4));
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
} /*g1bl_SurfNetDiameterSqd*/

/* //////////////////
jest obliczany w kroku met Newtona dolny prog na zmniejszenie ->
jesli nie zmniejsza dostatecz oblicza czy warto iterowac met Newtona
//////////////////////////////////////////////////////// */
boolean _g1bl_ComputeDeltaQd ( int n, const int *prof, double **hrows,
                               const double *grad, const double *dcoeff,
                               double *dq )
{
  void   *sp;
  double *aux;

  sp = pkv_GetScratchMemTop ();
  aux = pkv_GetScratchMemd ( 2*n );
  if ( !aux ) {
    PKV_SIGNALERROR ( LIB_G1BLENDING, ERRCODE_2, ERRMSG_2 );
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
} /*_g1bl_ComputeDeltaQd*/

/*

*/
boolean _g1bl_ShiftDecompHessiand ( int neqs, int hsize, int *prof, double *hessian,
                                    double *Lhessian, double **Lhrows, double nu )
{
  int i;

  memcpy ( Lhessian, hessian, hsize*sizeof(double) );
  for ( i = 0; i < neqs; i++ )
    Lhrows[i][i] += nu;
  return pkn_NRBSymCholeskyDecompd ( neqs, prof, Lhessian, Lhrows, NULL );
} /*_g1bl_ShiftDecompHessiand*/


/*
oblicza wartosc f(ni) z 4 punktu algorytmu -> krok LM
*/
double _g1bl_AuxNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int hsize, int *prof, double *hessian,
                          double *Lhessian, double **Lhrows, double nu,
                          int lastknotu, int lastknotv, int pitch,
                          point3d *cp, point3d *acp,
                          double *grad, double *dcoeff,
                          double tC, double *ftab )
{
  if ( _g1bl_ShiftDecompHessiand ( neqs, hsize, prof, hessian,
                                   Lhessian, Lhrows, nu ) ) {
    pkn_NRBLowerTrSolved ( neqs, prof, Lhessian, Lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( neqs, prof, Lhessian, Lhrows, 1, 1, dcoeff, 1, dcoeff );
    pkn_SubtractMatrixd ( lastknotu-DEG3, 3*(lastknotv-DEG3), pitch,
                          &cp[DEG*(pitch/3)+DEG].x,
                          3*(lastknotv-DEG3), dcoeff,
                          3*(lastknotv-DEG), &acp[DEG*(lastknotv-DEG)+DEG].x );
#ifdef _DEBUG
printf ( "F" );
#endif
    return g1bl_UFuncd ( nknots, qcoeff, Nitab, lastknotu, lastknotv,
                         3*(lastknotv-DEG), acp, NULL, tC, ftab );
  }
  else
    return MYINFINITY;
} /*_g1bl_AuxNuFuncd*/

