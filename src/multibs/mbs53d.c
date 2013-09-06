
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#undef CONST_
#define CONST_

void mbs_BSC2CoonsFindCornersd ( int spdimen,
          int degc00, int lastknotc00, const double *knotsc00, const double *c00,
          int degc01, int lastknotc01, const double *knotsc01, const double *c01,
          int degc02, int lastknotc02, const double *knotsc02, const double *c02,
          int degc10, int lastknotc10, const double *knotsc10, const double *c10,
          int degc11, int lastknotc11, const double *knotsc11, const double *c11,
          int degc12, int lastknotc12, const double *knotsc12, const double *c12,
          double *pcorners )
{
  double u0, u1;

  u0 = knotsc00[degc00];
  u1 = knotsc00[lastknotc00-degc00];
  mbs_multideBoorDer2d ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, u0,
      &pcorners[0], &pcorners[spdimen*12], &pcorners[spdimen*24] );
  mbs_multideBoorDer2d ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, u1,
      &pcorners[spdimen*6], &pcorners[spdimen*18], &pcorners[spdimen*30] );
  mbs_multideBoorDer2d ( degc10, lastknotc10, knotsc10, 1, spdimen, 0, c10, u0,
      &pcorners[spdimen*1], &pcorners[spdimen*13], &pcorners[spdimen*25] );
  mbs_multideBoorDer2d ( degc10, lastknotc10, knotsc10, 1, spdimen, 0, c10, u1,
      &pcorners[spdimen*7], &pcorners[spdimen*19], &pcorners[spdimen*31] );
  mbs_multideBoorDer2d ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, u0,
      &pcorners[spdimen*2], &pcorners[spdimen*14], &pcorners[spdimen*26] );
  mbs_multideBoorDer2d ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, u1,
      &pcorners[spdimen*8], &pcorners[spdimen*20], &pcorners[spdimen*32] );
  mbs_multideBoorDer2d ( degc11, lastknotc11, knotsc11, 1, spdimen, 0, c11, u0,
      &pcorners[spdimen*3], &pcorners[spdimen*15], &pcorners[spdimen*27] );
  mbs_multideBoorDer2d ( degc11, lastknotc11, knotsc11, 1, spdimen, 0, c11, u1,
      &pcorners[spdimen*9], &pcorners[spdimen*21], &pcorners[spdimen*33] );
  mbs_multideBoorDer2d ( degc02, lastknotc02, knotsc02, 1, spdimen, 0, c02, u0,
      &pcorners[spdimen*4], &pcorners[spdimen*16], &pcorners[spdimen*28] );
  mbs_multideBoorDer2d ( degc02, lastknotc02, knotsc02, 1, spdimen, 0, c02, u1,
      &pcorners[spdimen*10], &pcorners[spdimen*22], &pcorners[spdimen*34] );
  mbs_multideBoorDer2d ( degc12, lastknotc12, knotsc12, 1, spdimen, 0, c12, u0,
      &pcorners[spdimen*5], &pcorners[spdimen*17], &pcorners[spdimen*29] );
  mbs_multideBoorDer2d ( degc12, lastknotc12, knotsc12, 1, spdimen, 0, c12, u1,
      &pcorners[spdimen*11], &pcorners[spdimen*23], &pcorners[spdimen*35] );
} /*mbs_BSC2CoonsFindCornersd*/

/*
static void Verify ( int spdimen,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degc12, int lastknotc12, const double *knotsc12, const double *c12,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int degd12, int lastknotd12, const double *knotsd12, const double *d12 )
{
  void  *sp;
  double *cc, *cd, *dd;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMemd ( 36*spdimen );
  cd = pkv_GetScratchMemd ( 36*spdimen );
  dd = pkv_GetScratchMemd ( 36*spdimen );

  mbs_BSC2CoonsFindCornersd ( spdimen,
         degc00, lastknotc00, knotsc00, c00,
         degc01, lastknotc01, knotsc01, c01,
         degc02, lastknotc02, knotsc02, c02,
         degc10, lastknotc10, knotsc10, c10,
         degc11, lastknotc11, knotsc11, c11,
         degc12, lastknotc12, knotsc12, c12, cc );
  mbs_BSC2CoonsFindCornersd ( spdimen,
         degd00, lastknotd00, knotsd00, d00,
         degd01, lastknotd01, knotsd01, d01,
         degd02, lastknotd02, knotsd02, d02,
         degd10, lastknotd10, knotsd10, d10,
         degd11, lastknotd11, knotsd11, d11,
         degd12, lastknotd12, knotsd12, d12, cd );
  pkv_TransposeMatrixc ( 6, 6, spdimen*sizeof(double),
                         6*spdimen*sizeof(double), (char*)cd,
                         6*spdimen*sizeof(double), (char*)dd );
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 6; j++ )
      printf ( "%7.4f ", cc[(6*i+j)*spdimen] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 6; j++ )
      printf ( "%7.4f ", dd[(6*i+j)*spdimen] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 6; j++ )
      printf ( "%7.4f ", dd[(6*i+j)*spdimen]-cc[(6*i+j)*spdimen] );
    printf ( "\n" );
  }
  pkv_SetScratchMemTop ( sp );
} / *Verify*/

boolean mbs_BSC2CoonsToBSd ( int spdimen,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degc12, int lastknotc12, const double *knotsc12, const double *c12,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int degd12, int lastknotd12, const double *knotsd12, const double *d12,
      int *degreeu, int *lastuknot, double *uknots,
      int *degreev, int *lastvknot, double *vknots, double *p )
{
  void  *sp;
  int   degu, degv, lastukn, lastvkn, ptch, auxptch;
  double *auxuknots, *auxvknots, *auknots, *avknots;
  double *_cd00, *_cd01, *_cd02, *_cd10, *_cd11, *_cd12;
  double *pc, *bc;
  double *p1, *p2;
  int   i;

/*
Verify ( spdimen,
      degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
      degc02, lastknotc02, knotsc02, c02, degc10, lastknotc10, knotsc10, c10,
      degc11, lastknotc11, knotsc11, c11, degc12, lastknotc12, knotsc12, c12,
      degd00, lastknotd00, knotsd00, d00, degd01, lastknotd01, knotsd01, d01,
      degd02, lastknotd02, knotsd02, d02, degd10, lastknotd10, knotsd10, d10,
      degd11, lastknotd11, knotsd11, d11, degd12, lastknotd12, knotsd12, d12 );
*/
  sp = pkv_GetScratchMemTop ();
      /* find the degree and knot sequences of the final patch */
  degu = degv = 5;
        /* with respect to u */
  if ( !mbs_FindBSCommonKnotSequenced ( &degu, &lastukn, &auxuknots, 6,
               degc00, lastknotc00, knotsc00, degc01, lastknotc01, knotsc01,
               degc02, lastknotc02, knotsc02, degc10, lastknotc10, knotsc10,
               degc11, lastknotc11, knotsc11, degc12, lastknotc12, knotsc12 ) )
    goto failure;
  *degreeu = degu;
  *lastuknot = lastukn;
        /* with respect to v */
  if ( !mbs_FindBSCommonKnotSequenced ( &degv, &lastvkn, &auxvknots, 6,
               degd00, lastknotd00, knotsd00, degd01, lastknotd01, knotsd01,
               degd02, lastknotd02, knotsd02, degd10, lastknotd10, knotsd10,
               degd11, lastknotd11, knotsd11, degd12, lastknotd12, knotsd12 ) )
    goto failure;
  *degreev = degv;
  *lastvknot = lastvkn;

        /* allocate memory for the common representations of the curves */
  auknots = pkv_GetScratchMemd ( 24 );
  auxptch = spdimen*max ( lastukn-degu, lastvkn-degv );
  _cd00 = pkv_GetScratchMemd ( 6*auxptch );
  if ( !auknots || !_cd00 )
    goto failure;
  avknots = &auknots[12];

  pc = pkv_GetScratchMemd ( 36*spdimen );
  bc = pkv_GetScratchMemd ( 36*spdimen );
  p1 = pkv_GetScratchMemd ( (lastukn-degu)*(lastvkn-degv)*spdimen );
  p2 = pkv_GetScratchMemd ( (lastukn-degu)*(lastvkn-degv)*spdimen );
  if ( !pc || !bc || !p1 || !p2 )
    goto failure;

        /* construct the patch p1 */
  ptch = (lastukn-degu)*spdimen;
  _cd01 = &_cd00[ptch];  _cd02 = &_cd01[ptch];  _cd10 = &_cd02[ptch];
  _cd11 = &_cd10[ptch];  _cd12 = &_cd11[ptch];
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc00, lastknotc00, knotsc00,
     0, c00, degu, lastukn, auxuknots, 0, _cd00 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc01, lastknotc01, knotsc01,
     0, c01, degu, lastukn, auxuknots, 0, _cd01 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc02, lastknotc02, knotsc02,
     0, c02, degu, lastukn, auxuknots, 0, _cd02 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc10, lastknotc10, knotsc10,
     0, c10, degu, lastukn, auxuknots, 0, _cd10 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc11, lastknotc11, knotsc11,
     0, c11, degu, lastukn, auxuknots, 0, _cd11 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degc12, lastknotc12, knotsc12,
     0, c12, degu, lastukn, auxuknots, 0, _cd12 ) )
    goto failure;
  for ( i = 0; i < 6; i++ ) avknots[i] = auxvknots[degv];
  for ( i = 6; i < 12; i++ ) avknots[i] = auxvknots[lastvkn-degv];
  mbs_multiInterp2knHermiteBSd ( 1, ptch, 5, 11, avknots, 3, auxptch, _cd00,
           3, auxptch, _cd10, ptch, p1 );
  pkv_TransposeMatrixc ( 6, lastukn-degu, spdimen*sizeof(double),
                         (lastukn-degu)*spdimen*sizeof(double), (char*)p1,
                         6*spdimen*sizeof(double), (char*)p2 );
  if ( !mbs_multiAdjustBSCRepd ( lastukn-degu, spdimen,
            5, 11, avknots, 6*spdimen, p2, degv, lastvkn, auxvknots,
            (lastvkn-degv)*spdimen, p ) )
    goto failure;

        /* construct the patch p2 */
  ptch = (lastvkn-degv)*spdimen;
  _cd01 = &_cd00[ptch];  _cd02 = &_cd01[ptch];  _cd10 = &_cd02[ptch];
  _cd11 = &_cd10[ptch];  _cd12 = &_cd11[ptch];
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd00, lastknotd00, knotsd00,
     0, d00, degv, lastvkn, auxvknots, 0, _cd00 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd01, lastknotd01, knotsd01,
     0, d01, degv, lastvkn, auxvknots, 0, _cd01 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd02, lastknotd02, knotsd02,
     0, d02, degv, lastvkn, auxvknots, 0, _cd02 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd10, lastknotd10, knotsd10,
     0, d10, degv, lastvkn, auxvknots, 0, _cd10 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd11, lastknotd11, knotsd11,
     0, d11, degv, lastvkn, auxvknots, 0, _cd11 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, spdimen, degd12, lastknotd12, knotsd12,
     0, d12, degv, lastvkn, auxvknots, 0, _cd12 ) )
    goto failure;
  for ( i = 0; i < 6; i++ ) auknots[i] = auxuknots[degu];
  for ( i = 6; i < 12; i++ ) auknots[i] = auxuknots[lastukn-degu];
  mbs_multiInterp2knHermiteBSd ( 1, ptch, 5, 11, auknots, 3, auxptch, _cd00,
           3, auxptch, _cd10, ptch, p2 );
  if ( !mbs_multiAdjustBSCRepd ( 1, ptch, 5, 11, auknots, 0, p2,
            degu, lastukn, auxuknots, 0, p1 ) )
    goto failure;
  pkn_AddMatrixd ( 1, ptch*(lastukn-degu), 0, p, 0, p1, 0, p );

        /* construct the patch p3 */
  mbs_BSC2CoonsFindCornersd ( spdimen,
        degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
        degc02, lastknotc02, knotsc02, c02, degc10, lastknotc10, knotsc10, c10,
        degc11, lastknotc11, knotsc11, c11, degc12, lastknotc12, knotsc12, c12,
        pc );

  pkv_Selectd ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, pc, p1 );
  pkv_Selectd ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, &pc[6*spdimen], &p1[6*3*spdimen] );
  for ( i = 0; i < 6; i++ ) {
    pkv_Selectd ( 3, spdimen, 2*spdimen, spdimen, &p1[6*i*spdimen], &bc[6*i*spdimen] );
    pkv_Selectd ( 3, spdimen, 2*spdimen, spdimen, &p1[(6*i+1)*spdimen], &bc[(6*i+3)*spdimen] );
  }
  mbs_multiInterp2knHermiteBSd ( 1, 6*spdimen, 5, 11, auknots, 3, 0, bc,
                                 3, 0, &bc[6*3*spdimen], 0, p1 );
  mbs_multiInterp2knHermiteBSd ( 6, spdimen, 5, 11, avknots, 3, 6*spdimen, p1,
                                 3, 6*spdimen, &p1[3*spdimen], 6*spdimen, p2 );
  if ( !mbs_multiAdjustBSCRepd ( 6, spdimen, 5, 11, avknots, 6*spdimen, p2,
               degv, lastvkn, auxvknots, (lastvkn-degv)*spdimen, p1 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepd ( 1, ptch, 5, 11, auknots, 0, p1,
                                 degu, lastukn, auxuknots, 0, p2 ) )
    goto failure;
  pkn_SubtractMatrixd ( 1, ptch*(lastukn-degu), 0, p, 0, p2, 0, p );

  memcpy ( uknots, auxuknots, (lastukn+1)*sizeof(double) );
  memcpy ( vknots, auxvknots, (lastvkn+1)*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BSC2CoonsToBSd*/

/* ///////////////////////////////////////////////////////////////////////// */
void mbs_TabBSCurveDer3d ( int spdimen, int degree, int lastknot,
                           const double *knots, const double *cp,
                           int nkn, const double *kn, int ppitch,
                           double *p, double *dp, double *ddp, double *dddp )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    mbs_multideBoorDer3d ( degree, lastknot, knots, 1, spdimen, 0, cp,
                           kn[i], &p[j], &dp[j], &ddp[j], &dddp[j] );
} /*mbs_TabBSCurveDer3d*/

boolean _mbs_TabBSC2Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp )
{
  void  *sp;
  int   i, j, k, l;
  double *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemd ( 6*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 6*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 6; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(6*i+j)*spdimen+l] -= hv[6*i+k]*p[(6*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(6*i+k)*spdimen+l]*hv[6*j+k] +
                                      a[(6*j+k)*spdimen+l]*hu[6*i+k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBSC2Coonsd*/

boolean mbs_TabBSC2CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degc12, int lastknotc12, const double *knotsc12, const double *c12,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int degd12, int lastknotd12, const double *knotsd12, const double *d12,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void *sp;
  int  ku, kv;
  double *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = pkv_GetScratchMemd ( (24*(nknu+nknv)+36)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  mbs_BSC2CoonsFindCornersd ( spdimen,
          degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
          degc02, lastknotc02, knotsc02, c02, degc10, lastknotc10, knotsc10, c10,
          degc11, lastknotc11, knotsc11, c11, degc12, lastknotc12, knotsc12, c12,
          pcorners );

  mbs_TabBSCurveDer3d ( spdimen, degc00, lastknotc00, knotsc00, c00, nknu, knu,
      6*spdimen, &c[0], &dc[0], &ddc[0], &dddc[0] );
  mbs_TabBSCurveDer3d ( spdimen, degc10, lastknotc10, knotsc10, c10, nknu, knu,
      6*spdimen, &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degc01, lastknotc01, knotsc01, c01, nknu, knu,
      6*spdimen, &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degc11, lastknotc11, knotsc11, c11, nknu, knu,
      6*spdimen, &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degc02, lastknotc02, knotsc02, c02, nknu, knu,
      6*spdimen, &c[4*spdimen], &dc[4*spdimen], &ddc[4*spdimen], &dddc[4*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degc12, lastknotc12, knotsc12, c12, nknu, knu,
      6*spdimen, &c[5*spdimen], &dc[5*spdimen], &ddc[5*spdimen], &dddc[5*spdimen] );

  mbs_TabBSCurveDer3d ( spdimen, degd00, lastknotd00, knotsd00, d00, nknv, knv,
      6*spdimen, &d[0], &dd[0], &ddd[0], &dddd[0] );
  mbs_TabBSCurveDer3d ( spdimen, degd10, lastknotd10, knotsd10, d10, nknv, knv,
      6*spdimen, &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd01, lastknotd01, knotsd01, d01, nknv, knv,
      6*spdimen, &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd11, lastknotd11, knotsd11, d11, nknv, knv,
      6*spdimen, &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd02, lastknotd02, knotsd02, d02, nknv, knv,
      6*spdimen, &d[4*spdimen], &dd[4*spdimen], &ddd[4*spdimen], &dddd[4*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd12, lastknotd12, knotsd12, d12, nknv, knv,
      6*spdimen, &d[5*spdimen], &dd[5*spdimen], &ddd[5*spdimen], &dddd[5*spdimen] );

  if ( p )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, c, d, pcorners,
                               hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, dc, d, pcorners,
                               dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, c, dd, pcorners,
                               hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, ddc, d, pcorners,
                               ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, dc, dd, pcorners,
                               dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, c, ddd, pcorners,
                               hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, dddc, d, pcorners,
                               dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, ddc, dd, pcorners,
                               ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, dc, ddd, pcorners,
                               dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBSC2Coonsd ( spdimen, nknu, nknv, c, dddd, pcorners,
                               hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBSC2CoonsDer3d*/

boolean _mbs_TabBSC2Coons0d ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp )
{
  void  *sp;
  int   i, j, k, l;
  double *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemd ( 3*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 3*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(3*i+j)*spdimen+l] -= hv[6*i+2*k]*p[(3*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 3; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(3*i+k)*spdimen+l]*hv[6*j+2*k] +
                                      a[(3*j+k)*spdimen+l]*hu[6*i+2*k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBSC2Coons0d*/

boolean mbs_TabBSC2Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void   *sp;
  int    ku, kv;
  double *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 3*spdimen*nknu;
  kv = 3*spdimen*nknv;
  c = pkv_GetScratchMemd ( (24*(nknu+nknv)+9)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  mbs_multideBoorDer2d ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*3], &pcorners[spdimen*6] );
  mbs_multideBoorDer2d ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*4], &pcorners[spdimen*7] );
  mbs_multideBoorDer2d ( degc02, lastknotc02, knotsc02, 1, spdimen, 0, c02, 0.0,
      &pcorners[spdimen*2], &pcorners[spdimen*5], &pcorners[spdimen*8] );

  mbs_TabBSCurveDer3d ( spdimen, degc00, lastknotc00, knotsc00, c00, nknu, knu,
      3*spdimen, &c[0], &dc[0], &ddc[0], &dddc[0] );
  mbs_TabBSCurveDer3d ( spdimen, degc01, lastknotc01, knotsc01, c01, nknu, knu,
      3*spdimen, &c[1*spdimen], &dc[1*spdimen], &ddc[1*spdimen], &dddc[1*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degc02, lastknotc02, knotsc02, c02, nknu, knu,
      3*spdimen, &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd00, lastknotd00, knotsd00, d00, nknv, knv,
      3*spdimen, &d[0], &dd[0], &ddd[0], &dddd[0] );
  mbs_TabBSCurveDer3d ( spdimen, degd01, lastknotd01, knotsd01, d01, nknv, knv,
      3*spdimen, &d[1*spdimen], &dd[1*spdimen], &ddd[1*spdimen], &dddd[1*spdimen] );
  mbs_TabBSCurveDer3d ( spdimen, degd02, lastknotd02, knotsd02, d02, nknv, knv,
      3*spdimen, &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen] );

  if ( p )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, dddc, d, pcorners,
                                dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBSC2Coons0d ( spdimen, nknu, nknv, c, dddd, pcorners,
                                hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBSC2Coons0Der3d*/

