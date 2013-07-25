
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "msgpool.h"

#undef CONST_
#define CONST_

void mbs_BSC1CoonsFindCornersf ( int spdimen,
          int degc00, int lastknotc00, const float *knotsc00, const float *c00,
          int degc01, int lastknotc01, const float *knotsc01, const float *c01,
          int degc10, int lastknotc10, const float *knotsc10, const float *c10,
          int degc11, int lastknotc11, const float *knotsc11, const float *c11,
          float *pcorners )
{
  float u0, u1;

  u0 = knotsc00[degc00];
  u1 = knotsc00[lastknotc00-degc00];
  mbs_multideBoorDerf ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, u0,
      &pcorners[0], &pcorners[spdimen*8] );
  mbs_multideBoorDerf ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, u1,
      &pcorners[spdimen*4], &pcorners[spdimen*12] );
  mbs_multideBoorDerf ( degc10, lastknotc10, knotsc10, 1, spdimen, 0, c10, u0,
      &pcorners[spdimen*1], &pcorners[spdimen*9] );
  mbs_multideBoorDerf ( degc10, lastknotc10, knotsc10, 1, spdimen, 0, c10, u1,
      &pcorners[spdimen*5], &pcorners[spdimen*13] );
  mbs_multideBoorDerf ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, u0,
      &pcorners[spdimen*2], &pcorners[spdimen*10] );
  mbs_multideBoorDerf ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, u1,
      &pcorners[spdimen*6], &pcorners[spdimen*14] );
  mbs_multideBoorDerf ( degc11, lastknotc11, knotsc11, 1, spdimen, 0, c11, u0,
      &pcorners[spdimen*3], &pcorners[spdimen*11] );
  mbs_multideBoorDerf ( degc11, lastknotc11, knotsc11, 1, spdimen, 0, c11, u1,
      &pcorners[spdimen*7], &pcorners[spdimen*15] );
} /*mbs_BSC1CoonsFindCornersf*/

/*
static void Verify ( int spdimen,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11 )
{
  void  *sp;
  float *cc, *cd, *dd;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMemf ( 16*spdimen );
  cd = pkv_GetScratchMemf ( 16*spdimen );
  dd = pkv_GetScratchMemf ( 16*spdimen );

  mbs_BSC1CoonsFindCornersf ( spdimen,
         degc00, lastknotc00, knotsc00, c00,
         degc01, lastknotc01, knotsc01, c01,
         degc10, lastknotc10, knotsc10, c10,
         degc11, lastknotc11, knotsc11, c11, cc );
  mbs_BSC1CoonsFindCornersf ( spdimen,
         degd00, lastknotd00, knotsd00, d00,
         degd01, lastknotd01, knotsd01, d01,
         degd10, lastknotd10, knotsd10, d10,
         degd11, lastknotd11, knotsd11, d11, cd );
  pkv_TransposeMatrixc ( 4, 4, spdimen*sizeof(float),
                         4*spdimen*sizeof(float), (char*)cd,
                         4*spdimen*sizeof(float), (char*)dd );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%7.4f ", cc[(4*i+j)*spdimen] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%7.4f ", dd[(4*i+j)*spdimen] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%7.4f ", dd[(4*i+j)*spdimen]-cc[(4*i+j)*spdimen] );
    printf ( "\n" );
  }
  pkv_SetScratchMemTop ( sp );
} / *Verify*/

boolean mbs_BSC1CoonsToBSf ( int spdimen,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      int *degreeu, int *lastuknot, float *uknots,
      int *degreev, int *lastvknot, float *vknots, float *p )
{
  void  *sp;
  int   degu, degv, lastukn, lastvkn, ptch, auxptch;
  float *auxuknots, *auxvknots, *auknots, *avknots;
  float *_cd00, *_cd01, *_cd10, *_cd11;
  float *pc, *bc;
  float *p1, *p2;
  int   i;

/*
Verify ( spdimen,
      degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
      degc10, lastknotc10, knotsc10, c10, degc11, lastknotc11, knotsc11, c11,
      degd00, lastknotd00, knotsd00, d00, degd01, lastknotd01, knotsd01, d01,
      degd10, lastknotd10, knotsd10, d10, degd11, lastknotd11, knotsd11, d11 );
*/
  sp = pkv_GetScratchMemTop ();
      /* find the degree and knot sequences of the final patch */
  degu = degv = 3;
        /* with respect to u */
  if ( !mbs_FindBSCommonKnotSequencef ( &degu, &lastukn, &auxuknots, 4,
               degc00, lastknotc00, knotsc00, degc01, lastknotc01, knotsc01,
               degc10, lastknotc10, knotsc10, degc11, lastknotc11, knotsc11 ) )
    goto failure;
  *degreeu = degu;
  *lastuknot = lastukn;
        /* with respect to v */
  if ( !mbs_FindBSCommonKnotSequencef ( &degv, &lastvkn, &auxvknots, 4,
               degd00, lastknotd00, knotsd00, degd01, lastknotd01, knotsd01,
               degd10, lastknotd10, knotsd10, degd11, lastknotd11, knotsd11 ) )
    goto failure;
  *degreev = degv;
  *lastvknot = lastvkn;

        /* allocate memory for the common representations of the curves */
  auknots = pkv_GetScratchMemf ( 16 );
  auxptch = spdimen*max ( lastukn-degu, lastvkn-degv );
  _cd00 = pkv_GetScratchMemf ( 4*auxptch );
  if ( !auknots || !_cd00 )
    goto failure;
  avknots = &auknots[8];

  pc = pkv_GetScratchMemf ( 16*spdimen );
  bc = pkv_GetScratchMemf ( 16*spdimen );
  p1 = pkv_GetScratchMemf ( (lastukn-degu)*(lastvkn-degv)*spdimen );
  p2 = pkv_GetScratchMemf ( (lastukn-degu)*(lastvkn-degv)*spdimen );
  if ( !pc || !bc || !p1 || !p2 )
    goto failure;

        /* construct the patch p1 */
  ptch = (lastukn-degu)*spdimen;
  _cd01 = &_cd00[ptch];  _cd10 = &_cd01[ptch];  _cd11 = &_cd10[ptch];
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degc00, lastknotc00, knotsc00,
     0, c00, degu, lastukn, auxuknots, 0, _cd00 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degc01, lastknotc01, knotsc01,
     0, c01, degu, lastukn, auxuknots, 0, _cd01 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degc10, lastknotc10, knotsc10,
     0, c10, degu, lastukn, auxuknots, 0, _cd10 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degc11, lastknotc11, knotsc11,
     0, c11, degu, lastukn, auxuknots, 0, _cd11 ) )
    goto failure;
  for ( i = 0; i < 4; i++ ) avknots[i] = auxvknots[degv];
  for ( i = 4; i < 8; i++ ) avknots[i] = auxvknots[lastvkn-degv];
  mbs_multiInterp2knHermiteBSf ( 1, ptch, 3, 7, avknots, 2, auxptch, _cd00,
           2, auxptch, _cd10, ptch, p1 );
  pkv_TransposeMatrixc ( 4, lastukn-degu, spdimen*sizeof(float),
                         (lastukn-degu)*spdimen*sizeof(float), (char*)p1,
                         4*spdimen*sizeof(float), (char*)p2 );
  if ( !mbs_multiAdjustBSCRepf ( lastukn-degu, spdimen,
            3, 7, avknots, 4*spdimen, p2, degv, lastvkn, auxvknots,
            (lastvkn-degv)*spdimen, p ) )
    goto failure;

        /* construct the patch p2 */
  ptch = (lastvkn-degv)*spdimen;
  _cd01 = &_cd00[ptch];  _cd10 = &_cd01[ptch];  _cd11 = &_cd10[ptch];
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degd00, lastknotd00, knotsd00,
     0, d00, degv, lastvkn, auxvknots, 0, _cd00 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degd01, lastknotd01, knotsd01,
     0, d01, degv, lastvkn, auxvknots, 0, _cd01 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degd10, lastknotd10, knotsd10,
     0, d10, degv, lastvkn, auxvknots, 0, _cd10 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, spdimen, degd11, lastknotd11, knotsd11,
     0, d11, degv, lastvkn, auxvknots, 0, _cd11 ) )
    goto failure;
  for ( i = 0; i < 4; i++ ) auknots[i] = auxuknots[degu];
  for ( i = 4; i < 8; i++ ) auknots[i] = auxuknots[lastukn-degu];
  mbs_multiInterp2knHermiteBSf ( 1, ptch, 3, 7, auknots, 2, auxptch, _cd00,
           2, auxptch, _cd10, ptch, p2 );
  if ( !mbs_multiAdjustBSCRepf ( 1, ptch, 3, 7, auknots, 0, p2,
            degu, lastukn, auxuknots, 0, p1 ) )
    goto failure;
  pkn_AddMatrixf ( 1, ptch*(lastukn-degu), 0, p, 0, p1, 0, p );

        /* construct the patch p3 */
  mbs_BSC1CoonsFindCornersf ( spdimen,
        degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
        degc10, lastknotc10, knotsc10, c10, degc11, lastknotc11, knotsc11, c11,
        pc );

  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, pc, p1 );
  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, &pc[4*spdimen], &p1[4*2*spdimen] );
  for ( i = 0; i < 4; i++ ) {
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p1[4*i*spdimen], &bc[4*i*spdimen] );
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p1[(4*i+1)*spdimen], &bc[(4*i+2)*spdimen] );
  }
  mbs_multiInterp2knHermiteBSf ( 1, 4*spdimen, 3, 7, auknots, 2, 0, bc,
                                 2, 0, &bc[4*2*spdimen], 0, p1 );
  mbs_multiInterp2knHermiteBSf ( 4, spdimen, 3, 7, avknots, 2, 4*spdimen, p1,
                                 2, 4*spdimen, &p1[2*spdimen], 4*spdimen, p2 );
  if ( !mbs_multiAdjustBSCRepf ( 4, spdimen, 3, 7, avknots, 4*spdimen, p2,
               degv, lastvkn, auxvknots, (lastvkn-degv)*spdimen, p1 ) )
    goto failure;
  if ( !mbs_multiAdjustBSCRepf ( 1, ptch, 3, 7, auknots, 0, p1,
                                 degu, lastukn, auxuknots, 0, p2 ) )
    goto failure;
  pkn_SubtractMatrixf ( 1, ptch*(lastukn-degu), 0, p, 0, p2, 0, p );

  memcpy ( uknots, auxuknots, (lastukn+1)*sizeof(float) );
  memcpy ( vknots, auxvknots, (lastvkn+1)*sizeof(float) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BSC1CoonsToBSf*/

/* ///////////////////////////////////////////////////////////////////////// */
void mbs_TabBSCurveDer2f ( int spdimen, int degree, int lastknot,
                           const float *knots, const float *cp,
                           int nkn, const float *kn, int ppitch,
                           float *p, float *dp, float *ddp )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    mbs_multideBoorDer2f ( degree, lastknot, knots, 1, spdimen, 0, cp,
                           kn[i], &p[j], &dp[j], &ddp[j] );
} /*mbs_TabBSCurveDer2f*/

boolean _mbs_TabBSC1Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp )
{
  void  *sp;
  int   i, j, k, l;
  float *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 4*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 4*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(4*i+j)*spdimen+l] -= hv[4*i+k]*p[(4*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(4*i+k)*spdimen+l]*hv[4*j+k] +
                                      a[(4*j+k)*spdimen+l]*hu[4*i+k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBSC1Coonsf*/

boolean mbs_TabBSC1CoonsDer2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv )
{
  void *sp;
  int  ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 4*spdimen*nknu;
  kv = 4*spdimen*nknv;
  c = pkv_GetScratchMemf ( (18*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  mbs_BSC1CoonsFindCornersf ( spdimen,
          degc00, lastknotc00, knotsc00, c00, degc01, lastknotc01, knotsc01, c01,
          degc10, lastknotc10, knotsc10, c10, degc11, lastknotc11, knotsc11, c11,
          pcorners );

  mbs_TabBSCurveDer2f ( spdimen, degc00, lastknotc00, knotsc00, c00, nknu, knu,
      4*spdimen, &c[0], &dc[0], &ddc[0] );
  mbs_TabBSCurveDer2f ( spdimen, degc10, lastknotc10, knotsc10, c10, nknu, knu,
      4*spdimen, &c[spdimen], &dc[spdimen], &ddc[spdimen] );
  mbs_TabBSCurveDer2f ( spdimen, degc01, lastknotc01, knotsc01, c01, nknu, knu,
      4*spdimen, &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen] );
  mbs_TabBSCurveDer2f ( spdimen, degc11, lastknotc11, knotsc11, c11, nknu, knu,
      4*spdimen, &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen] );

  mbs_TabBSCurveDer2f ( spdimen, degd00, lastknotd00, knotsd00, d00, nknv, knv,
      4*spdimen, &d[0], &dd[0], &ddd[0] );
  mbs_TabBSCurveDer2f ( spdimen, degd10, lastknotd10, knotsd10, d10, nknv, knv,
      4*spdimen, &d[spdimen], &dd[spdimen], &ddd[spdimen] );
  mbs_TabBSCurveDer2f ( spdimen, degd01, lastknotd01, knotsd01, d01, nknv, knv,
      4*spdimen, &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen] );
  mbs_TabBSCurveDer2f ( spdimen, degd11, lastknotd11, knotsd11, d11, nknv, knv,
      4*spdimen, &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen] );

  if ( p )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                               hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                               dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                               hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                               ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                               dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBSC1Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                               hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBSC1CoonsDer2f*/

boolean _mbs_TabBSC1Coons0f ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp )
{
  void  *sp;
  int   i, j, k, l;
  float *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 2*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 2*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(2*i+j)*spdimen+l] -= hv[4*i+2*k]*p[(2*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(2*i+k)*spdimen+l]*hv[4*j+2*k] +
                                      a[(2*j+k)*spdimen+l]*hu[4*i+2*k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBSC1Coons0f*/

boolean mbs_TabBSC1Coons0Der2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv )
{
  void *sp;
  int  ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 3*spdimen*nknu;
  kv = 3*spdimen*nknv;
  c = pkv_GetScratchMemf ( (9*(nknu+nknv)+4)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  mbs_multideBoorDerf ( degc00, lastknotc00, knotsc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2] );
  mbs_multideBoorDerf ( degc01, lastknotc01, knotsc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3] );

  mbs_TabBSCurveDer2f ( spdimen, degc00, lastknotc00, knotsc00, c00, nknu, knu,
      2*spdimen, &c[0], &dc[0], &ddc[0] );
  mbs_TabBSCurveDer2f ( spdimen, degc01, lastknotc01, knotsc01, c01, nknu, knu,
      2*spdimen, &c[1*spdimen], &dc[1*spdimen], &ddc[1*spdimen] );
  mbs_TabBSCurveDer2f ( spdimen, degd00, lastknotd00, knotsd00, d00, nknv, knv,
      2*spdimen, &d[0], &dd[0], &ddd[0] );
  mbs_TabBSCurveDer2f ( spdimen, degd01, lastknotd01, knotsd01, d01, nknv, knv,
      2*spdimen, &d[1*spdimen], &dd[1*spdimen], &ddd[1*spdimen] );

  if ( p )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBSC1Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBSC1Coons0Der2f*/

