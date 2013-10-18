
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
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_BezC1CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     double *pcorners )
{
  if ( !mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 0.0,
                          &pcorners[0], &pcorners[spdimen*8] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 1.0,
                          &pcorners[spdimen*4], &pcorners[spdimen*12] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc10, 1, spdimen, 0, c10, 0.0,
                          &pcorners[spdimen*1], &pcorners[spdimen*9] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc10, 1, spdimen, 0, c10, 1.0,
                          &pcorners[spdimen*5], &pcorners[spdimen*13] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 0.0,
                          &pcorners[spdimen*2], &pcorners[spdimen*10] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 1.0,
                          &pcorners[spdimen*6], &pcorners[spdimen*14] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc11, 1, spdimen, 0, c11, 0.0,
                          &pcorners[spdimen*3], &pcorners[spdimen*11] ) )
    return false;
  if ( !mbs_multiBCHornerDerd ( degc11, 1, spdimen, 0, c11, 1.0,
                          &pcorners[spdimen*7], &pcorners[spdimen*15] ) )
    return false;
  return true;
} /*mbs_BezC1CoonsFindCornersd*/

static int FindMaxInt ( int a, int b, int c, int d, int e )
{
  a = max ( a, b );  a = max ( a, c );  a = max ( a, d );  a = max ( a, e );
  return a;
} /*FindMaxInt*/

boolean mbs_BezC1CoonsToBezd ( int spdimen,
                               int degc00, const double *c00,
                               int degc01, const double *c01,
                               int degc10, const double *c10,
                               int degc11, const double *c11,
                               int degd00, const double *d00,
                               int degd01, const double *d01,
                               int degd10, const double *d10,
                               int degd11, const double *d11,
                               int *n, int *m, double *p )
{
  void   *sp;
  int    degu, degv, d, du, dv;
  int    pitch;
  double *pc, *bc;
  double *p1, *p2, *p3; 
  int    i;

  sp = pkv_GetScratchMemTop ();
  *n = degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  *m = degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  pc = pkv_GetScratchMemd ( 16*spdimen );
  p1 = pkv_GetScratchMemd ( 3*(degu+1)*(degv+1)*spdimen );
  bc = pkv_GetScratchMemd ( 4*(max(degu,degv)+1)*spdimen );
  if ( !pc || !p1 || !bc )
    goto failure;
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];

        /* construct the patch p1 */
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc00, c00, degu-degc00,
                              0, &d, &bc[0] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc01, c01, degu-degc01,
                              0, &d, &bc[(degu+1)*spdimen] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc10, c10, degu-degc10,
                              0, &d, &bc[2*(degu+1)*spdimen] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc11, c11, degu-degc11,
                              0, &d, &bc[3*(degu+1)*spdimen] ) )
    goto failure;
  pitch = (degu+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezd ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degu+1)*spdimen], pitch, p1 ) )
    goto failure;
  pkv_TransposeMatrixc ( 4, degu+1, spdimen*sizeof(double),
                         (degu+1)*spdimen*sizeof(double), (char*)p1,
                         4*spdimen*sizeof(double), (char*)p );
  if ( !mbs_BCDegElevPd ( spdimen, degu, 3, p, 0, degv-3, &du, &dv, p ) )
    goto failure;

        /* construct the patch p2 */
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd00, d00, degv-degd00,
                              0, &d, &bc[0] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd01, d01, degv-degd01,
                              0, &d, &bc[(degv+1)*spdimen] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd10, d10, degv-degd10,
                              0, &d, &bc[2*(degv+1)*spdimen] ) )
    goto failure;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd11, d11, degv-degd11,
                              0, &d, &bc[3*(degv+1)*spdimen] ) )
    goto failure;
  pitch = (degv+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezd ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degv+1)*spdimen], pitch, p2 ) )
    goto failure;
  if ( !mbs_BCDegElevPd ( spdimen, 3, degv, p2, degu-3, 0, &du, &dv, p2 ) )
    goto failure;
  pkn_AddMatrixd ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  if ( !mbs_BezC1CoonsFindCornersd ( spdimen, degc00, c00, degc01, c01,
                                     degc10, c10, degc11, c11, pc ) )
    goto failure;

  pkv_Selectd ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, pc, p3 );
  pkv_Selectd ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, &pc[4*spdimen], &p3[4*2*spdimen] );
  for ( i = 0; i < 4; i++ ) {
    pkv_Selectd ( 2, spdimen, 2*spdimen, spdimen, &p3[4*i*spdimen], &bc[4*i*spdimen] );
    pkv_Selectd ( 2, spdimen, 2*spdimen, spdimen, &p3[(4*i+1)*spdimen], &bc[(4*i+2)*spdimen] );
  }
  if ( !mbs_multiInterp2knHermiteBezd ( 1, 4*spdimen, 3, 2, 0, bc,
                                        2, 0, &bc[4*2*spdimen], 0, p3 ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezd ( 4, spdimen, 3, 2, 4*spdimen, p3,
                                  2, 4*spdimen, &p3[2*spdimen], 4*spdimen, p2 ) )
    goto failure;
  if ( !mbs_BCDegElevPd ( spdimen, 3, 3, p2, degu-3, degv-3, &du, &dv, p3 ) )
    goto failure;
  pkn_SubtractMatrixd ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BezC1CoonsToBezd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean mbs_TabCubicHFuncDer2d ( double a, double b, int nkn, const double *kn,
                                 double *hfunc, double *dhfunc, double *ddhfunc )
{
  double HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                      0.0, 0.0, 1.0, 1.0,                  /* h10 */
                      0.0, (double)(1.0/3.0), 0.0, 0.0,    /* h01 */
                      0.0, 0.0, (double)(-1.0/3.0), 0.0};  /* h11 */
  int    i, j;
  double h, h_1, h_2;

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       if ( !mbs_multiBCHornerDer2d ( 3, 4, 1, 4, HFunc, kn[i],
                                      &hfunc[j], &dhfunc[j], &ddhfunc[j] ) )
         return false;
  }
  else {
    h = b - a;
    h_1 = 1.0/h;
    h_2 = h_1*h_1;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      if ( !mbs_multiBCHornerDer2d ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                                     &hfunc[j], &dhfunc[j], &ddhfunc[j] ) )
        return false;
      hfunc[j+2]   *= h;      hfunc[j+3]   *= h;
      dhfunc[j]    *= h_1;    dhfunc[j+1]  *= h_1;  
      ddhfunc[j]   *= h_2;    ddhfunc[j+1] *= h_2;  
      ddhfunc[j+2] *= h_1;    ddhfunc[j+3] *= h_1;  
    }
  }
  return true;
} /*mbs_TabCubicHFuncDer2d*/

boolean mbs_TabCubicHFuncDer3d ( double a, double b, int nkn, const double *kn,
                                 double *hfunc, double *dhfunc, double *ddhfunc,
                                 double *dddhfunc )
{
  double HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                     0.0, 0.0, 1.0, 1.0,                  /* h10 */
                     0.0, (double)(1.0/3.0), 0.0, 0.0,     /* h01 */
                     0.0, 0.0, (double)(-1.0/3.0), 0.0};   /* h11 */
  int   i, j;
  double h, h_1, h_2, h_3;

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       if ( !mbs_multiBCHornerDer3d ( 3, 4, 1, 4, HFunc, kn[i],
                            &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
         return false;
  }
  else {
    h = b - a;
    h_1 = (double)1.0/h;
    h_2 = h_1*h_1;
    h_3 = h_1*h_2;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      if ( !mbs_multiBCHornerDer3d ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                           &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
        return false;
      hfunc[j+2]    *= h;      hfunc[j+3]   *= h;
      dhfunc[j]     *= h_1;    dhfunc[j+1]  *= h_1;
      ddhfunc[j]    *= h_2;    ddhfunc[j+1] *= h_2;
      ddhfunc[j+2]  *= h_1;    ddhfunc[j+3] *= h_1;
      dddhfunc[j]   *= h_3;    dddhfunc[j+1] *= h_3;
      dddhfunc[j+2] *= h_2;    dddhfunc[j+3] *= h_2;
    }
  }
  return true;
} /*mbs_TabCubicHFuncDer3d*/

boolean mbs_TabBezCurveDer2d ( int spdimen, int degree, const double *cp,
                               int nkn, const double *kn,
                               int ppitch,
                               double *p, double *dp, double *ddp )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    if ( !mbs_multiBCHornerDer2d ( degree, 1, spdimen, 0, cp, kn[i],
                                   &p[j], &dp[j], &ddp[j] ) )
      return false;
  return true;
} /*mbs_TabBezCurveDer2d*/

boolean _mbs_TabBezC1Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp )
{
  void   *sp;
  int    i, j, k, l;
  double *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemd ( 4*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 4*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(4*i+j)*spdimen+l] -= hv[4*i+k]*p[(4*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(double) );
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
} /*_mbs_TabBezC1Coonsd*/

boolean mbs_TabBezC1CoonsDer2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv )
{
  void   *sp;
  int    ku, kv;
  double *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 4*spdimen*nknu;
  kv = 4*spdimen*nknv;
  c = pkv_GetScratchMemd ( (18*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  if ( !mbs_BezC1CoonsFindCornersd ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners ) )
    goto failure;

  if ( !mbs_TabBezCurveDer2d ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen] ) )
    goto failure;

  if ( !mbs_TabBezCurveDer2d ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen] ) )
    goto failure;

  if ( p )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1CoonsDer2d*/

boolean mbs_TabBezC1CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void   *sp;
  int    ku, kv;
  double *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = pkv_GetScratchMemd ( (24*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  if ( !mbs_BezC1CoonsFindCornersd ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners ) )
    goto failure;

  if ( !mbs_TabBezCurveDer3d ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen] ) )
    goto failure;

  if ( !mbs_TabBezCurveDer3d ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen] ) )
    goto failure;

  if ( p )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dddc, d, pcorners,
                                dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBezC1Coonsd ( spdimen, nknu, nknv, c, dddd, pcorners,
                                hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1CoonsDer3d*/

boolean _mbs_TabBezC1Coons0d (
                   int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp )
{
  void   *sp;
  int    i, j, k, l;
  double *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemd ( 2*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 2*nknv*spdimen*sizeof(double) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(2*i+j)*spdimen+l] -= hv[4*i+2*k]*p[(2*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(double) );
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
} /*_mbs_TabBezC1Coons0d*/

boolean mbs_TabBezC1Coons0Der2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv )
{
  void   *sp;
  int    ku, kv;
  double *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 2*spdimen*nknu;
  kv = 2*spdimen*nknv;
  c = pkv_GetScratchMemd ( (9*(nknu+nknv)+4)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  if ( !mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2] ) )
    goto failure;
  if ( !mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3] ) )
    goto failure;

  if ( !mbs_TabBezCurveDer2d ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer2d ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen] ) )
    goto failure;

  if ( p )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1Coons0Der2d*/

boolean mbs_TabBezC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void   *sp;
  int    ku, kv;
  double *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = pkv_GetScratchMemd ( (24*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  if ( !mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2] ) )
    goto failure;
  if ( !mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3] ) )
    goto failure;

  if ( !mbs_TabBezCurveDer3d ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0] ) )
    goto failure;
  if ( !mbs_TabBezCurveDer3d ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen] ) )
    goto failure;

  if ( p )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dddc, d, pcorners,
                                 dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                 ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                 dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, dddd, pcorners,
                                 hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1Coons0Der3d*/

