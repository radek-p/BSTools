
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2013                            */
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

#undef CONST_
#define CONST_

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

#define _DEBUG

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_TabBSFuncDer3d ( int deg, int lastknot, const double *knots,
                              int i0, int i1,
                              int n, const double *tkn, int *fkn, int *lkn,
                              double *b, double *bt, double *btt, double *bttt )
{
  /* evaluate B-spline basis functions at a number of points  deg - degree,  */
  /* lastknot - number of the last knot, knots - knots, i0, i1 - numbers of  */
  /* the first and the last function to evaluate, n - number of (domain)     */
  /* points, tkn - the points (ascending), fkn, lkn - arrays in which for    */
  /* each function the numbers of the first and the last point, at which the */
  /* function is nonzero, are stored, b, bt, btt, bttt - arrays in which the */
  /* values of the functions and their derivatives are to be stored          */
  
  void   *sp;
  int    i, j, l, fk, lk;
  double *p;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMemd ( lastknot-deg );
  if ( !p )
    goto failure;

  memset ( b, 0, (i1-i0+1)*n*sizeof(double) );
  memset ( bt, 0, (i1-i0+1)*n*sizeof(double) );
  memset ( btt, 0, (i1-i0+1)*n*sizeof(double) );
  memset ( bttt, 0, (i1-i0+1)*n*sizeof(double) );

        /* the i-th function is nonzero for the arguments between knot[i] */
        /* and knot[i+deg+1] */
  for ( fk = 0; tkn[fk] < knots[i0]; fk ++ ) ;
  for ( lk = fk+1; lk+1 < n && tkn[lk+1] <= knots[i0+deg+1]; lk ++ ) ;
  memset ( p, 0, (lastknot-deg)*sizeof(double) );
  for ( i = i0, l = 0;  i <= i1;  i ++, l ++ ) {
    p[i] = 1.0;
    while ( tkn[fk] < knots[i] ) fk ++;
    while ( lk+1 < n && tkn[lk+1] <= knots[i+deg+1] ) lk ++;
    fkn[l] = fk;
    lkn[l] = lk;
    for ( j = fk; j <= lk; j++ )
      mbs_deBoorDer3C1d ( deg, lastknot, knots, p, tkn[j],
                          &b[l*n+j], &bt[l*n+j], &btt[l*n+j], &bttt[l*n+j] );
    p[i] = 0.0;
 }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabBSFuncDer3d*/

static boolean _g1h_Q2TabSplD1LapGradd ( int nkn, const double *tkn,
                   const double *hfunc, const double *dhfunc,
                   const double *ddhfunc, const double *dddhfunc,
                   int nzc, int i0, int i1,
                   int lastomcknot, const double *omcknots, const double *fcomc,
                   int lastpvknot, const double *pvknots, const double *pv,
                   const double *trd,
                   vector2d *lgr )
{
  void   *sp;
  double *fp, bp[9];
  int    i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 8*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 8*nkn*sizeof(double) );
  ii = 8*i0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer3d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
           i1-i0+1, &tkn[i0], 8, &fp[ii], &fp[ii+1], &fp[ii+2], &fp[ii+3] );
  mbs_TabBSCurveDer3d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pv,
         i1-i0+1, &tkn[i0], 8, &fp[ii+4], &fp[ii+5], &fp[ii+6], &fp[ii+7] );

  memset ( lgr, 0, nkn*nkn*sizeof(vector2d) );
  for ( i = i0;  i <= i1;  i++, ii += 8 ) {
    for ( j = jj = 0, k = nkn*i, kk = 18*k;
          j < nkn;
          j++, jj += 4, k++, kk += 18 ) {
      switch ( nzc ) {
    default:
        bp[0] = fp[ii+1]*hfunc[jj]   + fp[ii+5]*hfunc[jj+2];
        bp[2] = fp[ii+2]*hfunc[jj]   + fp[ii+6]*hfunc[jj+2];
        bp[5] = fp[ii+3]*hfunc[jj]   + fp[ii+7]*hfunc[jj+2];
        bp[1] = fp[ii]*dhfunc[jj]    + fp[ii+4]*dhfunc[jj+2];
        bp[3] = fp[ii+1]*dhfunc[jj]  + fp[ii+5]*dhfunc[jj+2];
        bp[6] = fp[ii+2]*dhfunc[jj]  + fp[ii+6]*dhfunc[jj+2];
        bp[4] = fp[ii]*ddhfunc[jj]   + fp[ii+4]*ddhfunc[jj+2];
        bp[7] = fp[ii+1]*ddhfunc[jj] + fp[ii+5]*ddhfunc[jj+2];
        bp[8] = fp[ii]*dddhfunc[jj]  + fp[ii+4]*dddhfunc[jj+2];
        break;
    case 1:
        bp[0] = fp[ii+5]*hfunc[jj+2];
        bp[2] = fp[ii+6]*hfunc[jj+2];
        bp[5] = fp[ii+7]*hfunc[jj+2];
        bp[1] = fp[ii+4]*dhfunc[jj+2];
        bp[3] = fp[ii+5]*dhfunc[jj+2];
        bp[6] = fp[ii+6]*dhfunc[jj+2];
        bp[4] = fp[ii+4]*ddhfunc[jj+2];
        bp[7] = fp[ii+5]*ddhfunc[jj+2];
        bp[8] = fp[ii+4]*dddhfunc[jj+2];
        break;
      }
      lgr[k].x = pkn_ScalarProductd ( 9, &trd[kk], bp );
      lgr[k].y = pkn_ScalarProductd ( 9, &trd[kk+9], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabSplD1LapGradd*/

static boolean _g1h_Q2TabSplD2LapGradd ( int nkn, const double *tkn,
                   const double *hfunc, const double *dhfunc,
                   const double *ddhfunc, const double *dddhfunc,
                   int nzc, int j0, int j1,
                   int lastomcknot, const double *omcknots, const double *fcomc,
                   int lastpuknot, const double *puknots, const double *pu,
                   const double *trd,
                   vector2d *lgr )
{
  void   *sp;
  double *fp, bp[9];
  int    i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 8*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 8*nkn*sizeof(double) );
  jj = 8*j0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer3d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
           j1-j0+1, &tkn[j0], 8, &fp[jj], &fp[jj+1], &fp[jj+2], &fp[jj+3] );
  mbs_TabBSCurveDer3d ( 1, G1_CROSS01DEG, lastpuknot, puknots, pu,
         j1-j0+1, &tkn[j0], 8, &fp[jj+4], &fp[jj+5], &fp[jj+6], &fp[jj+7] );

  memset ( lgr, 0, nkn*nkn*sizeof(vector2d) );
  for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
    for ( j = j0, jj = 8*j0, k = nkn*i+j0, kk = 18*k;
          j <= j1;
          j++, jj += 8, k++, kk += 18 ) {
      switch ( nzc ) {
    default:
        bp[1] = fp[jj+1]*hfunc[ii]   + fp[jj+5]*hfunc[ii+2];
        bp[4] = fp[jj+2]*hfunc[ii]   + fp[jj+6]*hfunc[ii+2];
        bp[8] = fp[jj+3]*hfunc[ii]   + fp[jj+7]*hfunc[ii+2];
        bp[0] = fp[jj]*dhfunc[ii]    + fp[jj+4]*dhfunc[ii+2];
        bp[3] = fp[jj+1]*dhfunc[ii]  + fp[jj+5]*dhfunc[ii+2];
        bp[7] = fp[jj+2]*dhfunc[ii]  + fp[jj+6]*dhfunc[ii+2];
        bp[2] = fp[jj]*ddhfunc[ii]   + fp[jj+4]*ddhfunc[ii+2];
        bp[6] = fp[jj+1]*ddhfunc[ii] + fp[jj+5]*ddhfunc[ii+2];
        bp[5] = fp[jj]*dddhfunc[ii]  + fp[jj+4]*dddhfunc[ii+2];
        break;
    case 1:
        bp[1] = fp[jj+5]*hfunc[ii+2];
        bp[4] = fp[jj+6]*hfunc[ii+2];
        bp[8] = fp[jj+7]*hfunc[ii+2];
        bp[0] = fp[jj+4]*dhfunc[ii+2];
        bp[3] = fp[jj+5]*dhfunc[ii+2];
        bp[7] = fp[jj+6]*dhfunc[ii+2];
        bp[2] = fp[jj+4]*ddhfunc[ii+2];
        bp[6] = fp[jj+5]*ddhfunc[ii+2];
        bp[5] = fp[jj+4]*dddhfunc[ii+2];
        break;
      }
      lgr[k].x = pkn_ScalarProductd ( 9, &trd[kk], bp );
      lgr[k].y = pkn_ScalarProductd ( 9, &trd[kk+9], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabSplD2LapGradd*/

static boolean _g1h_Q2TabCurveLapCoeff0d ( const vector2d *di,
                                           int nkn, const double *tkn, double v,
                                           double *trd )
{
  void     *sp;
  vector2d *p, *pv, *pvv;
  vector2d d, du, dv, duu, duv, dvv, gx, gy, gxx, gxy, gyy;
  double   A21[6], A22[9];
  int      i, j;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMem ( 3*(G1H_FINALDEG+1)*sizeof(vector2d) );
  if ( !p )
    goto failure;
  pv = &p[G1H_FINALDEG+1];  pvv = &pv[G1H_FINALDEG+1];
  if ( !mbs_multiBCHornerDer2d ( G1H_FINALDEG, G1H_FINALDEG+1, 2, 2*(G1H_FINALDEG+1),
                           (double*)di, v, (double*)p, (double*)pv, (double*)pvv ) )
    goto failure;
  for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
    if ( !mbs_BCHornerDer2C2d ( G1H_FINALDEG, p, tkn[i], &d, &du, &duu ) )
      goto failure;
    if ( !mbs_BCHornerDerC2d ( G1H_FINALDEG, pv, tkn[i], &dv, &duv ) )
      goto failure;
    if ( !mbs_BCHornerC2d ( G1H_FINALDEG, pvv, tkn[i], &dvv ) )
      goto failure;
    if ( !pkn_f2iDerivatives2d ( du.x, du.y, dv.x, dv.y, duu.x, duu.y, duv.x, duv.y,
                                 dvv.x, dvv.y, &gx.x, &gy.x, &gxx.x, &gxy.x, &gyy.x ) )
      goto failure;
    pkn_Setup2DerA21Matrixd ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
    pkn_Setup2DerA22Matrixd ( gx.x, gx.y, gy.x, gy.y, A22 );
    trd[j+0] = A21[0]+A21[4];  trd[j+1] = A21[1]+A21[5];
    trd[j+2] = A22[0]+A22[6];  trd[j+3] = A22[1]+A22[7];  trd[j+4] = A22[2]+A22[8];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabCurveLapCoeff0d*/

static boolean _g1h_Q2TabCurveLapCoeff1d ( const vector2d *di,
                                           double u, int nkn, const double *tkn,
                                           double *trd )
{
  void     *sp;
  vector2d *p, *pu, *puu;
  vector2d d, du, dv, duu, duv, dvv, gx, gy, gxx, gxy, gyy;
  double   A21[6], A22[9];
  int      i, j;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMem ( 3*(G1H_FINALDEG+1)*sizeof(vector2d) );
  if ( !p )
    goto failure;
  pu = &p[G1H_FINALDEG+1];  puu = &pu[G1H_FINALDEG+1];
  if ( !mbs_multiBCHornerDer2d ( G1H_FINALDEG, 1, 2*(G1H_FINALDEG+1), 0,
                           (double*)di, u, (double*)p, (double*)pu, (double*)puu ) )
    goto failure;
  for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
    if ( !mbs_BCHornerDer2C2d ( G1H_FINALDEG, p, tkn[i], &d, &dv, &dvv ) )
      goto failure;
    if ( !mbs_BCHornerDerC2d ( G1H_FINALDEG, pu, tkn[i], &du, &duv ) )
      goto failure;
    if ( !mbs_BCHornerC2d ( G1H_FINALDEG, puu, tkn[i], &duu ) )
      goto failure;
    if ( !pkn_f2iDerivatives2d ( du.x, du.y, dv.x, dv.y, duu.x, duu.y, duv.x, duv.y,
                                 dvv.x, dvv.y, &gx.x, &gy.x, &gxx.x, &gxy.x, &gyy.x ) )
      goto failure;
    pkn_Setup2DerA21Matrixd ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
    pkn_Setup2DerA22Matrixd ( gx.x, gx.y, gy.x, gy.y, A22 );
    trd[j+0] = A21[0]+A21[4];  trd[j+1] = A21[1]+A21[5];
    trd[j+2] = A22[0]+A22[6];  trd[j+3] = A22[1]+A22[7];  trd[j+4] = A22[2]+A22[8];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabCurveLapCoeff1d*/

boolean _g1h_TabBSFuncDer2Jd ( int rr, int deg, int nk, int m2,
             int lastcknot, const double *cknots,
             double *atbs, double *atbst, double *atbstt0, double *atbstt1 )
{
  void   *sp;
  int    i, j, k, l;
  double *p;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMemd ( lastcknot-deg );
  if ( !p )
    goto failure;

  memset ( atbs, 0, (nk+2)*rr*sizeof(double) );
  memset ( atbst, 0, (nk+2)*rr*sizeof(double) );
  memset ( atbstt0, 0, (nk+2)*rr*sizeof(double) );
  memset ( atbstt1, 0, (nk+2)*rr*sizeof(double) );

  memset ( p, 0, (lastcknot-deg)*sizeof(double) );
  for ( i = 2; i < lastcknot-deg-2; i++ ) {
    p[i] = 1.0;
    for ( j = 0; j <= nk+1; j++ ) {
      k = m2*j;
      l = (i-2)*(nk+2)+j;
      if ( j > 0 )
        mbs_deBoorDer2C1d ( deg, 2*deg+1, &cknots[k-m2], &p[k-m2],
            cknots[deg-m2+1+k], &atbs[l], &atbst[l], &atbstt0[l] );
      if ( j <= nk )
        mbs_deBoorDer2C1d ( deg, 2*deg+1, &cknots[k], &p[k],
            cknots[deg-m2+1+k], &atbs[l], &atbst[l], &atbstt1[l] );
    }
    p[i] = 0.0;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabBSFuncDer2Jd*/

static void _g1h_TabSplCLapJumpd ( int njcurves, int nk, int m2, int nquad,
                   int fknj, int lknj,
                   const double *tbsj, const double *tbstj, const double *tbsttj,
                   const double *atbsj, const double *atbstj,
                   const double *atbstt0j, const double *atbstt1j,
                   int fknk, int lknk,
                   const double *tbsk, const double *tbstk, const double *tbsttk,
                   const double *atbsk, const double *atbstk,
                   const double *atbstt0k, const double *atbstt1k,
                   const double *trd, double *lapjc )
{
  int    i, j, k;
  double der0[5], der1[5];  /* pu, pv, puu, puv, pvv */

  memset ( lapjc, 0, njcurves*nquad*sizeof(double) );
        /* curve v == 0 */
  for ( i = fknj; i <= lknj; i++ ) {
    _g1h_TensDer2d ( tbsj[i], tbstj[i], tbsttj[i],
                     atbsk[0], atbstk[0], atbstt1k[0], der0 );
    lapjc[i] = pkn_ScalarProductd ( 5, der0, &trd[5*i] );
  }
        /* curve v == 1 */
  k = 2*nquad;
  for ( i = fknj; i <= lknj; i++ ) {
    _g1h_TensDer2d ( tbsj[i], tbstj[i], tbsttj[i],
                     atbsk[nk+1], atbstk[nk+1], atbstt0k[nk+1], der0 );
    lapjc[nquad+i] = -pkn_ScalarProductd ( 5, der0, &trd[5*(k+i)] );
  }

        /* curve u == 1*/
  k = 4*nquad;
  for ( i = fknk; i <= lknk; i++ ) {
    _g1h_TensDer2d ( atbsj[nk+1], atbstj[nk+1], atbstt0j[nk+1],
                     tbsk[i], tbstk[i], tbsttk[i], der0 );
    lapjc[2*nquad+i] = -pkn_ScalarProductd ( 5, der0, &trd[5*(k+i)] );
  }

  if ( G1H_FINALDEG-m2 == 1 ) {
        /* curves v == const */
    for ( j = 0; j < nk; j++ ) {
      k += nquad;
      for ( i = fknj; i <= lknj; i++ ) {
        _g1h_TensDer2d ( tbsj[i], tbstj[i], tbsttj[i],
                         atbsk[j+1], atbstk[j+1], atbstt0k[j+1], der0 );
        _g1h_TensDer2d ( tbsj[i], tbstj[i], tbsttj[i],
                         atbsk[j+1], atbstk[j+1], atbstt1k[j+1], der1 );
        pkn_SubtractMatrixd ( 1, 5, 0, der1, 0, der0, 0, der0 );
        lapjc[k+i] = pkn_ScalarProductd ( 5, der0, &trd[5*(k+i)] );
      }
    }
        /* curves u == const */
    for ( j = 0; j < nk; j++ ) {
      k +=  nquad;
      for ( i = fknk; i <= lknk; i++ ) {
        _g1h_TensDer2d ( atbsj[j+1], atbstj[j+1], atbstt0j[j+1],
                         tbsk[i], tbstk[i], tbsttk[i], der0 );
        _g1h_TensDer2d ( atbsj[j+1], atbstj[j+1], atbstt1j[j+1],
                         tbsk[i], tbstk[i], tbsttk[i], der0 );
        pkn_SubtractMatrixd ( 1, 5, 0, der1, 0, der0, 0, der0 );
        lapjc[k+i] = pkn_ScalarProductd ( 5, der0, &trd[5*(k+i)] );
      }
    }
  }
} /*_g1h_TabSplCLapJumpd*/

static void _g1h_TabSplCLapJump0d ( int nquad, int fkn, int lkn,
                                    const double *tbs, const double *tbst,
                                    const double *tbstt,
                                    double atbstt, const double *trd,
                                    double *lapjc )
{
  int    i, k;
  double der[5];

        /* curve u == 0 */
  memset ( lapjc, 0, nquad*sizeof(double) );
  for ( i = fkn, k = 5*fkn;  i <= lkn;  i++, k += 5 ) {
    _g1h_TensDer2d ( 0.0, 0.0, atbstt, tbs[i], tbst[i], tbstt[i], der );
    lapjc[i] = -pkn_ScalarProductd ( 5, der, &trd[k] );
  }
} /*_g1h_TabSplCLapJump0d*/

static void _g1h_TabSplCLapJump1d ( int nquad, int fkn, int lkn,
                   const double *tbs, const double *tbst, const double *tbstt,
                   double atbsttj,
                   const double *trd, double *lapjc )
{
  int    i, k;
  double der[5];

  memset ( lapjc, 0, nquad*sizeof(double) );
  for ( i = fkn, k = 5*fkn;  i <= lkn;  i++, k += 5 ) {
    _g1h_TensDer2d ( tbs[i], tbst[i], tbstt[i], 0.0, 0.0, atbsttj, der );
    lapjc[i] = pkn_ScalarProductd ( 5, der, &trd[k] );
  }
} /*_g1h_TabSplCLapJump1d*/

static void _g1h_TabSplCLapJump2d ( int nquad, double atbsttj, int fkn, int lkn,
                   const double *tbs, const double *tbst, const double *tbstt,
                   const double *trd, double *lapjc )
{
  int    i, k;
  double der[5];

  memset ( lapjc, 0, nquad*sizeof(double) );
  for ( i = fkn, k = 5*fkn;  i <= lkn;  i++, k += 5 ) {
    _g1h_TensDer2d ( 0.0, 0.0, atbsttj, tbs[i], tbst[i], tbstt[i], der );
    lapjc[i] = pkn_ScalarProductd ( 5, der, &trd[k] );
  }
} /*_g1h_TabSplCLapJump2d*/

static boolean _g1h_TabSplD1LapJumpd ( int njcurves, int nk, int m1,
             int nquad, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *ahfunc, const double *adhfunc, const double *addhfunc,
             int nzc, int i0, int i1,
             int lastomcknot, const double *omcknots, const double *fcomc,
             int lastpvknot, const double *pvknots,
             const double *pv, const double *pu,
             const double *trd, const double *trdii, const double *trdin,
             double *lapjd )
{
  void   *sp;
  double *fp, bp[5], f00, f01, f02, f02a, f10, f11, f12, f12a;
  int    i, ii, k, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 9*nquad );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 9*nquad*sizeof(double) );
  ii = 9*i0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer2d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            i1-i0+1, &tkn[i0], 9, &fp[ii], &fp[ii+1], &fp[ii+2] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pv,
          i1-i0+1, &tkn[i0], 9, &fp[ii+3], &fp[ii+4], &fp[ii+5] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pu,
          i1-i0+1, &tkn[i0], 9, &fp[ii+6], &fp[ii+7], &fp[ii+8] );

  memset ( lapjd, 0, njcurves*nquad*sizeof(double) );
        /* curve v == 0 */
  for ( i = i0;  i <= i1;  i++, ii += 9 ) {
    switch ( nzc ) {
  default:
      bp[0] = fp[ii+1]*ahfunc[0]  + fp[ii+4]*ahfunc[2];
      bp[1] = fp[ii]*adhfunc[0]   + fp[ii+3]*adhfunc[2];
      bp[2] = fp[ii+2]*ahfunc[0]  + fp[ii+5]*ahfunc[2];
      bp[3] = fp[ii+1]*adhfunc[0] + fp[ii+4]*adhfunc[2];
      bp[4] = fp[ii]*addhfunc[0]  + fp[ii+3]*addhfunc[2];
      break;
  case 1:
      bp[0] = fp[ii+4]*ahfunc[2];
      bp[1] = fp[ii+3]*adhfunc[2];
      bp[2] = fp[ii+5]*ahfunc[2];
      bp[3] = fp[ii+4]*adhfunc[2];
      bp[4] = fp[ii+3]*addhfunc[2];
      break;
    }
    lapjd[i] = pkn_ScalarProductd ( 5, &trd[5*i], bp );
    switch ( nzc ) {
  default:
      bp[0] = fp[ii]*adhfunc[0]   + fp[ii+6]*adhfunc[2];
      bp[1] = fp[ii+1]*ahfunc[0]  + fp[ii+7]*ahfunc[2];
      bp[2] = fp[ii]*addhfunc[0]  + fp[ii+6]*addhfunc[2];
      bp[3] = fp[ii+1]*adhfunc[0] + fp[ii+7]*adhfunc[2];
      bp[4] = fp[ii+2]*ahfunc[0]  + fp[ii+8]*ahfunc[2];
      break;
  case 1:
      bp[0] = fp[ii+6]*adhfunc[2];
      bp[1] = fp[ii+7]*ahfunc[2];
      bp[2] = fp[ii+6]*addhfunc[2];
      bp[3] = fp[ii+7]*adhfunc[2];
      bp[4] = fp[ii+8]*ahfunc[2];
      break;
    }
    lapjd[i] -= pkn_ScalarProductd ( 5, &trdii[5*i], bp );
  }
        /* curve v == 1 */
  for ( i = i0, ii = 9*i0;  i <= i1;  i++, ii += 9 ) {
    switch ( nzc ) {
  default:
      bp[0] = fp[ii+1]*ahfunc[4]  + fp[ii+4]*ahfunc[6];
      bp[1] = fp[ii]*adhfunc[4]   + fp[ii+3]*adhfunc[6];
      bp[2] = fp[ii+2]*ahfunc[4]  + fp[ii+5]*ahfunc[6];
      bp[3] = fp[ii+1]*adhfunc[4] + fp[ii+4]*adhfunc[6];
      bp[4] = fp[ii]*addhfunc[4]  + fp[ii+3]*addhfunc[6];
      break;
  case 1:
      bp[0] = fp[ii+4]*ahfunc[6];
      bp[1] = fp[ii+3]*adhfunc[6];
      bp[2] = fp[ii+5]*ahfunc[6];
      bp[3] = fp[ii+4]*adhfunc[6];
      bp[4] = fp[ii+3]*addhfunc[6];
      break;
    }
    lapjd[nquad+i] = -pkn_ScalarProductd ( 5, &trd[5*(2*nquad+i)], bp );
  }
        /* curve u == 1 */
  if ( fcomc[lastomcknot-G1_CROSS00DEG-3] || pv[lastpvknot-G1_CROSS01DEG-3] ) {
    if ( nzc == 0 )
      mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc, 1.0,
                          &f00, &f01, &f02 );
    else
      f00 = f01 = f02 = 0.0;
    mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv, 1.0,
                        &f10, &f11, &f12 );
    for ( i = ii = 0;  i < nquad;  i++, ii += 4 ) {
      switch ( nzc ) {
    default:
        bp[0] = f01*hfunc[ii]   + f11*hfunc[ii+2];
        bp[1] = f00*dhfunc[ii]  + f10*dhfunc[ii+2];
        bp[2] = f02*hfunc[ii]   + f12*hfunc[ii+2];
        bp[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
        bp[4] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
        break;
    case 1:
        bp[0] = f11*hfunc[ii+2];
        bp[1] = f10*dhfunc[ii+2];
        bp[2] = f12*hfunc[ii+2];
        bp[3] = f11*dhfunc[ii+2];
        bp[4] = f10*ddhfunc[ii+2];
        break;
      }
      lapjd[2*nquad+i] = -pkn_ScalarProductd ( 5, &trd[5*(4*nquad+i)], bp );
    }
  }
        /* "inner" curves u == const */
  if ( G1_CROSS00DEG-G1_BF01DEG-m1 < 1 ) {  /* there are nonzero jumps */
    for ( k = 0; k < nk; k++ ) {
      kk = G1_CROSS00DEG-1+(k+1)*m1;
      if ( nzc == 0 ) {
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                            omcknots[kk], &f00, &f01, &f02 );
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, kk+G1_CROSS00DEG, omcknots, fcomc,
                            omcknots[kk],
                            &f00, &f01, &f02a );
        f02 -= f02a;
      }
      else
        f02 = 0.0;
      kk = G1_CROSS01DEG-3+(k+1)*(m1+G1_BF01DEG);
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv,
                          pvknots[kk], &f10, &f11, &f12 );
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, kk+G1_CROSS01DEG, pvknots, pv,
                          pvknots[kk], &f10, &f11, &f12a );
      f12 -= f12a;
      for ( i = ii = 0;  i < nquad;  i++, ii += 4 ) {
        switch ( nzc ) {
      default:
          bp[2] = f02*hfunc[ii] + f12*hfunc[ii+2];
          break;
      case 1:
          bp[2] = f12*hfunc[ii+2];
          break;
        }
        lapjd[(3+nk+k)*nquad+i] = trdin[5*((nk+k)*nquad+i)+2]*bp[2];
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabSplD1LapJumpd*/

static boolean _g1h_TabSplD2LapJumpd ( int njcurves, int nk, int m1,
             int nquad, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *ahfunc, const double *adhfunc, const double *addhfunc,
             int nzc, int i0, int i1,
             int lastomcknot, const double *omcknots, const double *fcomc,
             int lastpvknot, const double *pvknots, const double *pu,
             const double *trd, const double *trdin, double *lapjd )
{
  void   *sp;
  double *fp, bp[5], f00, f01, f02, f02a, f10, f11, f12, f12a;
  int    i, ii, k, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 6*nquad );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 6*nquad*sizeof(double) );
  ii = 6*i0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer2d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            i1-i0+1, &tkn[i0], 6, &fp[ii], &fp[ii+1], &fp[ii+2] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pu,
          i1-i0+1, &tkn[i0], 6, &fp[ii+3], &fp[ii+4], &fp[ii+5] );

  memset ( lapjd, 0, njcurves*nquad*sizeof(double) );
        /* curve v == 0 */
  if ( nzc == 0 )
    mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                        0.0, &f00, &f01, &f02 );
  else
    f00 = f01 = f02 = 0.0;
  mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                      0.0, &f10, &f11, &f12 );
  for ( i = ii = kk = 0;  i < nquad;  i++, ii += 4, kk += 5 ) {
    switch ( nzc ) {
  default:
      bp[0] = f00*dhfunc[ii]  + f10*dhfunc[ii+2];
      bp[1] = f01*hfunc[ii]   + f11*hfunc[ii+2];
      bp[2] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
      bp[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
      bp[4] = f02*hfunc[ii]   + f12*hfunc[ii+2];
      break;
  case 1:
      bp[0] = f10*dhfunc[ii+2];
      bp[1] = f11*hfunc[ii+2];
      bp[2] = f10*ddhfunc[ii+2];
      bp[3] = f11*dhfunc[ii+2];
      bp[4] = f12*hfunc[ii+2];
      break;
    }
    lapjd[i] = pkn_ScalarProductd ( 5, bp, &trd[kk] );
  }
        /* curve v == 1 */
  if ( nzc == 0 )
    mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                        1.0, &f00, &f01, &f02 );
  else
    f00 = f01 = f02 = 0.0;
  mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                      1.0, &f10, &f11, &f12 );
  for ( i = ii = 0, kk = 10*nquad;  i < nquad;  i++, ii += 4, kk += 5 ) {
    switch ( nzc ) {
  default:
      bp[0] = f00*dhfunc[ii]  + f10*dhfunc[ii+2];
      bp[1] = f01*hfunc[ii]   + f11*hfunc[ii+2];
      bp[2] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
      bp[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
      bp[4] = f02*hfunc[ii]   + f12*hfunc[ii+2];
      break;
  case 1:
      bp[0] = f10*dhfunc[ii+2];
      bp[1] = f11*hfunc[ii+2];
      bp[2] = f10*ddhfunc[ii+2];
      bp[3] = f11*dhfunc[ii+2];
      bp[4] = f12*hfunc[ii+2];
      break;
    }
    lapjd[nquad+i] = -pkn_ScalarProductd ( 5, bp, &trd[kk] );
  }
        /* curve u == 1 */
  for ( i = i0, ii = 6*i0, kk = 5*(4*nquad+i0);
        i <= i1;
        i++, ii += 6, kk += 5 ) {
    switch ( nzc ) {
  default:
      bp[0] = fp[ii]*adhfunc[4]   + fp[ii+3]*adhfunc[6];
      bp[1] = fp[ii+1]*ahfunc[4]  + fp[ii+4]*ahfunc[6];
      bp[2] = fp[ii]*addhfunc[4]  + fp[ii+3]*addhfunc[6];
      bp[3] = fp[ii+1]*adhfunc[4] + fp[ii+4]*adhfunc[6];
      bp[4] = fp[ii+2]*ahfunc[4]  + fp[ii+5]*ahfunc[6];
      break;
  case 1:
      bp[0] = fp[ii+3]*adhfunc[6];
      bp[1] = fp[ii+4]*ahfunc[6];
      bp[2] = fp[ii+3]*addhfunc[6];
      bp[3] = fp[ii+4]*adhfunc[6];
      bp[4] = fp[ii+5]*ahfunc[6];
      break;
    }
    lapjd[2*nquad+i] = -pkn_ScalarProductd ( 5, bp, &trd[kk] );
  }
        /* "inner" curves v == const */
  if ( G1_CROSS00DEG-G1_BF01DEG-m1 < 1 ) {  /* there are nonzero jumps */
    for ( k = 0; k < nk; k++ ) {
      kk = G1_CROSS00DEG-1+(k+1)*m1;
      if ( nzc == 0 ) {
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                            omcknots[kk], &f00, &f01, &f02 );
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, kk+G1_CROSS00DEG, omcknots, fcomc,
                            omcknots[kk], &f00, &f01, &f02a );
        f02 -= f02a;
      }
      else
        f02 = 0.0;
      kk = G1_CROSS01DEG-3+(k+1)*(m1+G1_BF01DEG);
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                          pvknots[kk], &f10, &f11, &f12 );
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, kk+G1_CROSS01DEG, pvknots, pu,
                          pvknots[kk], &f10, &f11, &f12a );
      f12 -= f12a;
      for ( i = ii = 0;  i < nquad;  i++, ii += 4 ) {
        switch ( nzc ) {
      default:
          bp[4] = f02*hfunc[ii] + f12*hfunc[ii+2];
          break;
      case 1:
          bp[4] = f12*hfunc[ii+2];
          break;
        }
        lapjd[(3+k)*nquad+i] = trdin[5*(k*nquad+i)+4]*bp[4];
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabSplD2LapJumpd*/

static void _g1h_TabSplD3LapJumpd ( int nquad,
              const double *hfunc, const double *dhfunc, const double *ddhfunc,
              int nzc,
              int lastomcknot, const double *omcknots, const double *fcomc,
              int lastpvknot, const double *pvknots, const double *pv,
              const double *trd,
              double *lapjd )
{
  double bp[5], f00, f01, f02, f10, f11, f12;
  int    i, ii, kk;

  if ( nzc == 0 )
    mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                        0.0, &f00, &f01, &f02 );
  else
    f00 = f01 = f02 = 0.0;
  mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv,
                      0.0, &f10, &f11, &f12 );
  memset ( lapjd, 0, nquad*sizeof(double) );
  for ( i = ii = kk = 0;  i < nquad;  i++, ii += 4, kk += 5 ) {
    switch ( nzc ) {
  default:
      bp[0] = f01*hfunc[ii]   + f11*hfunc[ii+2];
      bp[1] = f00*dhfunc[ii]  + f10*dhfunc[ii+2];
      bp[2] = f02*hfunc[ii]   + f12*hfunc[ii+2];
      bp[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
      bp[4] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
      break;
  case 1:
      bp[0] = f11*hfunc[ii+2];
      bp[1] = f10*dhfunc[ii+2];
      bp[2] = f12*hfunc[ii+2];
      bp[3] = f11*dhfunc[ii+2];
      bp[4] = f10*ddhfunc[ii+2];
      break;
    }
    lapjd[i] = -pkn_ScalarProductd ( 5, bp, &trd[kk] );
  }
} /*_g1h_TabSplD3LapJumpd*/

static double _g1h_Q2SplIntegral0d ( int nquad, int ncurv,
                                     double *jac, double *lapj1, double *lapj2 )
{
  int    i;
  double s;

  for ( i = 0, s = 0.0;  i < ncurv*nquad;  i++ )
    s += lapj1[i]*lapj2[i]*jac[i];
  s /= (double)nquad;
  return (double)s;
} /*_g1h_Q2SplIntegral0d*/

boolean g1h_Q2SplComputeFormMatrixd ( GHoleDomaind *domain )
{
  void *sp;
  GHolePrivateRecd    *privateG;
  G1HolePrivateRecd   *privateG1;
  G1HoleSPrivateRecd  *privateS;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int      lastomcknot, lastpvknot, lastcknot;
  double   *omcknots, *pvknots, *cknots;
  int      nk, m1, m2, nquad, nquadsq, rr, r, R, s, S, m, nab, nabc;
  int      i, ii, j, jj, k, l, l1, fn;
  int      i00, i01, i10, i11, j00, j01, j10, j11, nzc;
  int      asize, bsize;
  double   *amat, *bmat, *Aij, *Aik;
  double   *tkn, *atkn, *jac, *trd, *trdin, *fcomc, *pu, *pv;
  double   *tbs, *tbst, *tbstt, *tbsttt, *atbs, *atbst, *atbstt0, *atbstt1;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *ec00, *ec01, *ec10, *ec11, *ed00, *ed01, *ed10, *ed11;
  double   *hfunc, *dhfunc, *ddhfunc, *dddhfunc, *ahfunc, *adhfunc, *addhfunc;
  double   x, qc, C;
  vector2d *lgr;
  int      *fkn, *lkn;
  unsigned short *support_b, supp, mask;
  int      option, ndata, *idata;
  double   *fdata;
  int      njcurves, bstsize, degu, degv;
  point2d  *di, *cdi;
  double   *lapja, *lapjb, *lapjc, *lapjd, *eicp1, *eicp2;
  boolean  jumpC, jumpD;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  privateS  = domain->SprivateG1;
  if ( !privateS )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = privateS->nsfunc_c;
  nfunc_d = privateS->nsfunc_d;
  support_b = privateG->support_b;
  nk = privateS->nk;
  m1 = privateS->m1;
  m2 = privateS->m2;

  cknots = privateS->cknots;
  lastcknot = privateS->lastcknot;
  omcknots = privateS->omcknots;
  lastomcknot = privateS->lastomcknot;
  pvknots = privateS->pvknots;
  lastpvknot = privateS->lastpvknot;

        /* compute the size of the arrays and allocate memory */
  rr = G1H_FINALDEG-3+nk*m2;
  r = rr*rr;    /* == nfunc_c/hole_k */
  s = 2*nk*m1;  /* == nfunc_d/hole_k */
  R = r+hole_k*s;
  S = R+nfunc_a;
  nab = nfunc_a+nfunc_b;
  nabc = nab+r;
  asize = pkn_Block3ArraySize ( hole_k-1, r, S );
  bsize = nfunc_b*(nfunc_a+nfunc_c+nfunc_d);
  amat = privateS->Q2SAMat = malloc ( asize*sizeof(double) );
  bmat = privateS->Q2SBMat = malloc ( bsize*sizeof(double) );
  if ( !amat || !bmat )
    goto failure;
  if ( !privateG->diam )
    privateG->diam = gh_DomainDiamd ( domain );

        /* fix the number of quadrature knots */
  nquad = G1_QUAD_FACTOR*(nk+1);
  nquadsq = nquad*nquad;
  m = rr*nquad;

        /* allocate the scratch memory */
  tkn = pkv_GetScratchMemd ( 17*nquad );
  jac = pkv_GetScratchMemd ( nquadsq );
  trd = pkv_GetScratchMemd ( max ( 18*nquadsq, 30*hole_k*nquad ) );
  fkn = (int*)pkv_GetScratchMem ( 2*rr*sizeof(int) );
  lgr = (vector2d*)pkv_GetScratchMem ( ((nabc+2*s)*nquadsq)*sizeof(vector2d) );
  tbs = pkv_GetScratchMemd ( 4*m );
  fcomc = pkv_GetScratchMemd ( lastomcknot+2*lastpvknot-
                               (G1_CROSS00DEG+2*G1_CROSS01DEG) );
  if ( !tkn || !jac || !trd || !fkn || !lgr || !tbs || !fcomc )
    goto failure;
  hfunc = &tkn[nquad];         dhfunc = &hfunc[4*nquad];
  ddhfunc = &dhfunc[4*nquad];  dddhfunc = &ddhfunc[4*nquad];
  tbst = &tbs[m];              tbstt = &tbst[m];
  tbsttt = &tbstt[m];
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  lkn = &fkn[rr];

        /* integrate the Laplacian gradients */
  memset ( amat, 0, asize*sizeof(double) );
  memset ( bmat, 0, bsize*sizeof(double) );

          /* prepare the evaluation of basis functions */
  _gh_PrepareTabKnotsd ( nquad, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3d ( 0.0, 1.0, nquad, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;
  if ( !_g1h_TabBSFuncDer3d ( G1H_FINALDEG, lastcknot, cknots,
                              2, G1H_FINALDEG+nk*m2-2, nquad, tkn, fkn, lkn,
                              tbs, tbst, tbstt, tbsttt ) )
    goto failure;

  for ( i = 0, ii = 1, mask = 0x0001;
        i < hole_k;
        i++, ii = (ii+1) % hole_k, mask = (unsigned short)(mask << 1) ) {
          /* evaluate the Jacobian and the matrix of derivatives of the */
          /* domain patch at the quadrature knots */
    _g1h_GetDiPatchCurvesd ( domain, i, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    _g1h_Q2TabDiPatchJac3d ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             jac, trd );
          /* evaluate the Laplacian gradient of the A block basis functions */
    for ( j = 0; j < nfunc_a; j++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_Q2TabLaplacianGrad0d ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                  fc00, fc01, fd00, fd01, trd, &lgr[j*nquadsq] );
    }
          /* evaluate the Laplacian gradient of the B block basis functions */
    for ( j = 0; j < nfunc_b; j++ ) {
      _g1h_GetBFBPatchCurvesd ( domain, j, i, &fc00, &fc01, &fc10, &fc11,
                                &fd00, &fd01, &fd10, &fd11 );
      _g1h_Q2TabLaplacianGradd ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                 fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                                 trd, &lgr[(nfunc_a+j)*nquadsq] );
    }
          /* evaluate the Laplacian gradient of the C block basis functions */
    for ( j = l1 = 0;  j < rr;  j++ )
      for ( k = 0;  k < rr;  k++, l1++ )
        _g2h_TabSplCLaplacianGradd ( nquad,
                  fkn[j], lkn[j], &tbs[j*nquad], &tbst[j*nquad],
                  &tbstt[j*nquad], &tbsttt[j*nquad],
                  fkn[k], lkn[k], &tbs[k*nquad], &tbst[k*nquad],
                  &tbstt[k*nquad], &tbsttt[k*nquad],
                  trd, &lgr[(nab+l1)*nquadsq] );
          /* evaluate the Laplacian gradient of the D block basis functions; */
          /* do it only for the functions nonzero in Omega_i; j is the       */
          /* function number in the subblock, k is the complete function     */
          /* number */
    for ( j = 0, k = nab+nfunc_c+i*s;
          j < s;
          j++, k++ ) {
      _g1h_GetSplDBasisCrossDerd ( domain, k, i, fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_Q2TabSplD1LapGradd ( nquad, tkn,
                     hfunc, dhfunc, ddhfunc, dddhfunc,
                     nzc, i00, i01, lastomcknot, omcknots, fcomc,
                     lastpvknot, pvknots, pv,
                     trd, &lgr[(nabc+j)*nquadsq] ) )
        goto failure;
    }
    for ( j = 0, k = nab+nfunc_c+ii*s;
          j < s;
          j++, k++ ) {
      _g1h_GetSplDBasisCrossDerd ( domain, k, ii, fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_Q2TabSplD2LapGradd ( nquad, tkn,
                     hfunc, dhfunc, ddhfunc, dddhfunc,
                     nzc, j00, j01, lastomcknot, omcknots, fcomc,
                     lastpvknot, pvknots, pu,
                     trd, &lgr[(nabc+s+j)*nquadsq] ) )
        goto failure;
    }
          /* block A x A */
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1,r,S,hole_k-1,hole_k-1)];
    for ( j = 0; j < nfunc_a; j++ )
      for ( k = 0; k <= j; k++ )
        Aij[pkn_SymMatIndex(R+j,R+k)] +=
                 _g2h_SplIntegrald ( nquad, jac,
                                     0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                                     0, nquad-1, 0, nquad-1, &lgr[k*nquadsq] );

          /* block A x B */
    Aij = &bmat[(nfunc_c+nfunc_d)*nfunc_b];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0; k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          Aij[l1] += _g2h_SplIntegrald ( nquad, jac,
                         0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                         0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );

          /* block A x C */
    if ( i < hole_k-1 )
      Aij = &amat[pkn_Block3FindBlockPos(hole_k-1,r,S,hole_k-1,i)+
                  (r+nfunc_d)*r];
    else
      Aij = &amat[pkn_Block3FindBlockPos(hole_k-1,r,S,hole_k-1,hole_k-1)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        x = _g2h_SplIntegrald ( nquad, jac,
                      0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                      i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
        if ( i < hole_k-1 )
          Aij[l1] = x;
        else
          Aij[pkn_SymMatIndex(R+j,k)] = x;
      }

          /* block A x D */
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1,r,S,hole_k-1,hole_k-1)];
    for ( j = 0; j < nfunc_a; j++ )
      for ( k = 0; k < s; k++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        l1 = pkn_SymMatIndex ( R+j, r+i*s+k );
        Aij[l1] += _g2h_SplIntegrald ( nquad, jac,
                       0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                       i00, i01, j00, j01, &lgr[(nabc+k)*nquadsq] );
        _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        l1 = pkn_SymMatIndex ( R+j, r+ii*s+k );
        Aij[l1] += _g2h_SplIntegrald ( nquad, jac,
                       0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                       i00, i01, j00, j01, &lgr[(nabc+s+k)*nquadsq] );
      }

          /* block B x C */
    Aij = &bmat[i*r*nfunc_b];
    for ( j = l1 = 0;  j < r;  j++ ) {
      i00 = fkn[j/rr];  i01 = lkn[j/rr];
      j00 = fkn[j%rr];  j01 = lkn[j%rr];
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          Aij[l1] = _g2h_SplIntegrald ( nquad, jac,
                        i00, i01, j00, j01, &lgr[(nab+j)*nquadsq],
                        0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
    }

          /* block B x D */
    Aij = &bmat[(nfunc_c+i*s)*nfunc_b];
    Aik = &bmat[(nfunc_c+ii*s)*nfunc_b];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i10, &i11, &j10, &j11 );
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask ) {
          Aij[l1] += _g2h_SplIntegrald ( nquad, jac,
                         i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                         0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
          Aik[l1] += _g2h_SplIntegrald ( nquad, jac,
                         i10, i11, j10, j11, &lgr[(nabc+s+j)*nquadsq],
                         0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
        }
    }

          /* block C x C */
    Aij = &amat[pkn_Block3FindBlockPos ( hole_k-1, r, S, i, i )];
    for ( j = l1 = 0;  j < r;  j++ ) {
      i00 = fkn[j/rr];  i01 = lkn[j/rr];
      j00 = fkn[j%rr];  j01 = lkn[j%rr];
      for ( k = 0;  k <= j;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        Aij[l1] = _g2h_SplIntegrald ( nquad, jac,
                          i00, i01, j00, j01, &lgr[(nab+j)*nquadsq],
                          i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
      }
    }

          /* block C x D */
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, r, S, hole_k-1, i)];
    for ( j = 0, l1 = (r+i*s)*r; j < s; j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        x = _g2h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
        if ( i < hole_k-1 )
          Aij[l1] = +x;
        else
          Aij[pkn_SymMatIndex(r+(hole_k-1)*s+j,k)] = +x;
      }
    }
    for ( j = 0, l1 = (r+ii*s)*r; j < s; j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        x = _g2h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
        if ( i < hole_k-1 )
          Aij[l1] = +x;
        else
          Aij[pkn_SymMatIndex(r+j,k)] = +x;
      }
    }

          /* block D x D */
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, r, S, hole_k-1, hole_k-1)];
    for ( j = 0; j < s; j++) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0; k <= j; k++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        Aij[pkn_SymMatIndex(r+i*s+j,r+i*s+k)] +=
                _g2h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                            i10, i11, j10, j11, &lgr[(nabc+k)*nquadsq] );
      }
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0; k <= j; k++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        Aij[pkn_SymMatIndex(r+ii*s+j,r+ii*s+k)] +=
                _g2h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                            i10, i11, j10, j11, &lgr[(nabc+s+k)*nquadsq] );
      }
      for ( k = 0; k < s; k++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        Aij[pkn_SymMatIndex(r+ii*s+j,r+i*s+k)] =
                _g2h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                            i10, i11, j10, j11, &lgr[(nabc+k)*nquadsq] );
      }
    }
  } /* for ( i = ... ) */
        /* multiply the sums of function values at the quadrature knots */
        /* by the quadrature coefficient */
  qc = 1.0/(double)nquadsq;
  for ( i = 0; i < asize; i++ )
    amat[i] *= qc;
  for ( i = 0; i < bsize; i++ )
    bmat[i] *= qc;

        /* integrate the Laplacian jumps */
  option = privateG1->GetOption ( domain, G1HQUERY_Q2_FORM_CONSTANT, 0,
                                  &ndata, &idata, &fdata );
  switch ( option ) {
case G1H_DEFAULT:
    privateS->C1s = C = 1.0;
                       /* no idea whether this default value makes any sense */
    break;
case G1H_Q2_USE_SUPPLIED_CONSTANT:
    privateS->C1s = C = *fdata;   
    break;
default:  
    goto failure;
  }
  C *= (double)(5+nk*m2)/privateG->diam;

  di = pkv_GetScratchMem ( (G1H_FINALDEG+1)*(G1H_FINALDEG+2)*sizeof(point2d) );
  if ( !di )
    goto failure;
  cdi = &di[(G1H_FINALDEG+1)*(G1H_FINALDEG+1)];
  eicp1 = (double*)di;  eicp2 = &eicp1[16];
  jumpC = (boolean)(G1H_FINALDEG-m2 < 2);
  jumpD = (boolean)(G1_CROSS00DEG-G1_BF01DEG-m1 < 1);
  njcurves = 3 + 2*nk;  /* number of curves to integrate the Laplacian jumps */
                        /* for one Omega_i */
  if ( jumpC || jumpD ) {
    if ( !(trdin = pkv_GetScratchMemd ( 10*nk*nquad )) )
      goto failure;
  }
  else
    trdin = NULL;
  lapja = &lgr[0].x;
  lapjb = &lapja[3*nquad*nfunc_a];   /* functions of block A and B have */
  lapjc = &lapjb[3*nquad*nfunc_b];   /* nonzero Laplacian jump only at  */
  lapjd = &lapjc[(njcurves*r+rr)*nquad];     /* the boundary of Omega_i */
  atkn = pkv_GetScratchMemd ( 26 );
  if ( !atkn )
    goto failure;
  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc ) )
    goto failure;
  bstsize = rr*(nk+2);
  atbs = pkv_GetScratchMemd ( 4*bstsize );
  if ( !atbs )
    goto failure;
  atbst = &atbs[bstsize];
  atbstt0 = &atbst[bstsize];
  atbstt1 = &atbstt0[bstsize];
  if ( !_g1h_TabBSFuncDer2Jd ( rr, G1H_FINALDEG, nk, m2, lastcknot, cknots,
                               atbs, atbst, atbstt0, atbstt1 ) )
    goto failure;

        /* compute the coefficients for the computation of Laplacians */
        /* this is done in a separate loop for all boundary curves of */
        /* the areas Omega_i */
  for ( i = 0; i < hole_k; i++ ) {
    jj = (i+1) % hole_k;
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    if ( !_g1h_TabCurveLapCoeff0d ( c00, c01, c10, c11, d00, d01, d10, d11,
                            nquad, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc, &trd[30*i*nquad],
                            &trd[(30*i+10)*nquad], &trd[(30*i+5)*nquad],
                            &trd[(30*i+20)*nquad] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchd ( domain, jj, 1, di ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( di, nquad, tkn, &trd[(30*i+15)*nquad] );
    if ( !_gh_FindDomSurrndPatchd ( domain, i, 2, di ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( di, nquad, tkn, &trd[(30*i+25)*nquad] );
  }

  for ( i = 0, mask = 0x0001;
        i < hole_k;
        i++, mask = (unsigned short)(mask << 1) ) {
    ii = (i+hole_k-1) % hole_k;
    jj = (i+1) % hole_k;
          /* compute the curve Jacobians */
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    if ( !_g1h_TabCurveJacobiand ( G1_CROSS00DEG, c00, nquad, tkn, jac ) )
      goto failure;
    if ( !_g1h_TabCurveJacobiand ( G1_CROSS10DEG, c10, nquad, tkn, &jac[nquad] ) )
      goto failure;
    if ( !_g1h_TabCurveJacobiand ( G1_CROSS10DEG, d10, nquad, tkn, &jac[2*nquad] ) )
      goto failure;
    if ( jumpC || jumpD ) {
          /* if necessary, compute the Jacobians for the curves inside */
          /* Omega_i */
      mbs_BezC1CoonsToBezd ( 2,
          G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
          3, (double*)c10, G1_CROSS11DEG, (double*)c11,
          G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
          3, (double*)d10, G1_CROSS11DEG, (double*)d11,
          &degu, &degv, (double*)di );
      for ( k = 0; k < nk; k++ ) {
        x = cknots[G1H_FINALDEG+1+k*m2];
        if ( !mbs_multiBCHornerd ( G1H_FINALDEG, G1H_FINALDEG+1, 2, 2*(G1H_FINALDEG+1),
                                   (double*)di, x, (double*)cdi ) )
          goto failure;
        if ( !_g1h_TabCurveJacobiand ( G1H_FINALDEG, cdi, nquad, tkn,
                                       &jac[(3+k)*nquad] ) )
          goto failure;
        if ( !_g1h_Q2TabCurveLapCoeff0d ( di, nquad, tkn, x,
                                          &trdin[5*k*nquad] ) )
          goto failure;
      }
      for ( k = 0; k < nk; k++ ) {
        x = cknots[G1H_FINALDEG+1+k*m2];
        if ( !mbs_multiBCHornerd ( G1H_FINALDEG, 1, 2*(G1H_FINALDEG+1), 0,
                                   (double*)di, x, (double*)cdi ) )
          goto failure;
        if ( !_g1h_TabCurveJacobiand ( G1H_FINALDEG, cdi, nquad, tkn,
                                       &jac[(3+nk+k)*nquad] ) )
          goto failure;
        if ( !_g1h_Q2TabCurveLapCoeff1d ( di, x, nquad, tkn,
                                          &trdin[5*(nk+k)*nquad] ) )
          goto failure;
      }
    }
          /* evaluate the Laplacian jumps of the basis functions */
          /* block A */
    for ( fn = 0; fn < nfunc_a; fn++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, ii, &ec00, &ec01, &ed00, &ed01 );
      _g1h_GetBFAPatchCurvesd ( domain, fn, i, &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1h_Q2TabLaplacianJump0d ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                          atkn, ahfunc, adhfunc, addhfunc,
                          ec00, ec01, ed00, ed01, &trd[(30*ii+5)*nquad],
                          fc00, fc01, fd00, fd01, &trd[30*i*nquad],
                          &trd[(30*i+10)*nquad], &trd[(30*i+20)*nquad],
                          &lapja[3*fn*nquad], &lapja[(3*fn+1)*nquad],
                          &lapja[(3*fn+2)*nquad] ) )
        goto failure;
    }

          /* block B */
    for ( fn = 0; fn < nfunc_b; fn++ ) {
      supp = _g1h_ExtendSupport ( hole_k, support_b[fn] );
      if ( supp & mask ) {
        _g1h_GetBFBPatchCurvesd ( domain, fn, ii,
                  &ec00, &ec01, &ec10, &ec11, &ed00, &ed01, &ed10, &ed11 );
        _g1h_GetBFBPatchCurvesd ( domain, fn, i,
                  &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        gh_GetDomSurrndBFuncd ( domain, fn, jj, 1, eicp1 );
        gh_GetDomSurrndBFuncd ( domain, fn, i, 2, eicp2 );
        if ( !_g1h_Q2TabLaplacianJumpd ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc,
                            ec00, ec01, ec10, ec11, ed00, ed01, ed10, ed11,
                            &trd[(30*ii+5)*nquad],
                            fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                            &trd[30*i*nquad], &trd[(30*i+10)*nquad],
                            &trd[(30*i+20)*nquad],
                            eicp1, &trd[(30*i+15)*nquad],
                            eicp2, &trd[(30*i+25)*nquad],
                            &lapjb[3*fn*nquad], &lapjb[(3*fn+1)*nquad],
                            &lapjb[(3*fn+2)*nquad] ) )
          goto failure;
      }
      else
        memset ( &lapjb[fn*3*nquad], 0, 3*nquad*sizeof(double) );
    }

          /* block C */
            /* subblock i */
    for ( j = l1 = 0; j < rr; j++ )
      for ( k = 0; k < rr; k++, l1++ )
        _g1h_TabSplCLapJumpd ( njcurves, nk, m2, nquad,
                   fkn[j], lkn[j], &tbs[j*nquad], &tbst[j*nquad], &tbstt[j*nquad],
                   &atbs[j*(nk+2)], &atbst[j*(nk+2)],
                   &atbstt0[j*(nk+2)], &atbstt1[j*(nk+2)],
                   fkn[k], lkn[k], &tbs[k*nquad], &tbst[k*nquad], &tbstt[k*nquad],
                   &atbs[k*(nk+2)], &atbst[k*(nk+2)],
                   &atbstt0[k*(nk+2)], &atbstt1[k*(nk+2)],
                   &trd[30*i*nquad], &lapjc[l1*njcurves*nquad] );
    if ( jumpC ) {
      for ( j = l1 = 0; j < rr; j++ )
        for ( k = 0; k < rr; k++, l1++ ) {
              /* curves v == const */
          for ( l = 0; l < nk; l++ ) {
            _g1h_TabSplCLapJump1d ( nquad,
                   fkn[j], lkn[j], &tbs[j*nquad], &tbst[j*nquad], &tbstt[j*nquad],
                   atbstt1[k*(nk+2)+l+1]-atbstt0[k*(nk+2)+l+1],
                   &trdin[l*5*nquad], &lapjc[(l1*njcurves+3+l)*nquad] );
          }
              /* curves u == const */
          for ( l = 0; l < nk; l++ ) {
            _g1h_TabSplCLapJump2d ( nquad,
                   atbstt1[j*(nk+2)+l+1]-atbstt0[j*(nk+2)+l+1],
                   fkn[k], lkn[k], &tbs[k*nquad], &tbst[k*nquad], &tbstt[k*nquad],
                   &trdin[(nk+l)*5*nquad], &lapjc[(l1*njcurves+3+nk+l)*nquad] );
          }
        }
    }
            /* subblock ii = i-1 */
    for ( j = 0;  j < rr;  j++ )
      _g1h_TabSplCLapJump0d ( nquad, fkn[j], lkn[j],
                              &tbs[j*nquad], &tbst[j*nquad], &tbstt[j*nquad],
                              atbstt1[0], &trd[(30*ii+5)*nquad],
                              &lapjc[(njcurves*r+j)*nquad] );

          /* block D */
    for ( j = 0; j < s; j++ ) {
            /* subblock i */
      _g1h_GetSplDBasisCrossDerd ( domain, nfunc_a+nfunc_b+nfunc_c+i*s+j, i,
                                   fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_TabSplD1LapJumpd ( njcurves, nk, m1, nquad,
                                    tkn, hfunc, dhfunc, ddhfunc,
                                    ahfunc, adhfunc, addhfunc, nzc, i00, i01,
                                    lastomcknot, omcknots, fcomc,
                                    lastpvknot, pvknots, pv, pu,
                                    &trd[30*i*nquad], &trd[(30*ii+5)*nquad],
                                    trdin, &lapjd[njcurves*j*nquad] ) )
        goto failure;

            /* subblock jj == i+1 */
      _g1h_GetSplDBasisCrossDerd ( domain, nfunc_a+nfunc_b+nfunc_c+jj*s+j, jj,
                                   fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, jj*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_TabSplD2LapJumpd ( njcurves, nk, m1, nquad,
                                    tkn, hfunc, dhfunc, ddhfunc,
                                    ahfunc, adhfunc, addhfunc, nzc, j00, j01,
                                    lastomcknot, omcknots, fcomc,
                                    lastpvknot, pvknots, pu,
                                    &trd[30*i*nquad], trdin,
                                    &lapjd[njcurves*(j+s)*nquad] ) )
        goto failure;
    }
    for ( j = 0; j < 2; j++ ) {
            /* subblock ii == i-1 --- only the first two functions are */
            /* relevant here */
      _g1h_GetSplDBasisCrossDerd ( domain, nfunc_a+nfunc_b+nfunc_c+ii*s+j, ii,
                                   fcomc, pv, pu );
      _g1h_TabSplD3LapJumpd ( nquad, hfunc, dhfunc, ddhfunc,
                              j, lastomcknot, omcknots, fcomc,
                              lastpvknot, pvknots, pv,
                              &trd[(30*ii+5)*nquad],
                              &lapjd[njcurves*(2*s+j)*nquad] );
    }

          /* compute the coefficients by integration */
          /* block A x A */
    Aij = &amat[pkn_Block3FindBlockPos ( hole_k-1, r, S, hole_k-1, hole_k-1 )];
    for ( fn = 0;  fn < nfunc_a;  fn++ ) {
      l1 = pkn_SymMatIndex ( R+fn, R );
      for ( j = 0;  j <= fn;  j++, l1++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac,
                                   &lapja[fn*3*nquad], &lapja[j*3*nquad] );
        Aij[l1] += C*x;
      }
    }

          /* block A x B */
    Aij = &bmat[(nfunc_c+nfunc_d)*nfunc_b];
    for ( fn = l1 = 0;  fn < nfunc_a;  fn++ )
      for ( j = 0;  j < nfunc_b;  j++, l1++ ) {
        supp = _g1h_ExtendSupport ( hole_k, support_b[j] );
        if ( supp & mask ) {
          x = _g1h_Q2SplIntegral0d ( nquad, 3, jac,
                                     &lapja[fn*3*nquad], &lapjb[j*3*nquad] );
          Aij[l1] += C*x;
        }
      }

          /* block A x C */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, i )];
    for ( fn = 0, l1 = (r+nfunc_d)*r;  fn < nfunc_a;  fn++ )
      for ( j = 0;  j < r;  j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapja[fn*3*nquad],
                                   &lapjc[j*njcurves*nquad] );
        if ( i < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+nfunc_d+fn,j)] += C*x;
      }
            /* subblock ii == i-1 */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, ii )];
    for ( fn = 0;  fn < nfunc_a;  fn++ )
      for ( j = 0, l1 = (r+nfunc_d+fn)*r;  j < rr;  j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapja[fn*3*nquad],
                                   &lapjc[(njcurves*r+j)*nquad] );
        if ( ii < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+nfunc_d+fn,j)] += C*x;
      }

          /* block A x D */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, hole_k-1 )];
    for ( fn = 0; fn < nfunc_a; fn++ ) {
      l1 = pkn_SymMatIndex ( r+nfunc_d+fn, r+i*s );
      for ( j = 0; j < s; j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapja[fn*3*nquad],
                                   &lapjd[njcurves*j*nquad] );
        Aij[l1++] += C*x;
      }
    }
    for ( fn = 0; fn < nfunc_a; fn++ ) {
      l1 = pkn_SymMatIndex ( r+nfunc_d+fn, r+jj*s );
      for ( j = 0; j < s; j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapja[fn*3*nquad],
                                   &lapjd[njcurves*(s+j)*nquad] );
        Aij[l1++] += C*x;
      }
    }
    for ( fn = 0; fn < nfunc_a; fn++ ) {
      l1 = pkn_SymMatIndex ( r+nfunc_d+fn, r+ii*s );
      for ( j = 0; j < 2; j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapja[fn*3*nquad],
                                   &lapjd[njcurves*(2*s+j)*nquad] );
        Aij[l1++] += C*x;
      }
    }

          /* block B x C */
    Aij = &bmat[i*r*nfunc_b];
    for ( j = l1 = 0;  j < r;  j++ )  
      for ( fn = 0; fn < nfunc_b; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjc[j*njcurves*nquad],
                                   &lapjb[fn*3*nquad] );
        Aij[l1++] += C*x;
      }
            /* subblock ii == i-1 */
    Aij = &bmat[ii*r*nfunc_b];
    for ( j = 0, l1 = 0;  j < rr;  j++ )
      for ( fn = 0; fn < nfunc_b; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+j)*nquad],
                                   &lapjb[fn*3*nquad] );
        Aij[l1++] += C*x;
      }

          /* block B x D */
    Aij = &bmat[nfunc_b*(nfunc_c+i*s)];
    for ( j = l1 = 0;  j < s;  j++ )
      for ( fn = 0;  fn < nfunc_b;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjb[fn*3*nquad],
                                   &lapjd[njcurves*j*nquad] );
        Aij[l1++] += C*x;
      }
    Aij = &bmat[nfunc_b*(nfunc_c+jj*s)];
    for ( j = l1 = 0;  j < s;  j++ )
      for ( fn = 0;  fn < nfunc_b;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjb[fn*3*nquad],
                                   &lapjd[njcurves*(s+j)*nquad] );
        Aij[l1++] += C*x;
      }
    Aij = &bmat[nfunc_b*(nfunc_c+ii*s)];
    for ( j = l1 = 0;  j < 2;  j++ )
      for ( fn = 0;  fn < nfunc_b;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjb[fn*3*nquad],
                                   &lapjd[njcurves*(2*s+j)*nquad] );
        Aij[l1++] += C*x;
      }

          /* block C x C */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, i, i )];
    for ( fn = 0; fn < r; fn++ )
      for ( j = 0;  j <= fn;  j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjc[fn*njcurves*nquad],
                                   &lapjc[j*njcurves*nquad] );
        Aij[pkn_SymMatIndex(fn,j)] += C*x;
      }
    if ( jumpC ) {
      for ( fn = 0; fn < r; fn++ )
        for ( j = 0; j <= fn; j++ ) {
          x = _g1h_Q2SplIntegral0d ( nquad, 2*nk, &jac[3*nquad],
                                     &lapjc[(fn*njcurves+3)*nquad],
                                     &lapjc[(j*njcurves+3)*nquad] );
          Aij[pkn_SymMatIndex(fn,j)] += C*x;
        }
    }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, ii, ii )];
    for ( fn = 0; fn < rr; fn++ )
      for ( j = 0; j <= fn; j++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+fn)*nquad],
                                   &lapjc[(njcurves*r+j)*nquad] );
        Aij[pkn_SymMatIndex(fn,j)] += C*x;
      }
    if ( i > 0 )
      Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, i, ii )];
    else
      Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, ii, i )];
    for ( fn = l1 = 0;  fn < rr;  fn++ )
      for ( j = 0;  j < r;  j += rr ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+fn)*nquad],
                                   &lapjc[j*njcurves*nquad] );
        if ( i > 0 )
          Aij[j*r+fn] += C*x;
        else
          Aij[fn*r+j] += C*x;
      }

          /* block C x D */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, i )];
    for ( j = 0, l1 = (r+i*s)*r;  j < s;  j++ )
      for ( fn = 0;  fn < r;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjc[fn*njcurves*nquad],
                                   &lapjd[njcurves*j*nquad] );
        if ( jumpC && jumpD )
          x += _g1h_Q2SplIntegral0d ( nquad, nk, &jac[(3+nk)*nquad],
                                      &lapjc[(fn*njcurves+3+nk)*nquad],
                                      &lapjd[(j*njcurves+3+nk)*nquad] );
        if ( i < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+i*s+j,fn)] += C*x;
      }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, i )];
    for ( j = 0; j < 2; j++ )
      for ( fn = 0, l1 = (r+ii*s+j)*r;  fn < r;  fn += rr ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[fn*njcurves*nquad],
                                   &lapjd[(j+2*s)*njcurves*nquad] );
        if ( i < hole_k-1 )
          { Aij[l1] += C*x;  l1 += rr; }
        else
          Aij[pkn_SymMatIndex(r+ii*s+j,fn)] += C*x;
      }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, i )];
    for ( j = 0, l1 = (r+jj*s)*r;  j < s;  j++ )
      for ( fn = 0;  fn < r;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjc[fn*njcurves*nquad],
                                   &lapjd[(s+j)*njcurves*nquad] );
        if ( jumpC && jumpD )
          x += _g1h_Q2SplIntegral0d ( nquad, nk, &jac[3*nquad],
                                      &lapjc[(fn*njcurves+3)*nquad],
                                      &lapjd[((s+j)*njcurves+3)*nquad] );
        if ( i < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+jj*s+j,fn)] += C*x;
      }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, ii )];
    for ( j = 0; j < s; j++ )
      for ( fn = 0, l1 = (r+i*s+j)*r;  fn < rr;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+fn)*nquad],
                                   &lapjd[j*njcurves*nquad] );
        if ( ii < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+i*s+j,fn)] += C*x;
      }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, ii )];
    for ( j = 0;  j < 2;  j++ )
      for ( fn = 0, l1 = (r+jj*s+j)*r;  fn < rr;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+fn)*nquad],
                                   &lapjd[njcurves*(s+j)*nquad] );
        if ( ii < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+jj*s+j,fn)] += C*x;
      }
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, ii )];
    for ( j = 0;  j < 2;  j++ )
      for ( fn = 0, l1 = (r+ii*s+j)*r;  fn < rr;  fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjc[(njcurves*r+fn)*nquad],
                                   &lapjd[njcurves*(2*s+j)*nquad] );
        if ( ii < hole_k-1 )
          Aij[l1++] += C*x;
        else
          Aij[pkn_SymMatIndex(r+ii*s+j,fn)] += C*x;
      }

          /* block D x D */
    Aij = &amat[pkn_Block3FindBlockPos( hole_k-1, r, S, hole_k-1, hole_k-1 )];
    for ( j = 0; j < 2; j++ )
      for ( fn = 0; fn <= j; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjd[(j+2*s)*njcurves*nquad],
                                   &lapjd[(fn+2*s)*njcurves*nquad] );
        Aij[pkn_SymMatIndex(r+ii*s+j,r+ii*s+fn)] += C*x;
      }
    for ( j = 0; j < s; j++ )
      for ( fn = 0; fn <= j; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjd[j*njcurves*nquad],
                                   &lapjd[fn*njcurves*nquad] );
        if ( jumpD )
          x += _g1h_Q2SplIntegral0d ( nquad, nk, &jac[(3+nk)*nquad],
                                      &lapjd[(fn*njcurves+3+nk)*nquad],
                                      &lapjd[(j*njcurves+3+nk)*nquad] );
        Aij[pkn_SymMatIndex(r+i*s+j,r+i*s+fn)] += C*x;
      }
    for ( j = 0; j < s; j++ )
      for ( fn = 0; fn <= j; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjd[(j+s)*njcurves*nquad],
                                   &lapjd[(fn+s)*njcurves*nquad] );
        if ( jumpD )
          x += _g1h_Q2SplIntegral0d ( nquad, nk, &jac[3*nquad],
                                      &lapjd[((fn+s)*njcurves+3)*nquad],
                                      &lapjd[((j+s)*njcurves+3)*nquad] );
        Aij[pkn_SymMatIndex(r+jj*s+j,r+jj*s+fn)] += C*x;
      }
    for ( j = 0; j < 2; j++ )
      for ( fn = 0; fn < s; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjd[(j+2*s)*njcurves*nquad],
                                   &lapjd[fn*njcurves*nquad] );
        Aij[pkn_SymMatIndex(r+ii*s+j,r+i*s+fn)] += C*x;
      }
    for ( j = 0; j < 2; j++ )
      for ( fn = 0; fn < s; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 1, jac, &lapjd[(j+2*s)*njcurves*nquad],
                                   &lapjd[(fn+s)*njcurves*nquad] );
        Aij[pkn_SymMatIndex(r+ii*s+j,r+jj*s+fn)] += C*x;
      }
    for ( j = 0; j < s; j++ )
      for ( fn = 0; fn < s; fn++ ) {
        x = _g1h_Q2SplIntegral0d ( nquad, 3, jac, &lapjd[j*njcurves*nquad],
                                   &lapjd[(fn+s)*njcurves*nquad] );
        Aij[pkn_SymMatIndex(r+i*s+j,r+jj*s+fn)] += C*x;
      }

  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2SplComputeFormMatrixd*/

boolean g1h_Q2SplDecomposeMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd   *privateG1;
  G1HoleSPrivateRecd  *privateS;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, bldim1, bldim2;
  int    size;
  double *lmat;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  privateS = domain->SprivateG1;
  if ( !privateS->Q2SLMat ) {
    if ( !privateS->Q2SAMat )
      if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
        goto failure;
    hole_k = domain->hole_k;
    nfunc_a = privateG1->nfunc_a;
    nfunc_c = privateS->nsfunc_c;
    nfunc_d = privateS->nsfunc_d;
    bldim1 = nfunc_c/hole_k;
    bldim2 = bldim1+nfunc_d+nfunc_a;

    size = pkn_Block3ArraySize ( hole_k-1, bldim1, bldim2 );
    lmat = privateS->Q2SLMat = malloc ( size*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateS->Q2SAMat, size*sizeof(double) );
    if ( !pkn_Block3CholeskyDecompMd ( hole_k-1, bldim1, bldim2, lmat ) ) {
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2SplDecomposeMatrixd*/

boolean g1h_Q2SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,              
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,              
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecd   *privateG1;
  G1HoleSPrivateRecd  *privateS;
  int     hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, bldim1, bldim2;
  double  *bmat, *lmat, *x, *fc00;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  privateS  = domain->SprivateG1;
  if ( !privateS )
    goto failure;
  if ( !privateS->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateS->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = privateS->nsfunc_c;
  nfunc_d = privateS->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  lmat = privateS->Q2SLMat;
  bmat = privateS->Q2SBMat;
  bldim1 = nfunc_c/hole_k;
  bldim2 = bldim1+nfunc_d+nfunc_a;

  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
    goto failure;

  pkn_Block3LowerTrMSolved ( hole_k-1, bldim1, bldim2,
                             lmat, spdimen, spdimen, x );
  pkn_Block3UpperTrMSolved ( hole_k-1, bldim1, bldim2,
                             lmat, spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2SplFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr, 
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *cmat, *rcmat, *lmat, *bmat, *aa;
  double s, t;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->SCmat || sprivate->splnconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  lmat = sprivate->Q2SLMat;
  bmat = sprivate->Q2SBMat;
  rcmat = sprivate->Q2SRCmat;
  if ( !rcmat ) {
    rcmat = sprivate->Q2SRCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
    if ( !rcmat )
      goto failure;
    cmat = pkv_GetScratchMemd ( (nbf+2)*nconstr );
    if ( !cmat ) {
      domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
      goto failure;
    }
    aa = &cmat[nbf*nconstr];
    pkv_TransposeMatrixd ( nconstr, nbf, nbf, sprivate->SCmat, nconstr, cmat );
    pkn_Block3LowerTrMSolved ( hole_k-1, nfunc_c/hole_k,
                               nfunc_c/hole_k+nfunc_d+nfunc_a, lmat,
                               nconstr, nconstr, cmat );
    if ( !pkn_QRDecomposeMatrixd ( nbf, nconstr, cmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < nconstr; i++ )
      for ( j = i; j < nconstr; j++ )
        rcmat[pkn_SymMatIndex(i,j)] = cmat[i*nconstr+j];
    for ( i = 0; i < nconstr; i++ ) {
      for ( j = 0, s = 0.0;  j < i;  j++ ) {
        t = rcmat[pkn_SymMatIndex(i,j)];
        s += t*t;
      }
      t = rcmat[pkn_SymMatIndex(i,i)];
      if ( t*t < 1.0e-10*s ) {
        domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
        goto failure;
      }
    }
    pkv_FreeScratchMemd ( (nbf+2)*nconstr );
  }
  cmat = sprivate->SCmat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
    goto failure;

  pkn_Block3LowerTrMSolved ( hole_k-1, nfunc_c/hole_k,
                             nfunc_c/hole_k+nfunc_d+nfunc_a, lmat,
                             spdimen, spdimen, x );
  pkn_Block3UpperTrMSolved ( hole_k-1, nfunc_c/hole_k,
                             nfunc_c/hole_k+nfunc_d+nfunc_a, lmat,
                             spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixd ( nconstr, nbf, nbf, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixd ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemd ( spdimen*nbf );
  if ( !y ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixd ( nconstr, nbf, nbf, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_Block3LowerTrMSolved ( hole_k-1, nfunc_c/hole_k,
                             nfunc_c/hole_k+nfunc_d+nfunc_a, lmat,
                             spdimen, spdimen, y );
  pkn_Block3UpperTrMSolved ( hole_k-1, nfunc_c/hole_k,
                             nfunc_c/hole_k+nfunc_d+nfunc_a, lmat,
                             spdimen, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nbf, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemd ( spdimen*nbf );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );


  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2SplFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp, *sp1;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *acmat, *rcmat, *lmat, *bmat, *aa, s, t;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nbf;
  int    d1, d2, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->ASCmat || sprivate->splnaconstr != naconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  lmat = sprivate->Q2SLMat;
  bmat = sprivate->Q2SBMat;
  rcmat = sprivate->Q2SARCmat;
  d1 = nfunc_c/hole_k;
  d2 = d1+nfunc_d+nfunc_a;
  if ( !rcmat ) {
    sp1 = pkv_GetScratchMemTop ();
    rcmat = sprivate->Q2SARCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
    if ( !rcmat )
      goto failure;

    if ( !(acmat = pkv_GetScratchMemd ( (nbf*spdimen+2)*naconstr )) ) {
      domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
      goto failure;
    }
    aa = &acmat[nbf*spdimen*naconstr];
    pkv_TransposeMatrixd ( naconstr, nbf*spdimen, nbf*spdimen,
                           sprivate->ASCmat, naconstr, acmat );
    for ( i = 0; i < spdimen; i++ )
      pkn_Block3LowerTrMSolved ( hole_k-1, d1, d2, lmat,
                                 naconstr, naconstr, &acmat[i*nbf*naconstr] );
    if ( !pkn_QRDecomposeMatrixd ( nbf*spdimen, naconstr, acmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < naconstr; i++ )
      for ( j = i; j < naconstr; j++ )
        rcmat[pkn_SymMatIndex(i,j)] = acmat[i*naconstr+j];
    for ( i = 0; i < naconstr; i++ ) {
      for ( j = 0, s = 0.0;  j < i;  j++ ) {
        t = rcmat[pkn_SymMatIndex(i,j)];
        s += t*t;
      }
      t = rcmat[pkn_SymMatIndex(i,i)];
      if ( t*t < 1.0e-10*s ) {
        domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
        goto failure;
      }
    }
    pkv_SetScratchMemTop ( sp1 );
  }
  acmat = sprivate->ASCmat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
    goto failure;

  pkn_Block3LowerTrMSolved ( hole_k-1, d1, d2, lmat, spdimen, spdimen, x );
  pkn_Block3UpperTrMSolved ( hole_k-1, d1, d2, lmat, spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixd ( nbf, spdimen, spdimen, x, nbf, y );
  pkn_MultMatrixd ( naconstr, spdimen*nbf, spdimen*nbf, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixd ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( naconstr, rcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolved ( naconstr, rcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixd ( naconstr, nbf*spdimen, nbf*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixd ( spdimen, nbf, nbf, y, spdimen, b );
  pkn_Block3LowerTrMSolved ( hole_k-1, d1, d2, lmat, spdimen, spdimen, b );
  pkn_Block3UpperTrMSolved ( hole_k-1, d1, d2, lmat, spdimen, spdimen, b );
  pkn_SubtractMatrixd ( 1, spdimen*nbf, 0, x, 0, b, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2SplFillHoleAltConstrd*/

