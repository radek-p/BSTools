
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"

#define _DEBUG

#ifdef DEBUG
#define QUAD_FACTOR 2
#else
#define QUAD_FACTOR 10
#endif

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2h_TabBSFuncDer3f ( int deg, int lastknot, const float *knots,
                              int i0, int i1,
                              int n, const float *tkn, int *fkn, int *lkn,
                              float *b, float *bt, float *btt, float *bttt )
{
  /* evaluate B-spline basis functions at a number of points  deg - degree,  */
  /* lastknot - number of the last knot, knots - knots, i0, i1 - numbers of  */
  /* the first and the last function to evaluate, n - number of (domain)     */
  /* points, tkn - the points (ascending), fkn, lkn - arrays in which for    */
  /* each function the numbers of the first and the last point, at which the */
  /* function is nonzero, are stored, b, bt, btt, bttt - arrays in which     */
  /* the values of the functions and their derivatives are to be stored      */
  
  void  *sp;
  int   i, j, l, fk, lk;
  float *p;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMemf ( lastknot-deg );
  if ( !p )
    goto failure;

  memset ( b, 0, (i1-i0+1)*n*sizeof(float) );
  memset ( bt, 0, (i1-i0+1)*n*sizeof(float) );
  memset ( btt, 0, (i1-i0+1)*n*sizeof(float) );
  memset ( bttt, 0, (i1-i0+1)*n*sizeof(float) );

        /* the i-th function is nonzero for the arguments between knot[i] */
        /* and knot[i+deg+1] */
  for ( fk = 0; tkn[fk] < knots[i0]; fk ++ ) ;
  for ( lk = fk+1; lk+1 < n && tkn[lk+1] <= knots[i0+deg+1]; lk ++ ) ;
  memset ( p, 0, (lastknot-deg)*sizeof(float) );
  for ( i = i0, l = 0;  i <= i1;  i ++, l ++ ) {
    p[i] = 1.0;
    while ( tkn[fk] < knots[i] ) fk ++;
    while ( lk+1 < n && tkn[lk+1] <= knots[i+deg+1] ) lk ++;
    fkn[l] = fk;
    lkn[l] = lk;
    for ( j = fk; j <= lk; j++ )
      mbs_deBoorDer3C1f ( deg, lastknot, knots, p, tkn[j],
                          &b[l*n+j], &bt[l*n+j], &btt[l*n+j], &bttt[l*n+j] );
    p[i] = 0.0;
 }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabBSFuncDer3f*/

void _g2h_TensDer3f ( float p, float pu, float puu, float puuu,
                      float q, float qv, float qvv, float qvvv,
                      float *pq )
{
  pq[0] = pu*q;    pq[1] = p*qv;
  pq[2] = puu*q;   pq[3] = pu*qv;   pq[4] = p*qvv;
  pq[5] = puuu*q;  pq[6] = puu*qv;  pq[7] = pu*qvv;  pq[8] = p*qvvv;
} /*_g2h_TensDer3f*/

void _g2h_TabSplCLaplacianGradf ( int nkn,
            int fknp, int lknp, const float *tp, const float *tpu,
            const float *tpuu, const float *tpuuu,
            int fknq, int lknq, const float *tq, const float *tqv,
            const float *tqvv, const float *tqvvv,
            const float *trd, vector2f *lgr )
{
  float pq[9];
  int   i, j, k, kk;

  memset ( lgr, 0, nkn*nkn*sizeof(vector2f) );
  for ( i = fknp; i <= lknp; i++ ) {
    for ( j = fknq, k = i*nkn+fknq, kk = 18*k;
          j <= lknq;
          j++, k++, kk += 18 ) {
      _g2h_TensDer3f ( tp[i], tpu[i], tpuu[i], tpuuu[i],
                       tq[j], tqv[j], tqvv[j], tqvvv[j], pq );
      lgr[k].x = (float)pkn_ScalarProductf ( 9, &trd[kk], pq );
      lgr[k].y = (float)pkn_ScalarProductf ( 9, &trd[kk+9], pq );
    }
  }
} /*_g2h_TabSplCLaplacianGradf*/

boolean _g2h_FuncDSuppf ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 )
{
  /* this function returns true if the function fn (in the block D) is  */
  /* nonzero in Omega_i; then the values assigned to *i0, *i1, *j0, *j1 */
  /* indicate the ranges of indexes of the quadrature knots, at which   */
  /* the function takes nonzero values                                  */
  int j, k, l, cn;

  cn = fn / (3*nk*m1);
  if ( cn == i ) {
    *j0 = 0;
    *j1 = l = (nk+1)*QUAD_FACTOR-1;
    fn %= 3*nk*m1;
    j = fn / 3;           /* which knot? */
    *nzc = k = fn % 3;    /* which order of cross derivative is nonzero? */
    if ( j <= 2 ) *i0 = 0;
    else *i0 = ((j-3)/m1+1)*QUAD_FACTOR;
    k = ((5-k+j)/m1+1)*QUAD_FACTOR-1;
    *i1 = min ( k, l );
  }
  else if ( (i+1) % hole_k == cn ) {
    *i0 = 0;
    *i1 = l = (nk+1)*QUAD_FACTOR-1;
    fn %= 3*nk*m1;
    j = fn / 3;           /* which knot? */
    *nzc = k = fn % 3;    /* which order of cross derivative is nonzero? */
    if ( j <= 2 ) *j0 = 0;
    else *j0 = ((j-3)/m1+1)*QUAD_FACTOR;
    k = ((5-k+j)/m1+1)*QUAD_FACTOR-1;
    *j1 = min ( k, l );
  }
  else
    return false;

  return true;
} /*_g2h_FuncDSuppf*/

static boolean _g2h_TabSplD1LaplacianGraphf ( int nkn, const float *tkn,
                   const float *hfunc, const float *dhfunc,
                   const float *ddhfunc, const float *dddhfunc,
                   int nzc, int i0, int i1,
                   int lastomcknot, const float *omcknots, const float *fcomc,
                   int lastpvknot, const float *pvknots, const float *pv,
                   int lastpvvknot, const float *pvvknots, const float *pvv,
                   const float *trd,
                   vector2f *lgr )
{
  void  *sp;
  float *fp, bp[9];
  int   i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemf ( 12*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 12*nkn*sizeof(float) );
  ii = 12*i0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer3f ( 1, G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
           i1-i0+1, &tkn[i0], 12, &fp[ii], &fp[ii+1], &fp[ii+2], &fp[ii+3] );
  if ( nzc <= 1 )
    mbs_TabBSCurveDer3f ( 1, G2_CROSS01DEG, lastpvknot, pvknots, pv,
           i1-i0+1, &tkn[i0], 12, &fp[ii+4], &fp[ii+5], &fp[ii+6], &fp[ii+7] );
  mbs_TabBSCurveDer3f ( 1, G2_CROSS02DEG, lastpvvknot, pvvknots, pvv,
         i1-i0+1, &tkn[i0], 12, &fp[ii+8], &fp[ii+9], &fp[ii+10], &fp[ii+11] );

  memset ( lgr, 0, nkn*nkn*sizeof(vector2f) );
  for ( i = i0;  i <= i1;  i++, ii += 12 ) {
    for ( j = jj = 0, k = nkn*i, kk = 18*k;
          j < nkn;
          j++, jj += 6, k++, kk += 18 ) {
      switch ( nzc ) {
    default:
        bp[0] = fp[ii+1]*hfunc[jj]   + fp[ii+5]*hfunc[jj+2]    + fp[ii+9]*hfunc[jj+4];
        bp[2] = fp[ii+2]*hfunc[jj]   + fp[ii+6]*hfunc[jj+2]    + fp[ii+10]*hfunc[jj+4];
        bp[5] = fp[ii+3]*hfunc[jj]   + fp[ii+7]*hfunc[jj+2]    + fp[ii+11]*hfunc[jj+4];
        bp[1] = fp[ii]*dhfunc[jj]    + fp[ii+4]*dhfunc[jj+2]   + fp[ii+8]*dhfunc[jj+4];
        bp[3] = fp[ii+1]*dhfunc[jj]  + fp[ii+5]*dhfunc[jj+2]   + fp[ii+9]*dhfunc[jj+4];
        bp[6] = fp[ii+2]*dhfunc[jj]  + fp[ii+6]*dhfunc[jj+2]   + fp[ii+10]*dhfunc[jj+4];
        bp[4] = fp[ii]*ddhfunc[jj]   + fp[ii+4]*ddhfunc[jj+2]  + fp[ii+8]*ddhfunc[jj+4];
        bp[7] = fp[ii+1]*ddhfunc[jj] + fp[ii+5]*ddhfunc[jj+2]  + fp[ii+9]*ddhfunc[jj+4];
        bp[8] = fp[ii]*dddhfunc[jj]  + fp[ii+4]*dddhfunc[jj+2] + fp[ii+8]*dddhfunc[jj+4];
        break;
    case 1:
        bp[0] = fp[ii+5]*hfunc[jj+2]    + fp[ii+9]*hfunc[jj+4];
        bp[2] = fp[ii+6]*hfunc[jj+2]    + fp[ii+10]*hfunc[jj+4];
        bp[5] = fp[ii+7]*hfunc[jj+2]    + fp[ii+11]*hfunc[jj+4];
        bp[1] = fp[ii+4]*dhfunc[jj+2]   + fp[ii+8]*dhfunc[jj+4];
        bp[3] = fp[ii+5]*dhfunc[jj+2]   + fp[ii+9]*dhfunc[jj+4];
        bp[6] = fp[ii+6]*dhfunc[jj+2]   + fp[ii+10]*dhfunc[jj+4];
        bp[4] = fp[ii+4]*ddhfunc[jj+2]  + fp[ii+8]*ddhfunc[jj+4];
        bp[7] = fp[ii+5]*ddhfunc[jj+2]  + fp[ii+9]*ddhfunc[jj+4];
        bp[8] = fp[ii+4]*dddhfunc[jj+2] + fp[ii+8]*dddhfunc[jj+4];
        break;
    case 2:
        bp[0] = fp[ii+9]*hfunc[jj+4];
        bp[2] = fp[ii+10]*hfunc[jj+4];
        bp[5] = fp[ii+11]*hfunc[jj+4];
        bp[1] = fp[ii+8]*dhfunc[jj+4];
        bp[3] = fp[ii+9]*dhfunc[jj+4];
        bp[6] = fp[ii+10]*dhfunc[jj+4];
        bp[4] = fp[ii+8]*ddhfunc[jj+4];
        bp[7] = fp[ii+9]*ddhfunc[jj+4];
        bp[8] = fp[ii+8]*dddhfunc[jj+4];
        break;
      }
      lgr[k].x = (float)pkn_ScalarProductf ( 9, &trd[kk], bp );
      lgr[k].y = (float)pkn_ScalarProductf ( 9, &trd[kk+9], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabSplD1LaplacianGraphf*/

static boolean _g2h_TabSplD2LaplacianGraphf ( int nkn, const float *tkn,
                   const float *hfunc, const float *dhfunc,
                   const float *ddhfunc, const float *dddhfunc,
                   int nzc, int j0, int j1,
                   int lastomcknot, const float *omcknots, const float *fcomc,
                   int lastpuknot, const float *puknots, const float *pu,
                   int lastpuuknot, const float *puuknots, const float *puu,
                   const float *trd,
                   vector2f *lgr )
{
  void  *sp;
  float *fp, bp[9];
  int   i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemf ( 12*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 12*nkn*sizeof(float) );
  jj = 12*j0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer3f ( 1, G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
           j1-j0+1, &tkn[j0], 12, &fp[jj], &fp[jj+1], &fp[jj+2], &fp[jj+3] );
  if ( nzc <= 1 )
    mbs_TabBSCurveDer3f ( 1, G2_CROSS01DEG, lastpuknot, puknots, pu,
           j1-j0+1, &tkn[j0], 12, &fp[jj+4], &fp[jj+5], &fp[jj+6], &fp[jj+7] );
  mbs_TabBSCurveDer3f ( 1, G2_CROSS02DEG, lastpuuknot, puuknots, puu,
         j1-j0+1, &tkn[j0], 12, &fp[jj+8], &fp[jj+9], &fp[jj+10], &fp[jj+11] );

  memset ( lgr, 0, nkn*nkn*sizeof(vector2f) );
  for ( i = ii = 0;  i < nkn;  i++, ii += 6 ) {
    for ( j = j0, jj = 12*j0, k = nkn*i+j0, kk = 18*k;
          j <= j1;
          j++, jj += 12, k++, kk += 18 ) {
      switch ( nzc ) {
    default:
        bp[1] = fp[jj+1]*hfunc[ii]   + fp[jj+5]*hfunc[ii+2]    + fp[jj+9]*hfunc[ii+4];
        bp[4] = fp[jj+2]*hfunc[ii]   + fp[jj+6]*hfunc[ii+2]    + fp[jj+10]*hfunc[ii+4];
        bp[8] = fp[jj+3]*hfunc[ii]   + fp[jj+7]*hfunc[ii+2]    + fp[jj+11]*hfunc[ii+4];
        bp[0] = fp[jj]*dhfunc[ii]    + fp[jj+4]*dhfunc[ii+2]   + fp[jj+8]*dhfunc[ii+4];
        bp[3] = fp[jj+1]*dhfunc[ii]  + fp[jj+5]*dhfunc[ii+2]   + fp[jj+9]*dhfunc[ii+4];
        bp[7] = fp[jj+2]*dhfunc[ii]  + fp[jj+6]*dhfunc[ii+2]   + fp[jj+10]*dhfunc[ii+4];
        bp[2] = fp[jj]*ddhfunc[ii]   + fp[jj+4]*ddhfunc[ii+2]  + fp[jj+8]*ddhfunc[ii+4];
        bp[6] = fp[jj+1]*ddhfunc[ii] + fp[jj+5]*ddhfunc[ii+2]  + fp[jj+9]*ddhfunc[ii+4];
        bp[5] = fp[jj]*dddhfunc[ii]  + fp[jj+4]*dddhfunc[ii+2] + fp[jj+8]*dddhfunc[ii+4];
        break;
    case 1:
        bp[1] = fp[jj+5]*hfunc[ii+2]    + fp[jj+9]*hfunc[ii+4];
        bp[4] = fp[jj+6]*hfunc[ii+2]    + fp[jj+10]*hfunc[ii+4];
        bp[8] = fp[jj+7]*hfunc[ii+2]    + fp[jj+11]*hfunc[ii+4];
        bp[0] = fp[jj+4]*dhfunc[ii+2]   + fp[jj+8]*dhfunc[ii+4];
        bp[3] = fp[jj+5]*dhfunc[ii+2]   + fp[jj+9]*dhfunc[ii+4];
        bp[7] = fp[jj+6]*dhfunc[ii+2]   + fp[jj+10]*dhfunc[ii+4];
        bp[2] = fp[jj+4]*ddhfunc[ii+2]  + fp[jj+8]*ddhfunc[ii+4];
        bp[6] = fp[jj+5]*ddhfunc[ii+2]  + fp[jj+9]*ddhfunc[ii+4];
        bp[5] = fp[jj+4]*dddhfunc[ii+2] + fp[jj+8]*dddhfunc[ii+4];
        break;
    case 2:
        bp[1] = fp[jj+9]*hfunc[ii+4];
        bp[4] = fp[jj+10]*hfunc[ii+4];
        bp[8] = fp[jj+11]*hfunc[ii+4];
        bp[0] = fp[jj+8]*dhfunc[ii+4];
        bp[3] = fp[jj+9]*dhfunc[ii+4];
        bp[7] = fp[jj+10]*dhfunc[ii+4];
        bp[2] = fp[jj+8]*ddhfunc[ii+4];
        bp[6] = fp[jj+9]*ddhfunc[ii+4];
        bp[5] = fp[jj+8]*dddhfunc[ii+4];
        break;
      }
      lgr[k].x = (float)pkn_ScalarProductf ( 9, &trd[kk], bp );
      lgr[k].y = (float)pkn_ScalarProductf ( 9, &trd[kk+9], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabSplD2LaplacianGraphf*/

#ifdef DEBUG
static boolean dumpspl = false;

float _g2h_SplIntegralf ( int nkn, float *jac,
                          int i00, int i01, int j00, int j01, vector2f *func0,
                          int i10, int i11, int j10, int j11, vector2f *func1 )
{
  FILE   *f;
  int    i, j, k;
  double s, ss;

  if ( dumpspl )
    f = fopen ( "intspl.txt", "w+" );
        /* two-level summing is supposed to decrease the effect */
        /* of rounding errors */
  i00 = max ( i00, i10 );  i01 = min ( i01, i11 );
  j00 = max ( j00, j10 );  j01 = min ( j01, j11 );
  for ( i = i00, ss = 0.0; i <= i01; i++ ) {
    for ( j = j00, k = nkn*i+j00, s = 0.0;  j <= j01;  j++, k++ ) {
      if ( dumpspl )
        fprintf ( f, "%d, %e, %e, %e, %e %e\n",
                  k, func0[k].x, func0[k].y, func1[k].x, func1[k].y, jac[k] );
      s += (func0[k].x*func1[k].x + func0[k].y*func1[k].y)*jac[k];
    }
    ss += s;
  }
  if ( dumpspl ) {
    fprintf ( f, "\n" );
    fclose ( f );
  }
  return ss;
} /*_g2h_SplIntegralf*/
#else
float _g2h_SplIntegralf ( int nkn, float *jac,
                          int i00, int i01, int j00, int j01, vector2f *func0,
                          int i10, int i11, int j10, int j11, vector2f *func1 )
{
  int    i, j, k;
  double s, ss;

        /* two-level summing is supposed to decrease the effect */
        /* of rounding errors */
  i00 = max ( i00, i10 );  i01 = min ( i01, i11 );
  j00 = max ( j00, j10 );  j01 = min ( j01, j11 );
  for ( i = i00, ss = 0.0; i <= i01; i++ ) {
    for ( j = j00, k = nkn*i+j00, s = 0.0;  j <= j01;  j++, k++ )
      s += (func0[k].x*func1[k].x + func0[k].y*func1[k].y)*jac[k];
    ss += s;
  }
  return (float)ss;
} /*_g2h_SplIntegralf*/
#endif
/* ///////////////////////////////////////////////////////////////////////// */
#ifdef DEBUG
void DumpJacobian ( int i, int nkn, const float *jac )
{
  FILE *f;
  int  j;

  if ( i == 0 )
    f = fopen ( "jac.txt", "w+" );
  else
    f = fopen ( "jac.txt", "a+" );

  fprintf ( f, "i:%d\n", i );
  for ( j = 0; j < nkn*nkn; j++ )
    fprintf ( f, "%4d: %f\n", j, jac[j] );
  fprintf ( f, "\n" );
  fclose ( f );
} /*DumpJacobian*/

void DumpLapGrad ( int i, int f0, int nfunc, int nkn, const vector2f *lgr )
{
  FILE *f;
  int  j, k, l;

  if ( i == 0 )
    f = fopen ( "lgr.txt", "w+" );
  else
    f = fopen ( "lgr.txt", "a+" );

  fprintf ( f, "i:%d\n", i );
  for ( j = f0, l = f0*nkn*nkn;  j < f0+nfunc;  j++ ) {
    fprintf ( f, "fn=%d\n", j );
    for ( k = 0;  k < nkn*nkn;  k++, l++ )
      fprintf ( f, "%4d: %e,%e\n", k, lgr[l].x, lgr[l].y );
    fprintf ( f, "\n" );
  }
  fprintf ( f, "\n" );
  fclose ( f );
} /*DumpLapGrad*/

void DumpTBS ( int nf, int nkn, const float *tp, const float *tpu,
               const float *tpuu, const float *tpuuu )
{
  FILE   *f;
  int    i, j, k, l;
  float  pq[9];

  f = fopen ( "tbs.txt", "w+" );
  for ( i = 0; i < nf; i++ )
    for ( j = 0; j < nf; j++ ) {
      fprintf ( f, "fn:%d\n", nf*i+j );
      for ( k = 0; k < nkn; k++ )
        for ( l = 0; l < nkn; l++ ) {
          _g2h_TensDer3f ( tp[i*nkn+k], tpu[i*nkn+k], tpuu[i*nkn+k], tpuuu[i*nkn+k],
                           tp[j*nkn+l], tpu[j*nkn+l], tpuu[j*nkn+l], tpuuu[j*nkn+l],
                           pq );
          fprintf ( f,
            " pu:%e, pv:%e, puu:%e, puv:%e, puu:%e, puuu:%e, puuv:%e, puvv:%e, pvvv:%e\n",
            pq[0], pq[1], pq[2], pq[3], pq[4], pq[5], pq[6], pq[7], pq[8] );
        }
      fprintf ( f, "\n" );
    }
  fclose ( f );
} /*DumpTBS*/
#endif
/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_ComputeSplFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf   *privateG;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nab, nabc;
  int      lastomcknot, lastpvknot, lastpvvknot;
  float    *omcknots, *pvknots, *pvvknots;
  int      rr, r, s, asize, bsize;
  float    *SAMat, *SBMat, *auxmat1, *auxmat2;
  unsigned short *support_b, mask;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float    *fcomc, *pv, *pvv, *pu, *puu;
  float    *trd, *jac, *tbs, *tbst, *tbstt, *tbsttt;
  int      *fkn, *lkn;
  vector2f *lgr;
  int      nk, m1, m2, m, nquad, nquadsq;
  int      i, ii, j, k, l1, l2, twok;
  int      i00, i01, j00, j01, i10, i11, j10, j11, nzc;
  double   qc;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  twok = 2*hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  support_b = privateG->support_b;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;

  omcknots = sprivate->omcknots;  
  lastomcknot = sprivate->lastomcknot;
  pvknots = sprivate->pvknots;
  lastpvknot = sprivate->lastpvknot;
  pvvknots = sprivate->pvvknots;
  lastpvvknot = sprivate->lastpvvknot;

        /* compute the size of the arrays and allocate memory */
  rr = G2H_FINALDEG-5+nk*m2;
  r = rr*rr;    /* == nfunc_c/hole_k */
  s = 3*nk*m1;  /* == nfunc_d/hole_k */
  nab = nfunc_a+nfunc_b;
  nabc = nab+r;
  asize = pkn_Block2ArraySize ( hole_k, r, s, nfunc_a );
  bsize = nfunc_b*(nfunc_a+nfunc_c+nfunc_d);
  SAMat = sprivate->SAMat = malloc ( asize*sizeof(float) );
  SBMat = sprivate->SBMat = malloc ( bsize*sizeof(float) );
  if ( !SAMat || !SBMat )
    goto failure;

        /* fix the number of quadrature knots, perhaps later to be */
        /* introduced by the application via OptionProc */
  nquad = QUAD_FACTOR*(nk+1);
  nquadsq = nquad*nquad;
  m = rr*nquad;

        /* allocate the scratch memory */
  tkn = pkv_GetScratchMemf ( 25*nquad );
  jac = pkv_GetScratchMemf ( nquadsq );
  trd = pkv_GetScratchMemf ( 18*nquadsq );
  fkn = pkv_GetScratchMem ( 2*rr*sizeof(int) );
  lgr = (vector2f*)pkv_GetScratchMem ( (nabc+2*s)*nquadsq*sizeof(vector2f) );
  tbs = pkv_GetScratchMemf ( 4*m );
  fcomc = pkv_GetScratchMemf ( lastomcknot+2*(lastpvknot+lastpvvknot)-
                               (G2_CROSS00DEG+2*(G2_CROSS01DEG+G2_CROSS02DEG)) );
  if ( !tkn || !jac || !trd || !fkn || !lgr || !tbs || !fcomc )
    goto failure;
  hfunc = &tkn[nquad];         dhfunc = &hfunc[6*nquad];
  ddhfunc = &dhfunc[6*nquad];  dddhfunc = &ddhfunc[6*nquad];
  tbst = &tbs[m];              tbstt = &tbst[m];
  tbsttt = &tbstt[m];
  pv = &fcomc[lastomcknot-G2_CROSS00DEG];
  pvv = &pv[lastpvknot-G2_CROSS01DEG];
  pu = &pvv[lastpvvknot-G2_CROSS02DEG];
  puu = &pu[lastpvknot-G2_CROSS01DEG];
  lkn = &fkn[rr];

        /* prepare the evaluation of basis functions */
  _gh_PrepareTabKnotsf ( nquad, privateG2->opt_quad, tkn );
  mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, nquad, tkn,
                             hfunc, dhfunc, ddhfunc, dddhfunc );

  if ( !_g2h_TabBSFuncDer3f ( G2H_FINALDEG, sprivate->lastcknot, sprivate->cknots,
                              3, G2H_FINALDEG+nk*m2-3, nquad, tkn, fkn, lkn,
                              tbs, tbst, tbstt, tbsttt ) )
    goto failure;

        /* compute and sum the integrals in subsequent areas omega_i */
  memset ( SAMat, 0, asize*sizeof(float) );
  memset ( SBMat, 0, bsize*sizeof(float) );
  for ( i = 0, ii = 1, mask = 0x0001;
        i < hole_k;
        i++, ii = (ii+1) % hole_k, mask = (unsigned short)(mask << 1) ) {

          /* evaluate the Jacobian and the matrix of derivatives of the */
          /* domain patch at the quadrature knots */
    _g2h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    _g2h_TabDiPatchJac3f ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                           c00, c01, c02, c10, c11, c12,
                           d00, d01, d02, d10, d11, d12,
                           jac, trd );

#ifdef DEBUG
DumpJacobian ( i, nquad, jac );
#endif
          /* evaluate the Laplacian gradient of the A block basis functions */
          /* at the quadrature knots */
    for ( j = 0; j < nfunc_a; j++ ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      _g2h_TabLaplacianGrad0f ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                fc00, fc01, fc02, fd00, fd01, fd02, trd,
                                &lgr[j*nquadsq] );
    }

          /* evaluate the Laplacian gradient of the B block basis functions */
          /* at the quadrature knots */
    for ( j = 0; j < nfunc_b; j++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                  &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                  &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
        _g2h_TabLaplacianGradf ( nquad, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                 fc00, fc01, fc02, fc10, fc11, fc12,
                                 fd00, fd01, fd02, fd10, fd11, fd12,
                                 trd, &lgr[(nfunc_a+j)*nquadsq] );
      }

          /* evaluate the Laplacian gradient of the C block basis functions */
          /* at the quadrature knots */
    for ( j = l1 = 0;  j < rr;  j++ )
      for ( k = 0;  k < rr;  k++, l1++ )
        _g2h_TabSplCLaplacianGradf ( nquad,
                  fkn[j], lkn[j], &tbs[j*nquad], &tbst[j*nquad],
                  &tbstt[j*nquad], &tbsttt[j*nquad],
                  fkn[k], lkn[k], &tbs[k*nquad], &tbst[k*nquad],
                  &tbstt[k*nquad], &tbsttt[k*nquad],
                  trd, &lgr[(nab+l1)*nquadsq] );

          /* evaluate the Laplacian gradient of the D block basis functions */
          /* at the quadrature knots; do it only for the functions nonzero  */
          /* in Omega_i; j is the function number within the subblock, k is */
          /* the complete function number */
    for ( j = 0, k = nab+nfunc_c+i*s;
          j < s;
          j++, k++ ) {
      _g2h_GetSplDBasisCrossDerf ( domain, k, i, fcomc, pv, pvv, pu, puu );
      _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g2h_TabSplD1LaplacianGraphf ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                   dddhfunc, nzc, i00, i01, lastomcknot, omcknots, fcomc,
                   lastpvknot, pvknots, pv, lastpvvknot, pvvknots, pvv,
                   trd, &lgr[(nabc+j)*nquadsq] ) )
        goto failure;
    }
    for ( j = 0, k = nab+nfunc_c+ii*s;
          j < s;
          j++, k++ ) {
      _g2h_GetSplDBasisCrossDerf ( domain, k, ii, fcomc, pv, pvv, pu, puu );
      _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g2h_TabSplD2LaplacianGraphf ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                   dddhfunc, nzc, j00, j01, lastomcknot, omcknots, fcomc,
                   lastpvknot, pvknots, pu, lastpvvknot, pvvknots, puu,
                   trd, &lgr[(nabc+s+j)*nquadsq] ) )
        goto failure;
    }

#ifdef DEBUG
DumpLapGrad ( i, nabc, 2*s, nquad, lgr );
#endif
          /* block A x A */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a, twok, twok)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k <= j;  k++, l1++ )
        auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                           0, nquad-1, 0, nquad-1, &lgr[k*nquadsq] );

          /* block A x B */
    auxmat1 = &SBMat[(nfunc_c+nfunc_d)*nfunc_b];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                             0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                             0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );

          /* block A x C */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, i)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] = _g2h_SplIntegralf ( nquad, jac,
                          0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                          i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
      }

          /* block A x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, hole_k+i)];
    auxmat2 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, hole_k+ii)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < s;  k++, l1++ ) {
        _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                           i00, i01, j00, j01, &lgr[(nabc+k)*nquadsq] );
        _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        auxmat2[l1] += _g2h_SplIntegralf ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lgr[j*nquadsq],
                           i00, i01, j00, j01, &lgr[(nabc+s+k)*nquadsq] );
      }

          /* block B x C */
    auxmat1 = &SBMat[i*r*nfunc_b];
    for ( j = l1 = 0;  j < r;  j++ ) {
      i00 = fkn[j/rr];  i01 = lkn[j/rr];
      j00 = fkn[j%rr];  j01 = lkn[j%rr];
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          auxmat1[l1] = _g2h_SplIntegralf ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nab+j)*nquadsq],
                            0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
    }

          /* block B x D */
    auxmat1 = &SBMat[(nfunc_c+i*s)*nfunc_b];
    auxmat2 = &SBMat[(nfunc_c+ii*s)*nfunc_b];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i10, &i11, &j10, &j11 );
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask ) {
          auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                             i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                             0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
          auxmat2[l1] += _g2h_SplIntegralf ( nquad, jac,
                             i10, i11, j10, j11, &lgr[(nabc+s+j)*nquadsq],
                             0, nquad-1, 0, nquad-1, &lgr[(nfunc_a+k)*nquadsq] );
        }
    }

          /* block C x C */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a, i, i)];
    for ( j = l1 = 0;  j < r;  j++ ) {
      i00 = fkn[j/rr];  i01 = lkn[j/rr];
      j00 = fkn[j%rr];  j01 = lkn[j%rr];
      for ( k = 0;  k <= j;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] = _g2h_SplIntegralf ( nquad, jac,
                          i00, i01, j00, j01, &lgr[(nab+j)*nquadsq],
                          i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
      }
    }

          /* block C x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+i, i)];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
      }
    }
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+ii, i)];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nab+k)*nquadsq] );
      }
    }

          /* block D x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+i, hole_k+i)];
    auxmat2 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+ii, hole_k+ii)];
    for ( j = l1 = l2 = 0;  j < s;  j++ ) {

      _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k <= j;  k++, l1++ ) {
        _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        auxmat1[l1] += _g2h_SplIntegralf ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nabc+k)*nquadsq] );
      }
      _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k <= j;  k++, l2++ ) {
        _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        auxmat2[l2] += _g2h_SplIntegralf ( nquad, jac,
                           i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                           i10, i11, j10, j11, &lgr[(nabc+s+k)*nquadsq] );
      }
    }
    if ( i == hole_k-1 ) {
            /* for i == hole_k-1 it is necessary to compute the transposition */
      auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                              2*hole_k-1, hole_k)];
      for ( k = l1 = 0;  k < s;  k++ ) {
        _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        for ( j = 0;  j < s;  j++, l1++ ) {
          _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                            &nzc, &i10, &i11, &j10, &j11 );
          auxmat1[l1] = _g2h_SplIntegralf ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nabc+k)*nquadsq],
                            i10, i11, j10, j11, &lgr[(nabc+s+j)*nquadsq] );
        }
      }
    }
    else {
      auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                              hole_k+ii, hole_k+i)];
      for ( j = l1 = 0;  j < s;  j++ ) {
        _g2h_FuncDSuppf ( hole_k, nk, m1, ii*s+j, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        for ( k = 0;  k < s;  k++, l1++ ) {
          _g2h_FuncDSuppf ( hole_k, nk, m1, i*s+k, i,
                            &nzc, &i10, &i11, &j10, &j11 );
          auxmat1[l1] = _g2h_SplIntegralf ( nquad, jac,
                            i00, i01, j00, j01, &lgr[(nabc+s+j)*nquadsq],
                            i10, i11, j10, j11, &lgr[(nabc+k)*nquadsq] );
        }
      }
    }
  }

        /* multiply the sums of function values at the quadrature knots */
        /* by the quadrature coefficient */
  qc = 1.0/(double)nquadsq;
  for ( i = 0; i < asize; i++ )
    SAMat[i] *= (float)qc;
  for ( i = 0; i < bsize; i++ )
    SBMat[i] *= (float)qc;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ComputeSplFormMatrixf*/

