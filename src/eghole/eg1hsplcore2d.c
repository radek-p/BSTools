
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

#ifdef DEBUG
#define QUAD_FACTOR 2
#else
#define QUAD_FACTOR 10
#endif

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_TabBSFuncDer2d ( int deg, int lastknot, const double *knots,
                              int i0, int i1,
                              int n, const double *tkn, int *fkn, int *lkn,
                              double *b, double *bt, double *btt )
{
  /* evaluate B-spline basis functions at a number of points  deg - degree,  */
  /* lastknot - number of the last knot, knots - knots, i0, i1 - numbers of  */
  /* the first and the last function to evaluate, n - number of (domain)     */
  /* points, tkn - the points (ascending), fkn, lkn - arrays in which for    */
  /* each function the numbers of the first and the last point, at which the */
  /* function is nonzero, are stored, b, bt, btt - arrays in which the values*/
  /* of the functions and their derivatives are to be stored                 */
  
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
      mbs_deBoorDer2C1d ( deg, lastknot, knots, p, tkn[j],
                          &b[l*n+j], &bt[l*n+j], &btt[l*n+j] );
    p[i] = 0.0;
 }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabBSFuncDer2d*/

void _g1h_TensDer2d ( double p, double pu, double puu,
                      double q, double qv, double qvv,
                      double *pq )
{
  pq[0] = pu*q;    pq[1] = p*qv;
  pq[2] = puu*q;   pq[3] = pu*qv;   pq[4] = p*qvv;
} /*_g1h_TensDerd*/

static void _g1h_TabSplCLaplaciand ( int nkn,
            int fknp, int lknp, const double *tp, const double *tpu,
            const double *tpuu,
            int fknq, int lknq, const double *tq, const double *tqv,
            const double *tqvv,
            const double *trd, double *lap )
{
  double pq[5];
  int    i, j, k, kk;

  memset ( lap, 0, nkn*nkn*sizeof(double) );
  for ( i = fknp; i <= lknp; i++ ) {
    for ( j = fknq, k = i*nkn+fknq, kk = 5*k;
          j <= lknq;
          j++, k++, kk += 5 ) {
      _g1h_TensDer2d ( tp[i], tpu[i], tpuu[i],
                       tq[j], tqv[j], tqvv[j], pq );
      lap[k] = pkn_ScalarProductd ( 5, &trd[kk], pq );
    }
  }
} /*_g1h_TabSplCLaplaciand*/

boolean _g1h_FuncDSuppd ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 )
{
  /* this function returns true if the function fn (in the block D) is  */
  /* nonzero in Omega_i; then the values assigned to *i0, *i1, *j0, *j1 */
  /* indicate the ranges of indexes of the quadrature knots, at which   */
  /* the function takes nonzero values                                  */
  int j, k, l, cn;

  cn = fn / (2*nk*m1);
  if ( cn == i ) {
    *j0 = 0;
    *j1 = l = (nk+1)*QUAD_FACTOR-1;
    fn %= 2*nk*m1;
    j = fn / 2;           /* which knot? */
    *nzc = k = fn % 2;    /* which order of cross derivative is nonzero? */
    if ( j <= 1 ) *i0 = 0;
    else *i0 = ((j-2)/m1+1)*QUAD_FACTOR;
    k = ((3-k+j)/m1+1)*QUAD_FACTOR-1;
    *i1 = min ( k, l );
  }
  else if ( (i+1) % hole_k == cn ) {
    *i0 = 0;
    *i1 = l = (nk+1)*QUAD_FACTOR-1;
    fn %= 2*nk*m1;
    j = fn / 2;           /* which knot? */
    *nzc = k = fn % 2;    /* which order of cross derivative is nonzero? */
    if ( j <= 1 ) *j0 = 0;
    else *j0 = ((j-2)/m1+1)*QUAD_FACTOR;
    k = ((3-k+j)/m1+1)*QUAD_FACTOR-1;
    *j1 = min ( k, l );
  }
  else
    return false;

  return true;
} /*_g1h_FuncDSuppd*/

static boolean _g1h_TabSplD1Laplaciand ( int nkn, const double *tkn,
                   const double *hfunc, const double *dhfunc,
                   const double *ddhfunc,
                   int nzc, int i0, int i1,
                   int lastomcknot, const double *omcknots, const double *fcomc,
                   int lastpvknot, const double *pvknots, const double *pv,
                   const double *trd,
                   double *lap )
{
  void   *sp;
  double *fp, bp[5];
  int    i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 6*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 6*nkn*sizeof(double) );
  ii = 6*i0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer2d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
           i1-i0+1, &tkn[i0], 6, &fp[ii], &fp[ii+1], &fp[ii+2] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pv,
         i1-i0+1, &tkn[i0], 6, &fp[ii+3], &fp[ii+4], &fp[ii+5] );

  memset ( lap, 0, nkn*nkn*sizeof(double) );
  for ( i = i0;  i <= i1;  i++, ii += 6 ) {
    for ( j = jj = 0, k = nkn*i, kk = 5*k;
          j < nkn;
          j++, jj += 4, k++, kk += 5 ) {
      switch ( nzc ) {
    default:
        bp[0] = fp[ii+1]*hfunc[jj]  + fp[ii+4]*hfunc[jj+2];
        bp[1] = fp[ii]*dhfunc[jj]   + fp[ii+3]*dhfunc[jj+2];
        bp[2] = fp[ii+2]*hfunc[jj]  + fp[ii+5]*hfunc[jj+2];
        bp[3] = fp[ii+1]*dhfunc[jj] + fp[ii+4]*dhfunc[jj+2];
        bp[4] = fp[ii]*ddhfunc[jj]  + fp[ii+3]*ddhfunc[jj+2];
        break;
    case 1:
        bp[0] = fp[ii+4]*hfunc[jj+2];
        bp[1] = fp[ii+3]*dhfunc[jj+2];
        bp[2] = fp[ii+5]*hfunc[jj+2];
        bp[3] = fp[ii+4]*dhfunc[jj+2];
        bp[4] = fp[ii+3]*ddhfunc[jj+2];
        break;
      }
      lap[k] = pkn_ScalarProductd ( 5, &trd[kk], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabSplD1Laplaciand*/

static boolean _g1h_TabSplD2Laplaciand ( int nkn, const double *tkn,
                   const double *hfunc, const double *dhfunc,
                   const double *ddhfunc,
                   int nzc, int j0, int j1,
                   int lastomcknot, const double *omcknots, const double *fcomc,
                   int lastpuknot, const double *puknots, const double *pu,
                   const double *trd,
                   double *lap )
{
  void   *sp;
  double *fp, bp[5];
  int    i, j, k, ii, jj, kk;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 6*nkn );
  if ( !fp )
    goto failure;

  memset ( fp, 0, 6*nkn*sizeof(double) );
  jj = 6*j0;
  if ( nzc == 0 )
    mbs_TabBSCurveDer2d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
           j1-j0+1, &tkn[j0], 6, &fp[jj], &fp[jj+1], &fp[jj+2] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpuknot, puknots, pu,
         j1-j0+1, &tkn[j0], 6, &fp[jj+3], &fp[jj+4], &fp[jj+5] );

  memset ( lap, 0, nkn*nkn*sizeof(double) );
  for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
    for ( j = j0, jj = 6*j0, k = nkn*i+j0, kk = 5*k;
          j <= j1;
          j++, jj += 6, k++, kk += 5 ) {
      switch ( nzc ) {
    default:
        bp[0] = fp[jj]*dhfunc[ii]   + fp[jj+3]*dhfunc[ii+2];
        bp[1] = fp[jj+1]*hfunc[ii]  + fp[jj+4]*hfunc[ii+2];
        bp[2] = fp[jj]*ddhfunc[ii]  + fp[jj+3]*ddhfunc[ii+2];
        bp[3] = fp[jj+1]*dhfunc[ii] + fp[jj+4]*dhfunc[ii+2];
        bp[4] = fp[jj+2]*hfunc[ii]  + fp[jj+5]*hfunc[ii+2];
        break;
    case 1:
        bp[0] = fp[jj+3]*dhfunc[ii+2];
        bp[1] = fp[jj+4]*hfunc[ii+2];
        bp[2] = fp[jj+3]*ddhfunc[ii+2];
        bp[3] = fp[jj+4]*dhfunc[ii+2];
        bp[4] = fp[jj+5]*hfunc[ii+2];
        break;
      }
      lap[k] = pkn_ScalarProductd ( 5, &trd[kk], bp );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabSplD2Laplaciand*/

#ifdef DEBUG
static boolean dumpspl = false;

double _g1h_SplIntegrald ( int nkn, double *jac,
                          int i00, int i01, int j00, int j01, double *func0,
                          int i10, int i11, int j10, int j11, double *func1 )
{
  FILE   *f;
  int    i, j, k;
  long double s, ss;

  if ( dumpspl )
    f = fopen ( "intspl.txt", "w+" );
        /* two-level summing is supposed to decrease the effect */
        /* of rounding errors */
  i00 = max ( i00, i10 );  i01 = min ( i01, i11 );
  j00 = max ( j00, j10 );  j01 = min ( j01, j11 );
  for ( i = i00, ss = 0.0; i <= i01; i++ ) {
    for ( j = j00, k = nkn*i+j00, s = 0.0;  j <= j01;  j++, k++ ) {
      if ( dumpspl )
        fprintf ( f, "%d, %e, %e, %e\n", k, func0[k], func1[k], jac[k] );
      s += func0[k]*func1[k]*jac[k];
    }
    ss += s;
  }
  if ( dumpspl ) {
    fprintf ( f, "\n" );
    fclose ( f );
  }
  return ss;
} /*_g1h_SplIntegrald*/
#else
static double _g1h_SplIntegrald ( int nkn, double *jac,
                          int i00, int i01, int j00, int j01, double *func0,
                          int i10, int i11, int j10, int j11, double *func1 )
{
  int    i, j, k;
  long double s, ss;

        /* two-level summing is supposed to decrease the effect */
        /* of rounding errors */
  i00 = max ( i00, i10 );  i01 = min ( i01, i11 );
  j00 = max ( j00, j10 );  j01 = min ( j01, j11 );
  for ( i = i00, ss = 0.0; i <= i01; i++ ) {
    for ( j = j00, k = nkn*i+j00, s = 0.0;  j <= j01;  j++, k++ )
      s += func0[k]*func1[k]*jac[k];
    ss += s;
  }
  return (double)ss;
} /*_g1h_SplIntegrald*/
#endif
/* ///////////////////////////////////////////////////////////////////////// */
#ifdef DEBUG
void DumpJacobian ( int i, int nkn, const double *jac )
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

void DumpLap ( int i, int f0, int nfunc, int nkn, const double *lgr )
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
      fprintf ( f, "%4d: %e\n", k, lgr[l] );
    fprintf ( f, "\n" );
  }
  fprintf ( f, "\n" );
  fclose ( f );
} /*DumpLap*/

void DumpTBS ( int nf, int nkn,
               const double *tp, const double *tpu, const double *tpuu )
{
  FILE   *f;
  int    i, j, k, l;
  double pq[5];

  f = fopen ( "tbs.txt", "w+" );
  for ( i = 0; i < nf; i++ )
    for ( j = 0; j < nf; j++ ) {
      fprintf ( f, "fn:%d\n", nf*i+j );
      for ( k = 0; k < nkn; k++ )
        for ( l = 0; l < nkn; l++ ) {
          _g1h_TensDer2f ( tp[i*nkn+k], tpu[i*nkn+k], tpuu[i*nkn+k],
                           tp[j*nkn+l], tpu[j*nkn+l], tpuu[j*nkn+l], pq );
          fprintf ( f,
            " pu:%e, pv:%e, puu:%e, puv:%e, pvv:%e\n",
            pq[0], pq[1], pq[2], pq[3], pq[4] );
        }
      fprintf ( f, "\n" );
    }
  fclose ( f );
} /*DumpTBS*/
#endif
/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_ComputeSplFormMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd   *privateG;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nab, nabc;
  int      lastomcknot, lastpvknot;
  double   *omcknots, *pvknots;
  int      rr, r, s, asize, bsize;
  double   *SAMat, *SBMat, *auxmat1, *auxmat2;
  unsigned short *support_b, mask;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *fcomc, *pv, *pu;
  double   *trd, *jac, *tbs, *tbst, *tbstt;
  int      *fkn, *lkn;
  double   *lap;
  int      nk, m1, m2, m, nquad, nquadsq;
  int      i, ii, j, k, l1, l2, twok;
  int      i00, i01, j00, j01, i10, i11, j10, j11, nzc;
  double   qc;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  twok = 2*hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
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

        /* compute the size of the arrays and allocate memory */
  rr = G1H_FINALDEG-3+nk*m2;
  r = rr*rr;    /* == nfunc_c/hole_k */
  s = 2*nk*m1;  /* == nfunc_d/hole_k */
  nab = nfunc_a+nfunc_b;
  nabc = nab+r;
  asize = pkn_Block2ArraySize ( hole_k, r, s, nfunc_a );
  bsize = nfunc_b*(nfunc_a+nfunc_c+nfunc_d);
  SAMat = sprivate->SAMat = malloc ( asize*sizeof(double) );
  SBMat = sprivate->SBMat = malloc ( bsize*sizeof(double) );
  if ( !SAMat || !SBMat )
    goto failure;

        /* fix the number of quadrature knots, perhaps later to be */
        /* introduced by the application via OptionProc */
  nquad = QUAD_FACTOR*(nk+1);
  nquadsq = nquad*nquad;
  m = rr*nquad;

        /* allocate the scratch memory */
  tkn = pkv_GetScratchMemd ( 13*nquad );
  jac = pkv_GetScratchMemd ( nquadsq );
  trd = pkv_GetScratchMemd ( 5*nquadsq );
  fkn = pkv_GetScratchMem ( 2*rr*sizeof(int) );
  lap = pkv_GetScratchMemd ( (nabc+2*s)*nquadsq );
  tbs = pkv_GetScratchMemd ( 3*m );
  fcomc = pkv_GetScratchMemd ( lastomcknot+2*lastpvknot-
                               (G1_CROSS00DEG+2*G1_CROSS01DEG) );
  if ( !tkn || !jac || !trd || !fkn || !lap || !tbs || !fcomc )
    goto failure;
  hfunc = &tkn[nquad];         dhfunc = &hfunc[4*nquad];
  ddhfunc = &dhfunc[4*nquad];
  tbst = &tbs[m];              tbstt = &tbst[m];
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  lkn = &fkn[rr];

        /* prepare the evaluation of basis functions */
  _gh_PrepareTabKnotsd ( nquad, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, nquad, tkn, hfunc, dhfunc, ddhfunc ) )
    goto failure;

  if ( !_g1h_TabBSFuncDer2d ( G1H_FINALDEG, sprivate->lastcknot, sprivate->cknots,
                              2, G1H_FINALDEG+nk*m2-2, nquad, tkn, fkn, lkn,
                              tbs, tbst, tbstt ) )
    goto failure;

        /* compute and sum the integrals in subsequent areas omega_i */
  memset ( SAMat, 0, asize*sizeof(double) );
  memset ( SBMat, 0, bsize*sizeof(double) );
  for ( i = 0, ii = 1, mask = 0x0001;
        i < hole_k;
        i++, ii = (ii+1) % hole_k, mask = (unsigned short)(mask << 1) ) {

          /* evaluate the Jacobian and the matrix of derivatives of the */
          /* domain patch at the quadrature knots */
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    _g1h_TabDiPatchJac2d ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                           c00, c01, c10, c11, d00, d01, d10, d11,
                           jac, trd );

#ifdef DEBUG
DumpJacobian ( i, nquad, jac );
#endif
          /* evaluate the Laplacian gradient of the A block basis functions */
          /* at the quadrature knots */
    for ( j = 0; j < nfunc_a; j++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabLaplacian0d ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                            fc00, fc01, fd00, fd01, trd,
                            &lap[j*nquadsq] );
    }

          /* evaluate the Laplacian gradient of the B block basis functions */
          /* at the quadrature knots */
    for ( j = 0; j < nfunc_b; j++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesd ( domain, j, i,
                                  &fc00, &fc01, &fc10, &fc11,
                                  &fd00, &fd01, &fd10, &fd11 );
        _g1h_TabLaplaciand ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                             fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                             trd, &lap[(nfunc_a+j)*nquadsq] );
      }

          /* evaluate the Laplacian gradient of the C block basis functions */
          /* at the quadrature knots */
    for ( j = l1 = 0;  j < rr;  j++ )
      for ( k = 0;  k < rr;  k++, l1++ )
        _g1h_TabSplCLaplaciand ( nquad,
                  fkn[j], lkn[j], &tbs[j*nquad], &tbst[j*nquad], &tbstt[j*nquad],
                  fkn[k], lkn[k], &tbs[k*nquad], &tbst[k*nquad], &tbstt[k*nquad],
                  trd, &lap[(nab+l1)*nquadsq] );

          /* evaluate the Laplacian gradient of the D block basis functions */
          /* at the quadrature knots; do it only for the functions nonzero  */
          /* in Omega_i; j is the function number within the subblock, k is */
          /* the complete function number */
    for ( j = 0, k = nab+nfunc_c+i*s;
          j < s;
          j++, k++ ) {
      _g1h_GetSplDBasisCrossDerd ( domain, k, i, fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_TabSplD1Laplaciand ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                   nzc, i00, i01, lastomcknot, omcknots, fcomc,
                   lastpvknot, pvknots, pv,
                   trd, &lap[(nabc+j)*nquadsq] ) )
        goto failure;
    }
    for ( j = 0, k = nab+nfunc_c+ii*s;
          j < s;
          j++, k++ ) {
      _g1h_GetSplDBasisCrossDerd ( domain, k, ii, fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1h_TabSplD2Laplaciand ( nquad, tkn, hfunc, dhfunc, ddhfunc,
                   nzc, j00, j01, lastomcknot, omcknots, fcomc,
                   lastpvknot, pvknots, pu,
                   trd, &lap[(nabc+s+j)*nquadsq] ) )
        goto failure;
    }

#ifdef DEBUG
DumpLap ( i, nabc, 2*s, nquad, lap );
#endif
          /* block A x A */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a, twok, twok)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k <= j;  k++, l1++ )
        auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lap[j*nquadsq],
                           0, nquad-1, 0, nquad-1, &lap[k*nquadsq] );

          /* block A x B */
    auxmat1 = &SBMat[(nfunc_c+nfunc_d)*nfunc_b];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                             0, nquad-1, 0, nquad-1, &lap[j*nquadsq],
                             0, nquad-1, 0, nquad-1, &lap[(nfunc_a+k)*nquadsq] );

          /* block A x C */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, i)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] = _g1h_SplIntegrald ( nquad, jac,
                          0, nquad-1, 0, nquad-1, &lap[j*nquadsq],
                          i10, i11, j10, j11, &lap[(nab+k)*nquadsq] );
      }

          /* block A x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, hole_k+i)];
    auxmat2 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            2*hole_k, hole_k+ii)];
    for ( j = l1 = 0;  j < nfunc_a;  j++ )
      for ( k = 0;  k < s;  k++, l1++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lap[j*nquadsq],
                           i00, i01, j00, j01, &lap[(nabc+k)*nquadsq] );
        _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        auxmat2[l1] += _g1h_SplIntegrald ( nquad, jac,
                           0, nquad-1, 0, nquad-1, &lap[j*nquadsq],
                           i00, i01, j00, j01, &lap[(nabc+s+k)*nquadsq] );
      }

          /* block B x C */
    auxmat1 = &SBMat[i*r*nfunc_b];
    for ( j = l1 = 0;  j < r;  j++ ) {
      i00 = fkn[j/rr];  i01 = lkn[j/rr];
      j00 = fkn[j%rr];  j01 = lkn[j%rr];
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask )
          auxmat1[l1] = _g1h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lap[(nab+j)*nquadsq],
                            0, nquad-1, 0, nquad-1, &lap[(nfunc_a+k)*nquadsq] );
    }

          /* block B x D */
    auxmat1 = &SBMat[(nfunc_c+i*s)*nfunc_b];
    auxmat2 = &SBMat[(nfunc_c+ii*s)*nfunc_b];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i10, &i11, &j10, &j11 );
      for ( k = 0;  k < nfunc_b;  k++, l1++ )
        if ( support_b[k] & mask ) {
          auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                             i00, i01, j00, j01, &lap[(nabc+j)*nquadsq],
                             0, nquad-1, 0, nquad-1, &lap[(nfunc_a+k)*nquadsq] );
          auxmat2[l1] += _g1h_SplIntegrald ( nquad, jac,
                             i10, i11, j10, j11, &lap[(nabc+s+j)*nquadsq],
                             0, nquad-1, 0, nquad-1, &lap[(nfunc_a+k)*nquadsq] );
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
        auxmat1[l1] = _g1h_SplIntegrald ( nquad, jac,
                          i00, i01, j00, j01, &lap[(nab+j)*nquadsq],
                          i10, i11, j10, j11, &lap[(nab+k)*nquadsq] );
      }
    }

          /* block C x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+i, i)];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lap[(nabc+j)*nquadsq],
                           i10, i11, j10, j11, &lap[(nab+k)*nquadsq] );
      }
    }
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+ii, i)];
    for ( j = l1 = 0;  j < s;  j++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k < r;  k++, l1++ ) {
        i10 = fkn[k/rr];  i11 = lkn[k/rr];
        j10 = fkn[k%rr];  j11 = lkn[k%rr];
        auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lap[(nabc+s+j)*nquadsq],
                           i10, i11, j10, j11, &lap[(nab+k)*nquadsq] );
      }
    }

          /* block D x D */
    auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+i, hole_k+i)];
    auxmat2 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                            hole_k+ii, hole_k+ii)];
    for ( j = l1 = l2 = 0;  j < s;  j++ ) {

      _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k <= j;  k++, l1++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        auxmat1[l1] += _g1h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lap[(nabc+j)*nquadsq],
                           i10, i11, j10, j11, &lap[(nabc+k)*nquadsq] );
      }
      _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                        &nzc, &i00, &i01, &j00, &j01 );
      for ( k = 0;  k <= j;  k++, l2++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+k, i,
                          &nzc, &i10, &i11, &j10, &j11 );
        auxmat2[l2] += _g1h_SplIntegrald ( nquad, jac,
                           i00, i01, j00, j01, &lap[(nabc+s+j)*nquadsq],
                           i10, i11, j10, j11, &lap[(nabc+s+k)*nquadsq] );
      }
    }
    if ( i == hole_k-1 ) {
            /* for i == hole_k-1 it is necessary to compute the transposition */
      auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                              2*hole_k-1, hole_k)];
      for ( k = l1 = 0;  k < s;  k++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        for ( j = 0;  j < s;  j++, l1++ ) {
          _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                            &nzc, &i10, &i11, &j10, &j11 );
          auxmat1[l1] = _g1h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lap[(nabc+k)*nquadsq],
                            i10, i11, j10, j11, &lap[(nabc+s+j)*nquadsq] );
        }
      }
    }
    else {
      auxmat1 = &SAMat[pkn_Block2FindBlockPos(hole_k, r, s, nfunc_a,
                                              hole_k+ii, hole_k+i)];
      for ( j = l1 = 0;  j < s;  j++ ) {
        _g1h_FuncDSuppd ( hole_k, nk, m1, ii*s+j, i,
                          &nzc, &i00, &i01, &j00, &j01 );
        for ( k = 0;  k < s;  k++, l1++ ) {
          _g1h_FuncDSuppd ( hole_k, nk, m1, i*s+k, i,
                            &nzc, &i10, &i11, &j10, &j11 );
          auxmat1[l1] = _g1h_SplIntegrald ( nquad, jac,
                            i00, i01, j00, j01, &lap[(nabc+s+j)*nquadsq],
                            i10, i11, j10, j11, &lap[(nabc+k)*nquadsq] );
        }
      }
    }
  }

        /* multiply the sums of function values at the quadrature knots */
        /* by the quadrature coefficient */
  qc = 1.0/(double)nquadsq;
  for ( i = 0; i < asize; i++ )
    SAMat[i] *= qc;
  for ( i = 0; i < bsize; i++ )
    SBMat[i] *= qc;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeSplFormMatrixd*/

