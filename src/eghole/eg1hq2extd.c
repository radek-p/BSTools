
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

/*#define DEBUG*/

/* ///////////////////////////////////////////////////////////////////////// */
static void _g1h_TabBezierPolynomialsDer3d ( int nkn, const double *kn,
                                             double *bezp, double *dbezp,
                                             double *ddbezp, double *dddbezp )
{
  int    i, j, k, b;
  double s, t, u, v, binom[G1H_FINALDEG-2], bpoly[G1H_FINALDEG];

        /* compute the binomial coefficients */
  binom[0] = 1.0;
  b = G1H_FINALDEG-3;
  for ( i = 1; i <= G1H_FINALDEG-3; i++ ) {
    binom[i] = (double)b;
    b = (b*(G1H_FINALDEG-3-i))/(i+1);
  }

  for ( j = k = 0;  j < nkn;  j++, k += G1H_FINALDEG-3 ) {
    t = kn[j];
    s = (double)(1.0-t);

        /* evaluate the polynomials of degree G1H_FINALDEG-3 at the knot t */
    bpoly[0] = bpoly[G1H_FINALDEG-1] = 0.0;
    memcpy ( &bpoly[1], binom, (G1H_FINALDEG-2)*sizeof(double) );
    u = t;  v = s;
    bpoly[2] *= u;  bpoly[G1H_FINALDEG-3] *= v;
    for ( i = 2; i <= G1H_FINALDEG-3; i++ )
      { u *= t;  v *= s;  bpoly[i+1] *= u;  bpoly[G1H_FINALDEG-2-i] *= v; }

        /* evaluate the third order derivatives */
    b = G1H_FINALDEG*(G1H_FINALDEG-1)*(G1H_FINALDEG-2);
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      dddbezp[k+i] =
        (double)(b*((bpoly[i]+3.0*bpoly[i+2])-(3.0*bpoly[i+1]+bpoly[i+3])));

        /* evaluate the polynomials of degree G1H_FINALDEG-2 */
        /* - de Casteljau algorithm */
    bpoly[0] = s*bpoly[1];
    for ( i = 1; i < G1H_FINALDEG-2; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];
    bpoly[G1H_FINALDEG-2] = t*bpoly[G1H_FINALDEG-2];

        /* evaluate the second order derivatives */
    b = G1H_FINALDEG*(G1H_FINALDEG-1);
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      ddbezp[k+i] = (double)(b*((bpoly[i]+bpoly[i+2])-2.0*bpoly[i+1]));

        /* evaluate the polynomials of degree G1H_FINALDEG-1 */
    for ( i = 0; i <= G1H_FINALDEG-3; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];

        /* evaluate the first order derivatives */
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      dbezp[k+i] = (double)G1H_FINALDEG*(bpoly[i]-bpoly[i+1]);

        /* evaluate the polynomials of degree G1H_FINALDEG */
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      bezp[k+i] = s*bpoly[i+1] + t*bpoly[i];
  }
} /*_g1h_TabBezierPolynomialsDer3d*/

static void _g1h_TabTensBezd ( int nkn, const double *bez1, const double *bez2,
                               double *tbez )
{
  int i, j, k, l, m;

/* !!! this procedure is identical to that in eg1hextd.c !!! */
  for ( i = m = 0; i < G1H_FINALDEG-3; i++ )
    for ( j = 0; j < G1H_FINALDEG-3; j++ )  
      for ( k = 0; k < nkn; k++ )
        for ( l = 0;  l < nkn;  l++, m++ )
          tbez[m] = bez1[k*(G1H_FINALDEG-3)+i]*bez2[l*(G1H_FINALDEG-3)+j];
} /*_g1h_TabTensBezd*/

boolean _g1h_TabTensBezPolyDer3d ( int nkn, const double *tkn,
             double *tbez, double *tbezu, double *tbezv,
             double *tbezuu, double *tbezuv, double *tbezvv,
             double *tbezuuu, double *tbezuuv, double *tbezuvv, double *tbezvvv )
{
  void   *sp;
  double *bezp, *dbezp, *ddbezp, *dddbezp;

  sp = pkv_GetScratchMemTop ();
  bezp = pkv_GetScratchMemd ( (G1H_FINALDEG-3)*nkn*4 );
  if ( !bezp )
    goto failure;
  dbezp = &bezp[(G1H_FINALDEG-3)*nkn];
  ddbezp = &dbezp[(G1H_FINALDEG-3)*nkn];
  dddbezp = &ddbezp[(G1H_FINALDEG-3)*nkn];

  _g1h_TabBezierPolynomialsDer3d ( G1_NQUAD, tkn, bezp, dbezp, ddbezp, dddbezp );
  if ( tbez )
    _g1h_TabTensBezd ( nkn, bezp, bezp, tbez );
  _g1h_TabTensBezd ( nkn, dbezp, bezp, tbezu );
  _g1h_TabTensBezd ( nkn, bezp, dbezp, tbezv );
  _g1h_TabTensBezd ( nkn, ddbezp, bezp, tbezuu );
  _g1h_TabTensBezd ( nkn, dbezp, dbezp, tbezuv );
  _g1h_TabTensBezd ( nkn, bezp, ddbezp, tbezvv );
  _g1h_TabTensBezd ( nkn, dddbezp, bezp, tbezuuu );
  _g1h_TabTensBezd ( nkn, ddbezp, dbezp, tbezuuv );
  _g1h_TabTensBezd ( nkn, dbezp, ddbezp, tbezuvv );
  _g1h_TabTensBezd ( nkn, bezp, dddbezp, tbezvvv );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabTensBezPolyDer3d*/

static void _g1h_TabLaplacianGrad00d ( int nkn,
                const double *tbezu, const double *tbezv,
                const double *tbezuu, const double *tbezuv, const double *tbezvv,
                const double *tbezuuu, const double *tbezuuv,
                const double *tbezuvv, const double *tbezvvv,
                const double *trd, vector2d *lapgrad )
{
/* !!! the same as _g2h_TabLaplacianGrad00d in the file eg2hextd.c !!! */
/* to be replaced in both places by one nonstatic function */
  int i, k;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 18 ) {
    lapgrad[i].x = trd[k+0]*tbezu[i] + trd[k+1]*tbezv[i] +
                   trd[k+2]*tbezuu[i] + trd[k+3]*tbezuv[i] +
                   trd[k+4]*tbezvv[i] +
                   trd[k+5]*tbezuuu[i] + trd[k+6]*tbezuuv[i] +
                   trd[k+7]*tbezuvv[i] + trd[k+8]*tbezvvv[i];
    lapgrad[i].y = trd[k+9]*tbezu[i] + trd[k+10]*tbezv[i] +
                   trd[k+11]*tbezuu[i] + trd[k+12]*tbezuv[i] +
                   trd[k+13]*tbezvv[i] +
                   trd[k+14]*tbezuuu[i] + trd[k+15]*tbezuuv[i] +
                   trd[k+16]*tbezuvv[i] + trd[k+17]*tbezvvv[i];
  }
} /*_g1h_TabLaplacianGrad00d*/

boolean _g1h_TabLaplacianJump00d ( int nkn, const double *tkn, int fni,
                    const double *trdc00, const double *trdc10,
                    const double *trdd00, const double *trdd10,
                    double *lapc00, double *lapc10, double *lapd00, double *lapd10 )
{
  void   *sp;
  double *b, p, pu, pv, puu, puv, pvv, lap;
  int    fi, fj, i, j;

  sp = pkv_GetScratchMemTop ();
  b = pkv_GetScratchMemd ( (G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !b )
    goto failure;
  memset ( b, 0, (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(double) );
  fi = fni / (G1H_FINALDEG-3) + 2;
  fj = fni % (G1H_FINALDEG-3) + 2;
  b[fi*(G1H_FINALDEG+1)+fj] = 1.0;

/* !!! this is not as it ought to be!, to be replaced later !!! */
  if ( fi == 2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 0.0, tkn[i],
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdd00[j]*pu + trdd00[j+1]*pv +
            trdd00[j+2]*puu + trdd00[j+3]*puv + trdd00[j+4]*pvv;
      lapd00[i] = -lap;
    }
  }
  else
    memset ( lapd00, 0, nkn*sizeof(double) );
  if ( fi == G1H_FINALDEG-2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 1.0, tkn[i],
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdd10[j]*pu + trdd10[j+1]*pv +
            trdd10[j+2]*puu + trdd10[j+3]*puv + trdd10[j+4]*pvv;
      lapd10[i] = -lap;
    }
  }
  else
    memset ( lapd10, 0, nkn*sizeof(double) );

  if ( fj == 2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[i], 0.0,
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdc00[j]*pu + trdc00[j+1]*pv +
            trdc00[j+2]*puu + trdc00[j+3]*puv + trdc00[j+4]*pvv;
      lapc00[i] = lap;
    }
  }
  else
    memset ( lapc00, 0, nkn*sizeof(double) );
  if ( fj == G1H_FINALDEG-2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[i], 1.0,
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdc10[j]*pu + trdc10[j+1]*pv +
            trdc10[j+2]*puu + trdc10[j+3]*puv + trdc10[j+4]*pvv;
      lapc10[i] = -lap;
    }
  }
  else
    memset ( lapc10, 0, nkn*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabLaplacianJump00d*/

static double _g1h_Q2Integral0d ( double *jac, double *lapj1, double *lapj2 )
{
  double s;
  int    i;

  for ( i = 0, s = 0.0;  i < G1_NQUAD; i++ )
    s += jac[i]*lapj1[i]*lapj2[i];
  s /= (double)G1_NQUAD;
  return (double)s;
} /*_g1h_Q2Integral0d*/

boolean g1h_Q2ExtComputeFormMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd    *privateG;
  G1HolePrivateRecd   *privateG1;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;
  int      i, j, k, l, n, m, fn;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  double   *atkn, *ahfunc, *adhfunc, *addhfunc;
  double   *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *tbezuuu, *tbezuuv, *tbezuvv, *tbezvvv;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *ec00, *ec01, *ec10, *ec11, *ed00, *ed01, *ed10, *ed11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *trd, *lapj, *lapjb, *lapjc, *jac, *amat, *bmat, *Aij;
  vector2d *lgr, *lgrb, *lgrc;
  unsigned short *support_b, supp;
  int      option, ndata, *idata, Asize, Bsize;
  double   *fdata, *eicp1, *eicp2, C;
  double   s;
  point2d  *sicp;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisd ( domain ) )
      goto failure;
  }
  hole_k    = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a   = privateG1->nfunc_a;
  nfunc_b   = privateG1->nfunc_b;
  nfunc_c   = G1_DBDIM*hole_k;
  support_b = privateG->support_b;

  Asize = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
  Bsize = (nfunc_a+nfunc_c)*nfunc_b;
  if ( !(amat = privateG1->Q2EAMat) )
    amat = privateG1->Q2EAMat = malloc ( Asize*sizeof(double) );
  if ( !(bmat = privateG1->Q2EBMat) )
    bmat = privateG1->Q2EBMat = malloc ( Bsize*sizeof(double) );
  if ( !amat || !bmat )
    goto failure;
  if ( !privateG->diam )
    privateG->diam = gh_DomainDiamd ( domain );
 
  n = hole_k*G1_NQUADSQ;
  m = G1_DBDIM*G1_NQUADSQ;

  tkn = pkv_GetScratchMemd ( 17*G1_NQUAD );
  jac = pkv_GetScratchMemd ( n );
  trd = pkv_GetScratchMemd ( 18*n );
  lgr = (vector2d*)pkv_GetScratchMem (
                       ((nfunc_a+1)*n+nfunc_c*G1_NQUADSQ)*sizeof(vector2d) );
  tbezu = pkv_GetScratchMemd ( 9*m );
  if ( !tkn || !jac || !trd || !lgr|| !tbezu )
    goto failure;
  hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[4*G1_NQUAD];
  ddhfunc = &dhfunc[4*G1_NQUAD];  dddhfunc = &ddhfunc[4*G1_NQUAD];
  tbezv = &tbezu[m];           tbezuu = &tbezv[m];
  tbezuv = &tbezuu[m];         tbezvv = &tbezuv[m];
  tbezuuu = &tbezvv[m];        tbezuuv = &tbezuuu[m];
  tbezuvv = &tbezuuv[m];       tbezvvv = &tbezuvv[m];
  lgrc = &lgr[nfunc_a*n];
  lgrb = &lgrc[nfunc_c*G1_NQUADSQ];

  _gh_PrepareTabKnotsd ( G1_NQUAD, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3d ( 0.0, 1.0, G1_NQUAD, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;
  _g1h_TabTensBezPolyDer3d ( G1_NQUAD, tkn, NULL, tbezu, tbezv, tbezuu, tbezuv,
                             tbezvv, tbezuuu, tbezuuv, tbezuvv, tbezvvv );

  memset ( amat, 0, Asize*sizeof(double) );
  memset ( bmat, 0, Bsize*sizeof(double) );

        /* integrate the Laplacian gradients */
  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_Q2TabDiPatchJac3d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_Q2TabLaplacianGrad0d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                  fc00, fc01, fd00, fd01,
                                  &trd[i*G1_NQUADSQ*18],
                                  &lgr[(fn*hole_k+i)*G1_NQUADSQ] );
    }
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < G1_DBDIM; j++ )
      _g1h_TabLaplacianGrad00d ( G1_NQUAD, &tbezu[j*G1_NQUADSQ], &tbezv[j*G1_NQUADSQ],
            &tbezuu[j*G1_NQUADSQ], &tbezuv[j*G1_NQUADSQ], &tbezvv[j*G1_NQUADSQ],
            &tbezuuu[j*G1_NQUADSQ], &tbezuuv[j*G1_NQUADSQ], &tbezuvv[j*G1_NQUADSQ],
            &tbezvvv[j*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18], &lgrc[(i*G1_DBDIM+j)*G1_NQUADSQ] );

  Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                     hole_k-1, hole_k-1)];
  for ( i = 0; i < nfunc_a; i++ )
    for ( j = 0; j <= i; j++ )
      Aij[pkn_SymMatIndex(G1_DBDIM+i,G1_DBDIM+j)] =
        _g2h_Integrald ( hole_k, G1_NQUADSQ, jac,
                         0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  for ( k = 0; k < hole_k; k++ ) {
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       k, k)];
    for ( i = l = 0;  i < G1_DBDIM;  i++ )
      for ( j = 0;  j <= i;  j++, l++ )
        Aij[l] = _g2h_Integral0d ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
                     &lgrc[k*m+i*G1_NQUADSQ], &lgrc[k*m+j*G1_NQUADSQ] );
  }

  for ( k = 0; k < hole_k; k++ ) {
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       hole_k-1, k)];
    for ( i = 0; i < G1_DBDIM; i++ )
      for ( j = 0; j < nfunc_a; j++ ) {
        if ( k < hole_k-1 )
          l = (j+G1_DBDIM)*G1_DBDIM+i;
        else
          l = pkn_SymMatIndex(G1_DBDIM+j,i);
        Aij[l] = _g2h_Integral0d ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
                          &lgrc[k*m+i*G1_NQUADSQ], &lgr[j*n+k*G1_NQUADSQ] );
      }
  }

        /* right-hand side matrix */
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        _g1h_GetBFBPatchCurvesd ( domain, j, k,
                      &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        _g1h_Q2TabLaplacianGradd ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                   fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                                   &trd[k*G1_NQUADSQ*18], &lgrb[k*G1_NQUADSQ] );
      }

    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        for ( i = 0; i < G1_DBDIM; i++ )
          bmat[(k*G1_DBDIM+i)*nfunc_b+j] = _g2h_Integral0d ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
                                     &lgrc[k*m+i*G1_NQUADSQ], &lgrb[k*G1_NQUADSQ] );
      }

    for ( i = 0; i < nfunc_a; i++ )
      bmat[(k*G1_DBDIM+i)*nfunc_b+j] = _g2h_Integrald ( hole_k, G1_NQUADSQ, jac,
                                    0xFFFF, &lgr[i*n], support_b[j], lgrb );
  }

        /* integrate the Laplacian jumps */
  option = privateG1->GetOption ( domain, G1HQUERY_Q2_FORM_CONSTANT, 0,
                                  &ndata, &idata, &fdata );
  switch ( option ) {
case G1H_DEFAULT:
    privateG1->C1e = C = 1.0;
                       /* no idea whether this default value makes any sense */
    break;
case G1H_Q2_USE_SUPPLIED_CONSTANT:
    privateG1->C1e = C = *fdata;   
    break;
default:  
    goto failure;
  }
  C *= 5.0/privateG->diam;  /* 3.0 = sqrt(G1_DBDIM)+1 */

        /* the arrays allocated for the previous stage are reused */
        /* their length is sufficient */
  lapj = &lgr[0].x;
  lapjb = &lapj[3*G1_NQUAD*hole_k*nfunc_a];
  lapjc = &lapjb[3*G1_NQUAD*hole_k];
        /* two additional arrays are however needed */
  atkn = pkv_GetScratchMemd ( 26 );
  sicp = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !atkn || !sicp )
    goto failure;
  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  eicp1 = (double*)sicp;  eicp2 = &eicp1[16];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
          /* compute the curve Jacobians */
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_TabCurveJacobiand ( G1_CROSS00DEG, c00, G1_NQUAD, tkn, &jac[3*i*G1_NQUAD] );
    _g1h_TabCurveJacobiand ( G1_CROSS10DEG, c10, G1_NQUAD, tkn, &jac[(3*i+1)*G1_NQUAD] );
    _g1h_TabCurveJacobiand ( G1_CROSS10DEG, d10, G1_NQUAD, tkn, &jac[(3*i+2)*G1_NQUAD] );
          /* compute the coefficients for computing the Laplacians */
    if ( !_g1h_TabCurveLapCoeff0d ( c00, c01, c10, c11, d00, d01, d10, d11,
                            G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc,   
                            &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                            &trd[(30*j+5)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchd ( domain, j, 1, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( sicp, G1_NQUAD, tkn, &trd[(30*i+15)*G1_NQUAD] );
    if ( !_gh_FindDomSurrndPatchd ( domain, i, 2, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( sicp, G1_NQUAD, tkn, &trd[(30*i+25)*G1_NQUAD] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ ) {
    _g1h_GetBFAPatchCurvesd ( domain, fn, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1h_Q2TabLaplacianJump0d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                          atkn, ahfunc, adhfunc, addhfunc,
                          ec00, ec01, ed00, ed01, &trd[(30*i+5)*G1_NQUAD],
                          fc00, fc01, fd00, fd01, &trd[30*i*G1_NQUAD],
                          &trd[(30*i+10)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD],
                          &lapj[3*(fn*hole_k+i)*G1_NQUAD],
                          &lapj[(3*(fn*hole_k+i)+1)*G1_NQUAD],
                          &lapj[(3*(fn*hole_k+i)+2)*G1_NQUAD] ) )
        goto failure;
      ec00 = fc00;  ec01 = fc01;  ed00 = fd00;  ed01 = fd01;
    }
  }
  for ( k = 0;  k < hole_k;  k++ ) {
    j = (k+1) % hole_k;
    for ( i = 0;  i < G1_DBDIM;  i++ ) {
      if ( !_g1h_TabLaplacianJump00d ( G1_NQUAD, tkn, i,
                    &trd[30*k*G1_NQUAD], &trd[(30*k+10)*G1_NQUAD],
                    &trd[(30*j+5)*G1_NQUAD], &trd[(30*k+20)*G1_NQUAD],
                    &lapjc[(k*G1_DBDIM+i)*4*G1_NQUAD],
                    &lapjc[((k*G1_DBDIM+i)*4+1)*G1_NQUAD],
                    &lapjc[((k*G1_DBDIM+i)*4+2)*G1_NQUAD],
                    &lapjc[((k*G1_DBDIM+i)*4+3)*G1_NQUAD] ) )
        goto failure;
    }
  }

  Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                     hole_k-1, hole_k-1)];
  for ( i = 0; i < nfunc_a; i++ ) {
    l = pkn_SymMatIndex ( i+G1_DBDIM, G1_DBDIM );
    for ( j = 0; j <= i; j++ )
      Aij[l+j] += C*_g1h_Q2Integrald ( hole_k, G1_NQUAD, jac,
                                       0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                       0xFFFF, &lapj[j*hole_k*3*G1_NQUAD] );
  }

  for ( k = 0; k < hole_k; k++ ) {
    j = (k+1) % hole_k;
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       k, k)];
    for ( fn = l = 0;  fn < G1_DBDIM;  fn++ )
      for ( i = 0;  i <= fn;  i++, l++ ) {
        s = _g1h_Q2Integral0d ( &jac[k*3*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+fn)*4*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+i)*4*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[(k*3+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+fn)*4+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+1)*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[(k*3+2)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+fn)*4+3)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+3)*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[j*3*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+fn)*4+2)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+2)*G1_NQUAD] );
        Aij[l] += C*s;
      }
  }

  for ( k = 0; k < hole_k; k++ ) {
    j = (k+1) % hole_k;
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       hole_k-1, k)];
    for ( i = 0; i < G1_DBDIM; i++ )
      for ( fn = 0; fn < nfunc_a; fn++ ) {
        if ( k < hole_k-1 )
          l = (fn+G1_DBDIM)*G1_DBDIM+i;
        else
          l = pkn_SymMatIndex(G1_DBDIM+fn,i);
        s = _g1h_Q2Integral0d ( &jac[k*3*G1_NQUAD],
                                &lapj[(fn*hole_k+k)*3*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+i)*4*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[(k*3+1)*G1_NQUAD],
                                 &lapj[((fn*hole_k+k)*3+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+1)*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[(k*3+2)*G1_NQUAD],
                                 &lapj[((fn*hole_k+k)*3+2)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+3)*G1_NQUAD] );
        s += _g1h_Q2Integral0d ( &jac[j*3*G1_NQUAD],
                                 &lapj[((fn*hole_k+j)*3)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+2)*G1_NQUAD] );
        Aij[l] += C*s;
      }
  }

  for ( k = 0, j = hole_k-1;  k < hole_k;  j = k++ ) {
    if ( k == 0 )
      Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, j, k)];
    else
      Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, k, j)];
    for ( fn = 0; fn < G1_DBDIM; fn++ )
      for ( i = 0; i < G1_DBDIM; i++ ) {
        s = _g1h_Q2Integral0d (&jac[k*3*G1_NQUAD],
                               &lapjc[(k*G1_DBDIM+fn)*4*G1_NQUAD],
                               &lapjc[((j*G1_DBDIM+i)*4+2)*G1_NQUAD] );
        if ( k == 0 )
          l = i*G1_DBDIM+fn;
        else
          l = fn*G1_DBDIM+i;
        Aij[l] += C*s;
      }
  }
        /* right-hand side matrix */
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    supp = _g1h_ExtendSupport ( hole_k, support_b[fn] );
    for ( i = 0; i < hole_k; i++ ) {
      j = (i+1) % hole_k;
      if ( supp & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesd ( domain, fn, (i+hole_k-1) % hole_k,
                  &ec00, &ec01, &ec10, &ec11, &ed00, &ed01, &ed10, &ed11 );
        _g1h_GetBFBPatchCurvesd ( domain, fn, i,
                  &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        gh_GetDomSurrndBFuncd ( domain, fn, j, 1, eicp1 );
        gh_GetDomSurrndBFuncd ( domain, fn, i, 2, eicp2 );
        if ( !_g1h_Q2TabLaplacianJumpd ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                           atkn, ahfunc, adhfunc, addhfunc,
                           ec00, ec01, ec10, ec11, ed00, ed01, ed10, ed11,
                           &trd[(30*i+5)*G1_NQUAD],
                           fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                           &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                           &trd[(30*i+20)*G1_NQUAD],
                           eicp1, &trd[(30*i+15)*G1_NQUAD],
                           eicp2, &trd[(30*i+25)*G1_NQUAD],
                           &lapjb[3*i*G1_NQUAD],
                           &lapjb[(3*i+1)*G1_NQUAD],
                           &lapjb[(3*i+2)*G1_NQUAD] ) )
          goto failure;
      }
      else
        memset ( &lapjb[i*3*G1_NQUAD], 0, 3*G1_NQUAD*sizeof(double) );
    }
    for ( i = 0; i < hole_k; i++ )
      if ( supp & (0x0001 << i) ) {
        k = (i+1) % hole_k;
        for ( j = 0; j < G1_DBDIM; j++ ) {
          s = _g1h_Q2Integral0d ( &jac[i*3*G1_NQUAD],
                                  &lapjc[(i*G1_DBDIM+j)*4*G1_NQUAD],
                                  &lapjb[i*3*G1_NQUAD] );
          s += _g1h_Q2Integral0d ( &jac[(i*3+1)*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+1)*G1_NQUAD],
                                   &lapjb[(3*i+1)*G1_NQUAD] );
          s += _g1h_Q2Integral0d ( &jac[(i*3+2)*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+3)*G1_NQUAD],
                                   &lapjb[(3*i+2)*G1_NQUAD] );
          s += _g1h_Q2Integral0d ( &jac[k*3*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+2)*G1_NQUAD],
                                   &lapjb[k*3*G1_NQUAD] );
          bmat[(i*G1_DBDIM+j)*nfunc_b+fn] += C*s;
        }
      }
    for ( i = 0; i < nfunc_a; i++ )
      bmat[(i+nfunc_c)*nfunc_b+fn] += C*_g1h_Q2Integrald ( hole_k, G1_NQUAD, jac,
                                           0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                           supp, lapjb );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtComputeFormMatrixd*/

boolean g1h_Q2ExtDecomposeMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *lmat;
  int    hole_k, nfunc_a, size;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2EAMat )
    if ( !g1h_Q2ExtComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat ) {
    hole_k = domain->hole_k;
    nfunc_a = privateG1->nfunc_a;
    size = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
    lmat = privateG1->Q2ELMat = malloc ( size*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG1->Q2EAMat, size*sizeof(double) );
    if ( !pkn_Block3CholeskyDecompMd ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       lmat ) ) {
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtDecomposeMatrixd*/

boolean g1h_Q2ExtFillHoled ( GHoleDomaind *domain,
                             int spdimen, CONST_ double *hole_cp,
                             double *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const double *cp,
                                                void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  double *x, *fc00, *Bi, *Bk, *Lii;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2EAMat )
    if ( !g1h_Q2ExtComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat )
    if ( !g1h_Q2ExtDecomposeMatrixd ( domain ) )
      goto failure;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  Bi = privateG1->Q2EBMat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSided ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  Lii = privateG1->Q2ELMat;
  pkn_Block3LowerTrMSolved ( hole_k-1, G1_DBDIM, nfunc_a+G1_DBDIM, Lii,
                             spdimen, spdimen, x );
  pkn_Block3UpperTrMSolved ( hole_k-1, G1_DBDIM, nfunc_a+G1_DBDIM, Lii,
                             spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g1h_OutputExtPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtFillHoled*/

