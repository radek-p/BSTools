
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/*#define DEBUG*/

/* ///////////////////////////////////////////////////////////////////////// */
static void _g1h_TabBezierPolynomialsDer3f ( int nkn, const float *kn,
                                             float *bezp, float *dbezp,
                                             float *ddbezp, float *dddbezp )
{
  int   i, j, k, b;
  float s, t, u, v, binom[G1H_FINALDEG-2], bpoly[G1H_FINALDEG];

        /* compute the binomial coefficients */
  binom[0] = 1.0;
  b = G1H_FINALDEG-3;
  for ( i = 1; i <= G1H_FINALDEG-3; i++ ) {
    binom[i] = (float)b;
    b = (b*(G1H_FINALDEG-3-i))/(i+1);
  }

  for ( j = k = 0;  j < nkn;  j++, k += G1H_FINALDEG-3 ) {
    t = kn[j];
    s = (float)(1.0-t);

        /* evaluate the polynomials of degree G1H_FINALDEG-3 at the knot t */
    bpoly[0] = bpoly[G1H_FINALDEG-1] = 0.0;
    memcpy ( &bpoly[1], binom, (G1H_FINALDEG-2)*sizeof(float) );
    u = t;  v = s;
    bpoly[2] *= u;  bpoly[G1H_FINALDEG-3] *= v;
    for ( i = 2; i <= G1H_FINALDEG-3; i++ )
      { u *= t;  v *= s;  bpoly[i+1] *= u;  bpoly[G1H_FINALDEG-2-i] *= v; }

        /* evaluate the third order derivatives */
    b = G1H_FINALDEG*(G1H_FINALDEG-1)*(G1H_FINALDEG-2);
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      dddbezp[k+i] =
        (float)(b*((bpoly[i]+3.0*bpoly[i+2])-(3.0*bpoly[i+1]+bpoly[i+3])));

        /* evaluate the polynomials of degree G1H_FINALDEG-2 */
        /* - de Casteljau algorithm */
    bpoly[0] = s*bpoly[1];
    for ( i = 1; i < G1H_FINALDEG-2; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];
    bpoly[G1H_FINALDEG-2] = t*bpoly[G1H_FINALDEG-2];

        /* evaluate the second order derivatives */
    b = G1H_FINALDEG*(G1H_FINALDEG-1);
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      ddbezp[k+i] = (float)(b*((bpoly[i]+bpoly[i+2])-2.0*bpoly[i+1]));

        /* evaluate the polynomials of degree G1H_FINALDEG-1 */
    for ( i = 0; i <= G1H_FINALDEG-3; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];

        /* evaluate the first order derivatives */
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      dbezp[k+i] = (float)G1H_FINALDEG*(bpoly[i]-bpoly[i+1]);

        /* evaluate the polynomials of degree G1H_FINALDEG */
    for ( i = 0; i <= G1H_FINALDEG-4; i++ )
      bezp[k+i] = s*bpoly[i+1] + t*bpoly[i];
  }
} /*_g1h_TabBezierPolynomialsDer3f*/

static void _g1h_TabTensBezf ( int nkn, const float *bez1, const float *bez2,
                               float *tbez )
{
  int i, j, k, l, m;

/* !!! this procedure is identical to that in eg1hextf.c !!! */
  for ( i = m = 0; i < G1H_FINALDEG-3; i++ )
    for ( j = 0; j < G1H_FINALDEG-3; j++ )  
      for ( k = 0; k < nkn; k++ )
        for ( l = 0;  l < nkn;  l++, m++ )
          tbez[m] = bez1[k*(G1H_FINALDEG-3)+i]*bez2[l*(G1H_FINALDEG-3)+j];
} /*_g1h_TabTensBezf*/

boolean _g1h_TabTensBezPolyDer3f ( int nkn, const float *tkn,
             float *tbez, float *tbezu, float *tbezv,
             float *tbezuu, float *tbezuv, float *tbezvv,
             float *tbezuuu, float *tbezuuv, float *tbezuvv, float *tbezvvv )
{
  void  *sp;
  float *bezp, *dbezp, *ddbezp, *dddbezp;

  sp = pkv_GetScratchMemTop ();
  bezp = pkv_GetScratchMemf ( (G1H_FINALDEG-3)*nkn*4 );
  if ( !bezp )
    goto failure;
  dbezp = &bezp[(G1H_FINALDEG-3)*nkn];
  ddbezp = &dbezp[(G1H_FINALDEG-3)*nkn];
  dddbezp = &ddbezp[(G1H_FINALDEG-3)*nkn];

  _g1h_TabBezierPolynomialsDer3f ( G1_NQUAD, tkn, bezp, dbezp, ddbezp, dddbezp );
  if ( tbez )
    _g1h_TabTensBezf ( nkn, bezp, bezp, tbez );
  _g1h_TabTensBezf ( nkn, dbezp, bezp, tbezu );
  _g1h_TabTensBezf ( nkn, bezp, dbezp, tbezv );
  _g1h_TabTensBezf ( nkn, ddbezp, bezp, tbezuu );
  _g1h_TabTensBezf ( nkn, dbezp, dbezp, tbezuv );
  _g1h_TabTensBezf ( nkn, bezp, ddbezp, tbezvv );
  _g1h_TabTensBezf ( nkn, dddbezp, bezp, tbezuuu );
  _g1h_TabTensBezf ( nkn, ddbezp, dbezp, tbezuuv );
  _g1h_TabTensBezf ( nkn, dbezp, ddbezp, tbezuvv );
  _g1h_TabTensBezf ( nkn, bezp, dddbezp, tbezvvv );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabTensBezPolyDer3f*/

static void _g1h_TabLaplacianGrad00f ( int nkn,
                const float *tbezu, const float *tbezv,
                const float *tbezuu, const float *tbezuv, const float *tbezvv,
                const float *tbezuuu, const float *tbezuuv,
                const float *tbezuvv, const float *tbezvvv,
                const float *trd, vector2f *lapgrad )
{
/* !!! the same as _g2h_TabLaplacianGrad00f in the file eg2hextf.c !!! */
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
} /*_g1h_TabLaplacianGrad00f*/

boolean _g1h_TabLaplacianJump00f ( int nkn, const float *tkn, int fni,
                    const float *trdc00, const float *trdc10,
                    const float *trdd00, const float *trdd10,
                    float *lapc00, float *lapc10, float *lapd00, float *lapd10 )
{
  void  *sp;
  float *b, p, pu, pv, puu, puv, pvv, lap;
  int   fi, fj, i, j;

  sp = pkv_GetScratchMemTop ();
  b = pkv_GetScratchMemf ( (G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !b )
    goto failure;
  memset ( b, 0, (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(float) );
  fi = fni / (G1H_FINALDEG-3) + 2;
  fj = fni % (G1H_FINALDEG-3) + 2;
  b[fi*(G1H_FINALDEG+1)+fj] = 1.0;

/* !!! this is not as it ought to be!, to be replaced later !!! */
  if ( fi == 2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 0.0, tkn[i],
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdd00[j]*pu + trdd00[j+1]*pv +
            trdd00[j+2]*puu + trdd00[j+3]*puv + trdd00[j+4]*pvv;
      lapd00[i] = -lap;
    }
  }
  else
    memset ( lapd00, 0, nkn*sizeof(float) );
  if ( fi == G1H_FINALDEG-2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 1.0, tkn[i],
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdd10[j]*pu + trdd10[j+1]*pv +
            trdd10[j+2]*puu + trdd10[j+3]*puv + trdd10[j+4]*pvv;
      lapd10[i] = -lap;
    }
  }
  else
    memset ( lapd10, 0, nkn*sizeof(float) );

  if ( fj == 2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[i], 0.0,
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdc00[j]*pu + trdc00[j+1]*pv +
            trdc00[j+2]*puu + trdc00[j+3]*puv + trdc00[j+4]*pvv;
      lapc00[i] = lap;
    }
  }
  else
    memset ( lapc00, 0, nkn*sizeof(float) );
  if ( fj == G1H_FINALDEG-2 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[i], 1.0,
                                 &p, &pu, &pv, &puu, &puv, &pvv ) )
        goto failure;
      lap = trdc10[j]*pu + trdc10[j+1]*pv +
            trdc10[j+2]*puu + trdc10[j+3]*puv + trdc10[j+4]*pvv;
      lapc10[i] = -lap;
    }
  }
  else
    memset ( lapc10, 0, nkn*sizeof(float) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabLaplacianJump00f*/

static float _g1h_Q2Integral0f ( float *jac, float *lapj1, float *lapj2 )
{
  double s;
  int    i;

  for ( i = 0, s = 0.0;  i < G1_NQUAD; i++ )
    s += jac[i]*lapj1[i]*lapj2[i];
  s /= (double)G1_NQUAD;
  return (float)s;
} /*_g1h_Q2Integral0f*/

boolean g1h_Q2ExtComputeFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf    *privateG;
  G1HolePrivateRecf   *privateG1;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;
  int      i, j, k, l, n, m, fn;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  float    *atkn, *ahfunc, *adhfunc, *addhfunc;
  float    *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *tbezuuu, *tbezuuv, *tbezuvv, *tbezvvv;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *ec00, *ec01, *ec10, *ec11, *ed00, *ed01, *ed10, *ed11;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float    *trd, *lapj, *lapjb, *lapjc, *jac, *amat, *bmat, *Aij;
  vector2f *lgr, *lgrb, *lgrc;
  unsigned short *support_b, supp;
  int      option, ndata, *idata, Asize, Bsize;
  float    *fdata, *eicp1, *eicp2, C;
  float    s;
  point2f  *sicp;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisf ( domain ) )
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
    amat = privateG1->Q2EAMat = malloc ( Asize*sizeof(float) );
  if ( !(bmat = privateG1->Q2EBMat) )
    bmat = privateG1->Q2EBMat = malloc ( Bsize*sizeof(float) );
  if ( !amat || !bmat )
    goto failure;
  if ( !privateG->diam )
    privateG->diam = gh_DomainDiamf ( domain );

  n = hole_k*G1_NQUADSQ;
  m = G1_DBDIM*G1_NQUADSQ;

  tkn = pkv_GetScratchMemf ( 17*G1_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 18*n );
  lgr = (vector2f*)pkv_GetScratchMem (
                       ((nfunc_a+1)*n+nfunc_c*G1_NQUADSQ)*sizeof(vector2f) );
  tbezu = pkv_GetScratchMemf ( 9*m );
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

  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3f ( 0.0, 1.0, G1_NQUAD, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;
  _g1h_TabTensBezPolyDer3f ( G1_NQUAD, tkn, NULL, tbezu, tbezv, tbezuu, tbezuv,
                             tbezvv, tbezuuu, tbezuuv, tbezuvv, tbezvvv );

  memset ( amat, 0, Asize*sizeof(float) );
  memset ( bmat, 0, Bsize*sizeof(float) );

        /* integrate the Laplacian gradients */
  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_Q2TabDiPatchJac3f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_Q2TabLaplacianGrad0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                  fc00, fc01, fd00, fd01,
                                  &trd[i*G1_NQUADSQ*18],
                                  &lgr[(fn*hole_k+i)*G1_NQUADSQ] );
    }
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < G1_DBDIM; j++ )
      _g1h_TabLaplacianGrad00f ( G1_NQUAD, &tbezu[j*G1_NQUADSQ], &tbezv[j*G1_NQUADSQ],
            &tbezuu[j*G1_NQUADSQ], &tbezuv[j*G1_NQUADSQ], &tbezvv[j*G1_NQUADSQ],
            &tbezuuu[j*G1_NQUADSQ], &tbezuuv[j*G1_NQUADSQ], &tbezuvv[j*G1_NQUADSQ],
            &tbezvvv[j*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18], &lgrc[(i*G1_DBDIM+j)*G1_NQUADSQ] );

  Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                     hole_k-1, hole_k-1)];
  for ( i = 0; i < nfunc_a; i++ )
    for ( j = 0; j <= i; j++ )
      Aij[pkn_SymMatIndex(G1_DBDIM+i,G1_DBDIM+j)] =
        _g2h_Integralf ( hole_k, G1_NQUADSQ, jac,
                         0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  for ( k = 0; k < hole_k; k++ ) {
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       k, k)];
    for ( i = l = 0;  i < G1_DBDIM;  i++ )
      for ( j = 0;  j <= i;  j++, l++ )
        Aij[l] = _g2h_Integral0f ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
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
        Aij[l] = _g2h_Integral0f ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
                          &lgrc[k*m+i*G1_NQUADSQ], &lgr[j*n+k*G1_NQUADSQ] );
      }
  }

        /* right-hand side matrix */
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, k,
                      &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        _g1h_Q2TabLaplacianGradf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                   fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                                   &trd[k*G1_NQUADSQ*18], &lgrb[k*G1_NQUADSQ] );
      }

    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        for ( i = 0; i < G1_DBDIM; i++ )
          bmat[(k*G1_DBDIM+i)*nfunc_b+j] = _g2h_Integral0f ( G1_NQUADSQ, &jac[k*G1_NQUADSQ],
                                     &lgrc[k*m+i*G1_NQUADSQ], &lgrb[k*G1_NQUADSQ] );
      }

    for ( i = 0; i < nfunc_a; i++ )
      bmat[(k*G1_DBDIM+i)*nfunc_b+j] = _g2h_Integralf ( hole_k, G1_NQUADSQ, jac,
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
  C *= (float)(5.0/privateG->diam);  /* 3.0 = sqrt(G1_DBDIM)+1 */

        /* the arrays allocated for the previous stage are reused */
        /* their length is sufficient */
  lapj = &lgr[0].x;
  lapjb = &lapj[3*G1_NQUAD*hole_k*nfunc_a];
  lapjc = &lapjb[3*G1_NQUAD*hole_k];
        /* two additional arrays are however needed */
  atkn = pkv_GetScratchMemf ( 26 );
  sicp = (point2f*)pkv_GetScratchMem ( 16*sizeof(point2f) );
  if ( !atkn || !sicp )
    goto failure;
  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  eicp1 = (float*)sicp;  eicp2 = &eicp1[16];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  if ( !mbs_TabCubicHFuncDer2f ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
          /* compute the curve Jacobians */
    _g1h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS00DEG, c00, G1_NQUAD, tkn,
                                   &jac[3*i*G1_NQUAD] ) )
      goto failure;
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS10DEG, c10, G1_NQUAD, tkn,
                                   &jac[(3*i+1)*G1_NQUAD] ) )
      goto failure;
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS10DEG, d10, G1_NQUAD, tkn,
                                   &jac[(3*i+2)*G1_NQUAD] ) )
      goto failure;
          /* compute the coefficients for computing the Laplacians */
    if ( !_g1h_TabCurveLapCoeff0f ( c00, c01, c10, c11, d00, d01, d10, d11,
                            G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc,   
                            &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                            &trd[(30*j+5)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchf ( domain, j, 1, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1f ( sicp, G1_NQUAD, tkn, &trd[(30*i+15)*G1_NQUAD] );
    if ( !_gh_FindDomSurrndPatchf ( domain, i, 2, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1f ( sicp, G1_NQUAD, tkn, &trd[(30*i+25)*G1_NQUAD] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ ) {
    _g1h_GetBFAPatchCurvesf ( domain, fn, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1h_Q2TabLaplacianJump0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
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
      if ( !_g1h_TabLaplacianJump00f ( G1_NQUAD, tkn, i,
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
      Aij[l+j] += C*_g1h_Q2Integralf ( hole_k, G1_NQUAD, jac,
                                       0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                       0xFFFF, &lapj[j*hole_k*3*G1_NQUAD] );
  }

  for ( k = 0; k < hole_k; k++ ) {
    j = (k+1) % hole_k;
    Aij = &amat[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                       k, k)];
    for ( fn = l = 0;  fn < G1_DBDIM;  fn++ )
      for ( i = 0;  i <= fn;  i++, l++ ) {
        s = _g1h_Q2Integral0f ( &jac[k*3*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+fn)*4*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+i)*4*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[(k*3+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+fn)*4+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+1)*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[(k*3+2)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+fn)*4+3)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+3)*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[j*3*G1_NQUAD],
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
        s = _g1h_Q2Integral0f ( &jac[k*3*G1_NQUAD],
                                &lapj[(fn*hole_k+k)*3*G1_NQUAD],
                                &lapjc[(k*G1_DBDIM+i)*4*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[(k*3+1)*G1_NQUAD],
                                 &lapj[((fn*hole_k+k)*3+1)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+1)*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[(k*3+2)*G1_NQUAD],
                                 &lapj[((fn*hole_k+k)*3+2)*G1_NQUAD],
                                 &lapjc[((k*G1_DBDIM+i)*4+3)*G1_NQUAD] );
        s += _g1h_Q2Integral0f ( &jac[j*3*G1_NQUAD],
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
        s = _g1h_Q2Integral0f (&jac[k*3*G1_NQUAD],
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
        _g1h_GetBFBPatchCurvesf ( domain, fn, (i+hole_k-1) % hole_k,
                  &ec00, &ec01, &ec10, &ec11, &ed00, &ed01, &ed10, &ed11 );
        _g1h_GetBFBPatchCurvesf ( domain, fn, i,
                  &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        gh_GetDomSurrndBFuncf ( domain, fn, j, 1, eicp1 );
        gh_GetDomSurrndBFuncf ( domain, fn, i, 2, eicp2 );
        if ( !_g1h_Q2TabLaplacianJumpf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
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
        memset ( &lapjb[i*3*G1_NQUAD], 0, 3*G1_NQUAD*sizeof(float) );
    }
    for ( i = 0; i < hole_k; i++ )
      if ( supp & (0x0001 << i) ) {
        k = (i+1) % hole_k;
        for ( j = 0; j < G1_DBDIM; j++ ) {
          s = _g1h_Q2Integral0f ( &jac[i*3*G1_NQUAD],
                                  &lapjc[(i*G1_DBDIM+j)*4*G1_NQUAD],
                                  &lapjb[i*3*G1_NQUAD] );
          s += _g1h_Q2Integral0f ( &jac[(i*3+1)*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+1)*G1_NQUAD],
                                   &lapjb[(3*i+1)*G1_NQUAD] );
          s += _g1h_Q2Integral0f ( &jac[(i*3+2)*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+3)*G1_NQUAD],
                                   &lapjb[(3*i+2)*G1_NQUAD] );
          s += _g1h_Q2Integral0f ( &jac[k*3*G1_NQUAD],
                                   &lapjc[((i*G1_DBDIM+j)*4+2)*G1_NQUAD],
                                   &lapjb[k*3*G1_NQUAD] );
          bmat[(i*G1_DBDIM+j)*nfunc_b+fn] += C*s;
        }
      }
    for ( i = 0; i < nfunc_a; i++ )
      bmat[(i+nfunc_c)*nfunc_b+fn] += C*_g1h_Q2Integralf ( hole_k, G1_NQUAD, jac,
                                           0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                           supp, lapjb );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtComputeFormMatrixf*/

boolean g1h_Q2ExtDecomposeMatrixf ( GHoleDomainf *domain )
{
  void *sp;
  G1HolePrivateRecf *privateG1;
  float *lmat;
  int   hole_k, nfunc_a, size;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2EAMat )
    if ( !g1h_Q2ExtComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat ) {
    hole_k = domain->hole_k;
    nfunc_a = privateG1->nfunc_a;
    size = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
    lmat = privateG1->Q2ELMat = malloc ( size*sizeof(float) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG1->Q2EAMat, size*sizeof(float) );
    if ( !pkn_Block3CholeskyDecompMf ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
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
} /*g1h_Q2ExtDecomposeMatrixf*/

boolean g1h_Q2ExtFillHolef ( GHoleDomainf *domain,
                             int spdimen, CONST_ float *hole_cp,
                             float *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const float *cp,
                                                void *usrptr ) )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  float *x, *fc00, *Bi, *Bk, *Lii;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2EAMat )
    if ( !g1h_Q2ExtComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat )
    if ( !g1h_Q2ExtDecomposeMatrixf ( domain ) )
      goto failure;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  Bi = privateG1->Q2EBMat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  Lii = privateG1->Q2ELMat;
  pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, nfunc_a+G1_DBDIM, Lii,
                             spdimen, spdimen, x );
  pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, nfunc_a+G1_DBDIM, Lii,
                             spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtFillHolef*/

