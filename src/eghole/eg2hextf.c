
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

#undef CONST_
#define CONST_

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


boolean _g2h_GetExtBlockAddressesf ( GHoleDomainf *domain,
                                     float **Aii, float **Aki, float **Akk,
                                     float **Bi, float **Bk, float **Lii )
{
  G2HolePrivateRecf *privateG2;
  int  hole_k, nfunc_a, nfunc_b;

  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;

  if ( !privateG2->EAmat ) privateG2->EAmat =
        malloc ( pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a )*sizeof(float) );
  if ( !privateG2->EBmat ) privateG2->EBmat =
        malloc ( (hole_k*G2_DBDIM+nfunc_a)*nfunc_b*sizeof(float) );
  if ( !privateG2->ELmat ) privateG2->ELmat =
        malloc ( pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a )*sizeof(float) );
  if ( !privateG2->EAmat || !privateG2->EBmat || !privateG2->ELmat )
    return false;

  *Aii = privateG2->EAmat;
  *Aki = &privateG2->EAmat[pkn_Block1FindBlockPos( hole_k, G2_DBDIM, nfunc_a, hole_k, 0 )];
  *Akk = &privateG2->EAmat[pkn_Block1FindBlockPos( hole_k, G2_DBDIM, nfunc_a, hole_k, hole_k )];
  *Bi  = privateG2->EBmat;
  *Bk  = &(*Bi)[hole_k*G2_DBDIM*nfunc_b];
  *Lii = privateG2->ELmat;

  return true;
} /*_g2h_GetExtBlockAddressesf*/

static void _g2h_TabBezierPolynomialsDer3f ( int nkn, const float *kn,
                                             float *bezp, float *dbezp,
                                             float *ddbezp, float *dddbezp )
{
  int   i, j, k, b;
  float s, t, u, v, binom[G2H_FINALDEG-2], bpoly[G2H_FINALDEG-2];

        /* compute the binomial coefficients */
  binom[0] = 1.0;
  b = G2H_FINALDEG-3;
  for ( i = 1; i <= G2H_FINALDEG-3; i++ ) {
    binom[i] = (float)b;
    b = (b*(G2H_FINALDEG-3-i))/(i+1);
  }

  for ( j = k = 0;  j < nkn;  j++, k += G2H_FINALDEG-5 ) {
    t = kn[j];
    s = (float)(1.0-t);

        /* evaluate the polynomials of degree G2H_FINALDEG-3 at the knot t */
    memcpy ( bpoly, binom, (G2H_FINALDEG-2)*sizeof(float) );
    u = t;  v = s;
    bpoly[1] *= u;  bpoly[G2H_FINALDEG-4] *= v;
    for ( i = 2; i <= G2H_FINALDEG-3; i++ )
      { u *= t;  v *= s;  bpoly[i] *= u;  bpoly[G2H_FINALDEG-3-i] *= v; }

        /* evaluate the third order derivatives */
    b = G2H_FINALDEG*(G2H_FINALDEG-1)*(G2H_FINALDEG-2);
    for ( i = 0; i <= G2H_FINALDEG-6; i++ )
      dddbezp[k+i] =
        (float)(b*((bpoly[i]+3.0*bpoly[i+2])-(3.0*bpoly[i+1]+bpoly[i+3])));

        /* evaluate the polynomials of degree G2H_FINALDEG-2 */
        /* - de Casteljau algorithm */
    for ( i = 0; i <= G2H_FINALDEG-4; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];

        /* evaluate the second order derivatives */
    b = G2H_FINALDEG*(G2H_FINALDEG-1);
    for ( i = 0; i <= G2H_FINALDEG-6; i++ )
      ddbezp[k+i] = (float)(b*((bpoly[i]+bpoly[i+2])-2.0*bpoly[i+1]));

        /* evaluate the polynomials of degree G2H_FINALDEG-1 */
    for ( i = 0; i <= G2H_FINALDEG-5; i++ )
      bpoly[i] = s*bpoly[i+1] + t*bpoly[i];

        /* evaluate the first order derivatives */
    for ( i = 0; i <= G2H_FINALDEG-6; i++ )
      dbezp[k+i] = (float)G2H_FINALDEG*(bpoly[i]-bpoly[i+1]);

        /* evaluate the polynomials of degree G2H_FINALDEG */
    for ( i = 0; i <= G2H_FINALDEG-6; i++ )
      bezp[k+i] = s*bpoly[i+1] + t*bpoly[i];
  }
} /*_g2h_TabBezierPolynomialsDer3f*/

static void _g2h_TabTensBezf ( int nkn, const float *bez1, const float *bez2,
                               float *tbez )
{
  int i, j, k, l, m;

  for ( i = m = 0; i < G2H_FINALDEG-5; i++ )
    for ( j = 0; j < G2H_FINALDEG-5; j++ )
      for ( k = 0; k < nkn; k++ )
        for ( l = 0;  l < nkn;  l++, m++ )
          tbez[m] = bez1[k*(G2H_FINALDEG-5)+i]*bez2[l*(G2H_FINALDEG-5)+j];
} /*_g2h_TabTensBezf*/

boolean _g2h_TabTensBezPolyDer3f ( int nkn, const float *tkn,
             float *tbez, float *tbezu, float *tbezv,
             float *tbezuu, float *tbezuv, float *tbezvv,
             float *tbezuuu, float *tbezuuv, float *tbezuvv, float *tbezvvv )
{
  void  *sp;
  float *bezp, *dbezp, *ddbezp, *dddbezp;

  sp = pkv_GetScratchMemTop ();
  bezp = pkv_GetScratchMemf ( (G2H_FINALDEG-5)*nkn*4 );
  if ( !bezp )
    goto failure;
  dbezp = &bezp[(G2H_FINALDEG-5)*nkn];
  ddbezp = &dbezp[(G2H_FINALDEG-5)*nkn];
  dddbezp = &ddbezp[(G2H_FINALDEG-5)*nkn];

  _g2h_TabBezierPolynomialsDer3f ( G2_NQUAD, tkn, bezp, dbezp, ddbezp, dddbezp );
  if ( tbez )
    _g2h_TabTensBezf ( nkn, bezp, bezp, tbez );
  _g2h_TabTensBezf ( nkn, dbezp, bezp, tbezu );
  _g2h_TabTensBezf ( nkn, bezp, dbezp, tbezv );
  _g2h_TabTensBezf ( nkn, ddbezp, bezp, tbezuu );
  _g2h_TabTensBezf ( nkn, dbezp, dbezp, tbezuv );
  _g2h_TabTensBezf ( nkn, bezp, ddbezp, tbezvv );
  _g2h_TabTensBezf ( nkn, dddbezp, bezp, tbezuuu );
  _g2h_TabTensBezf ( nkn, ddbezp, dbezp, tbezuuv );
  _g2h_TabTensBezf ( nkn, dbezp, ddbezp, tbezuvv );
  _g2h_TabTensBezf ( nkn, bezp, dddbezp, tbezvvv );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabTensBezPolyDer3f*/

static void _g2h_TabLaplacianGrad00f ( int nkn,
                const float *tbezu, const float *tbezv,
                const float *tbezuu, const float *tbezuv, const float *tbezvv,
                const float *tbezuuu, const float *tbezuuv,
                const float *tbezuvv, const float *tbezvvv,
                const float *trd, vector2f *lapgrad )
{
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
} /*_g2h_TabLaplacianGrad00f*/

float _g2h_Integral0f ( int nknsq, const float *jac,
                        const vector2f *func1, const vector2f *func2 )
{
  int    i;
  double s;

  s = 0.0;
  for ( i = 0; i < nknsq; i++ )
    s += (func1[i].x*func2[i].x + func1[i].y*func2[i].y)*jac[i];
  s /= (double)(nknsq);
  return (float)s;
} /*_g2h_Integral0f*/

boolean g2h_DecomposeExtMatrixf ( GHoleDomainf *domain )
{
  void  *sp;
  float *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int   hole_k, nfunc_a;
  G2HolePrivateRecf *privateG2;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG2 ) {
    if ( !g2h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  if ( !privateG2->EAmat )
    if ( !g2h_ComputeExtFormMatrixf ( domain ) )
      goto failure;
  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  memcpy ( Lii, Aii, pkn_Block1ArraySize( hole_k, G2_DBDIM, nfunc_a )*sizeof(float) );
  if ( !pkn_Block1CholeskyDecompMf ( hole_k, G2_DBDIM, nfunc_a, Lii ) ) {
    domain->error_code = G2H_ERROR_NONPOSITIVE_EXT_MATRIX;
    goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_DecomposeExtMatrixf*/

boolean g2h_ComputeExtFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  int      i, j, k, l, n, m;
  float    *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  float    *aki;
  float    *tkn, *trd, *jac, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  float    *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *tbezuuu, *tbezuuv, *tbezuvv, *tbezvvv;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;
  vector2f *lgr, *lgrb, *lgrc;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  unsigned short    *support_b;


  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG2 ) {
    if ( !g2h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  support_b = privateG->support_b;
  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  n = hole_k*G2_NQUADSQ;
  m = G2_DBDIM*G2_NQUADSQ;

        /* allocate memory */
  tkn = pkv_GetScratchMemf ( 25*G2_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 18*n );
  lgr = (vector2f*)pkv_GetScratchMem (
                       ((nfunc_a+1)*n+nfunc_c*G2_NQUADSQ)*sizeof(vector2f) );
  tbezu = pkv_GetScratchMemf ( 9*m );
  if ( !tkn || !jac || !trd || !lgr || !tbezu )
    goto failure;
  hfunc = &tkn[G2_NQUAD];         dhfunc = &hfunc[6*G2_NQUAD];
  ddhfunc = &dhfunc[6*G2_NQUAD];  dddhfunc = &ddhfunc[6*G2_NQUAD];
  tbezv = &tbezu[m];           tbezuu = &tbezv[m];
  tbezuv = &tbezuu[m];         tbezvv = &tbezuv[m];
  tbezuuu = &tbezvv[m];        tbezuuv = &tbezuuu[m];
  tbezuvv = &tbezuuv[m];       tbezvvv = &tbezuvv[m];
  lgrc = &lgr[nfunc_a*n];
  lgrb = &lgrc[nfunc_c*G2_NQUADSQ];
  _gh_PrepareTabKnotsf ( G2_NQUAD, privateG2->opt_quad, tkn );

        /* prepare the evaluation of basis functions */
  if ( !mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, G2_NQUAD, tkn,
                                   hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;
  _g2h_TabTensBezPolyDer3f ( G2_NQUAD, tkn, NULL, tbezu, tbezv, tbezuu, tbezuv,
                             tbezvv, tbezuuu, tbezuuv, tbezuvv, tbezvvv );

  for ( i = 0; i < hole_k; i++ ) {
    _g2h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    _g2h_TabDiPatchJac3f ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                           c00, c01, c02, c10, c11, c12,
                           d00, d01, d02, d10, d11, d12,
                           &jac[i*G2_NQUADSQ], &trd[i*G2_NQUADSQ*18] );
  }

        /* evaluate the Laplacian gradient of the main basis functions */
  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      _g2h_TabLaplacianGrad0f ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                fc00, fc01, fc02, fd00, fd01, fd02,
                                &trd[i*G2_NQUADSQ*18],
                                &lgr[(j*hole_k+i)*G2_NQUADSQ] );
    }

        /* evaluate the Laplacian gradient of the extended basis functions */
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < G2_DBDIM; j++ )
      _g2h_TabLaplacianGrad00f ( G2_NQUAD, &tbezu[j*G2_NQUADSQ], &tbezv[j*G2_NQUADSQ],
            &tbezuu[j*G2_NQUADSQ], &tbezuv[j*G2_NQUADSQ], &tbezvv[j*G2_NQUADSQ],
            &tbezuuu[j*G2_NQUADSQ], &tbezuuv[j*G2_NQUADSQ], &tbezuvv[j*G2_NQUADSQ],
            &tbezvvv[j*G2_NQUADSQ], &trd[i*G2_NQUADSQ*18], &lgrc[(i*G2_DBDIM+j)*G2_NQUADSQ] );

        /* compute the form matrix coefficients */
          /* diagonal blocks first */
  for ( k = 0; k < hole_k; k++ )
    for ( i = l = 0;  i < G2_DBDIM;  i++ )
      for ( j = 0;  j <= i;  j++, l++ )
        Aii[k*G2_DIAGBLSIZE+l] = _g2h_Integral0f ( G2_NQUADSQ, &jac[k*G2_NQUADSQ],
                          &lgrc[k*m+i*G2_NQUADSQ], &lgrc[k*m+j*G2_NQUADSQ] );
  for ( i = l = 0;  i < nfunc_a;  i++ )
    for ( j = 0;  j <= i;  j++, l++ )
      Akk[l] = _g2h_Integralf ( hole_k, G2_NQUADSQ, jac,
                                0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );

          /* off-diagonal blocks */
  for ( k = 0; k < hole_k; k++ ) {
    aki = &Aki[k*G2_DBDIM*nfunc_a];
    for ( i = 0; i < G2_DBDIM; i++ )
      for ( j = 0; j < nfunc_a; j++ )
        aki[j*G2_DBDIM+i] = _g2h_Integral0f ( G2_NQUADSQ, &jac[k*G2_NQUADSQ],
                          &lgrc[k*m+i*G2_NQUADSQ], &lgr[j*n+k*G2_NQUADSQ] );
  }

        /* now the right-hand side matrix */
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        _g2h_GetBFBPatchCurvesf ( domain, j, k,
                                  &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                  &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
        _g2h_TabLaplacianGradf ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                 fc00, fc01, fc02, fc10, fc11, fc12,
                                 fd00, fd01, fd02, fd10, fd11, fd12,
                                 &trd[k*G2_NQUADSQ*18], &lgrb[k*G2_NQUADSQ] );
      }
          /* scalar products with the extended basis functions */
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        for ( i = 0; i < G2_DBDIM; i++ )
          Bi[(k*G2_DBDIM+i)*nfunc_b+j] = _g2h_Integral0f ( G2_NQUADSQ, &jac[k*G2_NQUADSQ],
                                      &lgrc[k*m+i*G2_NQUADSQ], &lgrb[k*G2_NQUADSQ] );
      }
      else
        for ( i = 0; i < G2_DBDIM; i++ )
          Bi[(k*G2_DBDIM+i)*nfunc_b+j] = 0.0;

          /* scalar products with the main basis functions */
    for ( i = 0; i < nfunc_a; i++ )
      Bk[i*nfunc_b+j] = _g2h_Integralf ( hole_k, G2_NQUADSQ, jac,
                          0xFFFF, &lgr[i*n], support_b[j], lgrb );
  }

  pkv_SetScratchMemTop ( sp );
  if ( !g2h_DecomposeExtMatrixf ( domain ) )
    return false;

  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ComputeExtFormMatrixf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2h_SetExtRightSidef ( GHoleDomainf *domain,
                                const float *Bi, const float *Bk,
                                int spdimen, CONST_ float *hole_cp,
                                float *fc00, float *b )
{
  void    *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  unsigned char *bfcpn;
  float  *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  float         *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float  *x, *y, *cp;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c;
  int    i, j, k, l;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  bfcpn = privateG->bfcpn;

  x = pkv_GetScratchMemf ( spdimen*nfunc_b );
  y = pkv_GetScratchMemf ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G2_CROSSDEGSUM+6)*2*hole_k*spdimen*sizeof(float) );
  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy ( &x[j*spdimen], cp, spdimen*sizeof(float) );
    for ( i = k = 0;  i < hole_k;  i++, k+= spdimen ) {
      _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,  
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fc10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fc11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS12DEG+1, 1, 1, bc12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fc12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fc12[k*(G2_CROSS12DEG+1)] );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],  
                       0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],  
                       0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],  
                       0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fd10[k*(G2_CROSS10DEG+1)],  
                       0, y, 0, &fd10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fd11[k*(G2_CROSS11DEG+1)],  
                       0, y, 0, &fd11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS12DEG+1, 1, 1, bd12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fd12[k*(G2_CROSS12DEG+1)],  
                       0, y, 0, &fd12[k*(G2_CROSS12DEG+1)] );
    }
  }

  if ( b ) {
    pkn_MultMatrixf ( nfunc_c, nfunc_b, nfunc_b, Bi,
                      spdimen, spdimen, x, spdimen, b );
    pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bk,
                      spdimen, spdimen, x, spdimen, &b[hole_k*G2_DBDIM*spdimen] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_SetExtRightSidef*/

boolean _g2h_OutputExtPatchesf ( GHoleDomainf *domain,
                    int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                    void (*outpatch) ( int n, int m, const float *cp,
                                       void *usrptr ) )
{
  void    *sp;
  G2HolePrivateRecf *privateG2;
  float  *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;
  float         *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  float  *y, *cp;
  int    hole_k, nfunc_a, nfunc_c, degu, degv;
  int    i, j, k, l, m;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;

  y = pkv_GetScratchMemf ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !y )
    goto failure;

  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[(nfunc_c+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );

      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
    }
  }  

  for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
    mbs_BezC2CoonsToBezf ( spdimen,
                   G2_CROSS00DEG, &fc00[k*(G2_CROSS00DEG+1)],
                   G2_CROSS01DEG, &fc01[k*(G2_CROSS01DEG+1)],
                   G2_CROSS02DEG, &fc02[k*(G2_CROSS02DEG+1)],
                   G2_CROSS10DEG, &fc10[k*(G2_CROSS10DEG+1)],
                   G2_CROSS11DEG, &fc11[k*(G2_CROSS11DEG+1)],
                   G2_CROSS12DEG, &fc12[k*(G2_CROSS12DEG+1)],
                   G2_CROSS00DEG, &fd00[k*(G2_CROSS00DEG+1)],
                   G2_CROSS01DEG, &fd01[k*(G2_CROSS01DEG+1)],
                   G2_CROSS02DEG, &fd02[k*(G2_CROSS02DEG+1)],
                   G2_CROSS10DEG, &fd10[k*(G2_CROSS10DEG+1)],
                   G2_CROSS11DEG, &fd11[k*(G2_CROSS11DEG+1)],
                   G2_CROSS12DEG, &fd12[k*(G2_CROSS12DEG+1)],
                   &degu, &degv, y );

    for ( l = j = 0;  l < G2H_FINALDEG-5;  l++ )
      for ( m = 0;  m < G2H_FINALDEG-5;  m++, j++ )
        pkn_SubtractMatrixf ( 1, spdimen, 0, &y[((l+3)*(G2H_FINALDEG+1)+(m+3))*spdimen],
                              0, &x[(i*G2_DBDIM+j)*spdimen],
                              0, &y[((l+3)*(G2H_FINALDEG+1)+(m+3))*spdimen] );

    outpatch ( degu, degv, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_OutputExtPatchesf*/

boolean g2h_ExtFillHolef ( GHoleDomainf *domain,
                           int spdimen, CONST_ float *hole_cp,
                           float *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) )
{
  void    *sp;
  G2HolePrivateRecf *privateG2;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  float  *fc00, *x;
  int    hole_k, nfunc_a, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->EAmat )
    if ( !g2h_DecomposeExtMatrixf ( domain ) )
      goto failure;

  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  if ( !_g2h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g2h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ExtFillHolef*/

boolean g2h_ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void    *sp;
  G2HolePrivateRecf *privateG2;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii, *cmat, *rcmat;
  float  *fc00, *b, *x, *y;
  int    hole_k, nfunc_a, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->ECmat || !privateG2->ERCmat ||
       privateG2->extnconstr != nconstr ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  cmat = privateG2->ECmat;
  rcmat = privateG2->ERCmat;
  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixf ( nconstr, nbf, nbf, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixf ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemf ( spdimen*nbf );
  if ( !y ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixf ( nconstr, nbf, nbf, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_Block1LowerTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, y );
  pkn_Block1UpperTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemf ( spdimen*nbf );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g2h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ExtFillHoleConstrf*/

boolean g2h_ExtFillHoleAltConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void    *sp;
  G2HolePrivateRecf *privateG2;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii, *acmat, *arcmat;
  float  *fc00, *b, *x, *y;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->AECmat || !privateG2->AERCmat ||
       privateG2->extnaconstr != naconstr || privateG2->extacdim != spdimen ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  acmat = privateG2->AECmat;
  arcmat = privateG2->AERCmat;
  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixf ( nbf, spdimen, spdimen, x, nbf, y );
  pkn_MultMatrixf ( naconstr, spdimen*nbf, spdimen*nbf, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixf ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixf ( naconstr, nbf*spdimen, nbf*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixf ( spdimen, nbf, nbf, y, spdimen, b );
  pkn_Block1LowerTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, b );
  pkn_Block1UpperTrMSolvef ( hole_k, G2_DBDIM, nfunc_a, Lii, spdimen, spdimen, b );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, b, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g2h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ExtFillHoleAltConstrf*/

