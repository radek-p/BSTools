
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"


boolean _g1h_GetExtBlockAddressesf ( GHoleDomainf *domain,
                                     float **Aii, float **Aki, float **Akk,
                                     float **Bi, float **Bk, float **Lii )
{
  G1HolePrivateRecf *privateG1;
  int  hole_k, nfunc_a, nfunc_b;

  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;

  if ( !privateG1->EAmat ) privateG1->EAmat =
        malloc ( pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a )*sizeof(float) );
  if ( !privateG1->EBmat ) privateG1->EBmat =
        malloc ( (hole_k*G1_DBDIM+nfunc_a)*nfunc_b*sizeof(float) );
  if ( !privateG1->ELmat ) privateG1->ELmat =
        malloc ( pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a )*sizeof(float) );
  if ( !privateG1->EAmat || !privateG1->EBmat || !privateG1->ELmat )
    return false;

  *Aii = privateG1->EAmat;
  *Aki = &privateG1->EAmat[pkn_Block1FindBlockPos( hole_k, G1_DBDIM, nfunc_a, hole_k, 0)];
  *Akk = &privateG1->EAmat[pkn_Block1FindBlockPos( hole_k, G1_DBDIM, nfunc_a, hole_k, hole_k)];
  *Bi  = privateG1->EBmat;
  *Bk  = &privateG1->EBmat[hole_k*G1_DBDIM*nfunc_b];
  *Lii = privateG1->ELmat;

  return true;
} /*_g1h_GetExtBlockAddressesf*/

static void _g1h_TabBezierPolynomialsDer2f ( int nkn, const float *kn,
                                             float *bezp, float *dbezp,
                                             float *ddbezp )
{
  int   i, j, k, b;
  float s, t, u, v, binom[G1H_FINALDEG-1], bpoly[G1H_FINALDEG-1];

        /* compute the binomial coefficients */
  binom[0] = 1.0;
  b = G1H_FINALDEG-2;
  for ( i = 1; i <= G1H_FINALDEG-2; i++ ) {
    binom[i] = (float)b;
    b = (b*(G1H_FINALDEG-2-i))/(i+1);
  }

  for ( j = k = 0;  j < nkn;  j++, k += G1H_FINALDEG-3 ) {
    t = kn[j];
    s = (float)(1.0-t);

        /* evaluate the polynomials of degree G1H_FINALDEG-2 at the knot t */
    memcpy ( bpoly, binom, (G1H_FINALDEG-1)*sizeof(float) );
    u = t;  v = s;
    bpoly[1] *= u;  bpoly[G1H_FINALDEG-3] *= v;
    for ( i = 2; i <= G1H_FINALDEG-2; i++ )
      { u *= t;  v *= s;  bpoly[i] *= u;  bpoly[G1H_FINALDEG-2-i] *= v; }

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
} /*_g1h_TabBezierPolynomialsDer2f*/

static void _g1h_TabTensBezf ( int nkn, const float *bez1, const float *bez2,
                               float *tbez )
{
  int i, j, k, l, m;

  for ( i = m = 0; i < G1H_FINALDEG-3; i++ )
    for ( j = 0; j < G1H_FINALDEG-3; j++ )
      for ( k = 0; k < nkn; k++ )
        for ( l = 0;  l < nkn;  l++, m++ )
          tbez[m] = bez1[k*(G1H_FINALDEG-3)+i]*bez2[l*(G1H_FINALDEG-3)+j];
} /*_g1h_TabTensBezf*/

boolean _g1h_TabTensBezPolyDer2f ( int nkn, const float *tkn,
                     float *tbez, float *tbezu, float *tbezv,
                     float *tbezuu, float *tbezuv, float *tbezvv )
{
  void  *sp;
  float *bezp, *dbezp, *ddbezp;

  sp = pkv_GetScratchMemTop ();
  bezp = pkv_GetScratchMemf ( (G1H_FINALDEG-3)*nkn*3 );
  if ( !bezp )
    goto failure;
  dbezp = &bezp[(G1H_FINALDEG-3)*nkn];
  ddbezp = &dbezp[(G1H_FINALDEG-3)*nkn];

  _g1h_TabBezierPolynomialsDer2f ( nkn, tkn, bezp, dbezp, ddbezp );
  if ( tbez )
    _g1h_TabTensBezf ( nkn, bezp, bezp, tbez );
  _g1h_TabTensBezf ( nkn, dbezp, bezp, tbezu );
  _g1h_TabTensBezf ( nkn, bezp, dbezp, tbezv );
  _g1h_TabTensBezf ( nkn, ddbezp, bezp, tbezuu );
  _g1h_TabTensBezf ( nkn, dbezp, dbezp, tbezuv );
  _g1h_TabTensBezf ( nkn, bezp, ddbezp, tbezvv );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabTensBezPolyDer2f*/

static void _g1h_TabLaplacian00f ( int nkn,
                const float *tbezu, const float *tbezv,
                const float *tbezuu, const float *tbezuv, const float *tbezvv,
                const float *trd, float *lap )
{
  int i, k;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 5 )
    lap[i] = trd[k+0]*tbezu[i] + trd[k+1]*tbezv[i] + trd[k+2]*tbezuu[i] +
             trd[k+3]*tbezuv[i] + trd[k+4]*tbezvv[i];
} /*_g1h_TabLaplacian00f*/

static float _g1h_Integral0f ( const float *jac,
                               const float *func1, const float *func2 )
{
  int    i;
  double s;

  s = 0.0;
  for ( i = 0; i < G1_NQUADSQ; i++ )
    s += func1[i]*func2[i]*jac[i];
  s /= (double)(G1_NQUADSQ);
  return (float)s;
} /*_g1h_Integral0f*/

boolean g1h_DecomposeExtMatrixf ( GHoleDomainf *domain )
{
  void  *sp;
  float *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int   hole_k, nfunc_a;
  G1HolePrivateRecf *privateG1;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  if ( !privateG1->EAmat )
    if ( !g1h_ComputeExtFormMatrixf ( domain ) )
      goto failure;
  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  memcpy ( Lii, Aii, pkn_Block1ArraySize( hole_k, G1_DBDIM, nfunc_a )*sizeof(float) );
  if ( !pkn_Block1CholeskyDecompMf ( hole_k, G1_DBDIM, nfunc_a, Lii ) ) {
    domain->error_code = G1H_ERROR_NONPOSITIVE_EXT_MATRIX;
    goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_DecomposeExtMatrixf*/

boolean g1h_ComputeExtFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  int      i, j, k, l, n, m;
  float    *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  float    *aki;
  float    *tkn, *trd, *jac, *hfunc, *dhfunc, *ddhfunc;
  float    *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;
  float    *lgr, *lgrb, *lgrc;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  unsigned short    *support_b;


  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  support_b = privateG->support_b;
  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  n = hole_k*G1_NQUADSQ;
  m = G1_DBDIM*G1_NQUADSQ;

        /* allocate memory */
  tkn = pkv_GetScratchMemf ( 25*G1_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 5*n );
  lgr = pkv_GetScratchMemf ( (nfunc_a+1)*n+nfunc_c*G1_NQUADSQ );
  tbezu = pkv_GetScratchMemf ( 5*m );
  if ( !tkn || !jac || !trd || !lgr || !tbezu )
    goto failure;
  hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[4*G1_NQUAD];
  ddhfunc = &dhfunc[4*G1_NQUAD];
  tbezv = &tbezu[m];           tbezuu = &tbezv[m];
  tbezuv = &tbezuu[m];         tbezvv = &tbezuv[m];
  lgrc = &lgr[nfunc_a*n];
  lgrb = &lgrc[nfunc_c*G1_NQUADSQ];
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );

        /* prepare the evaluation of basis functions */
  if ( !mbs_TabCubicHFuncDer2f ( 0.0, 1.0, G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc ) )
    goto failure;
  _g1h_TabTensBezPolyDer2f ( G1_NQUAD, tkn, NULL, tbezu, tbezv,
                             tbezuu, tbezuv, tbezvv );

  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesf ( domain, i, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    _g1h_TabDiPatchJac2f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                           c00, c01, c10, c11, d00, d01, d10, d11,
                           &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*5] );
  }

        /* evaluate the Laplacian of the main basis functions */
  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, j, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabLaplacian0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            fc00, fc01, fd00, fd01, &trd[i*G1_NQUADSQ*5],
                            &lgr[(j*hole_k+i)*G1_NQUADSQ] );
    }

        /* evaluate the Laplacian of the extended basis functions */
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < G1_DBDIM; j++ )
      _g1h_TabLaplacian00f ( G1_NQUAD, &tbezu[j*G1_NQUADSQ], &tbezv[j*G1_NQUADSQ],
            &tbezuu[j*G1_NQUADSQ], &tbezuv[j*G1_NQUADSQ], &tbezvv[j*G1_NQUADSQ],
            &trd[i*G1_NQUADSQ*5], &lgrc[(i*G1_DBDIM+j)*G1_NQUADSQ] );

        /* compute the form matrix coefficients */
          /* diagonal blocks first */
  for ( k = 0; k < hole_k; k++ )
    for ( i = l = 0;  i < G1_DBDIM;  i++ )
      for ( j = 0;  j <= i;  j++, l++ )
        Aii[k*G1_DIAGBLSIZE+l] = _g1h_Integral0f ( &jac[k*G1_NQUADSQ],
                          &lgrc[k*m+i*G1_NQUADSQ], &lgrc[k*m+j*G1_NQUADSQ] );
  for ( i = l = 0;  i < nfunc_a;  i++ )
    for ( j = 0;  j <= i;  j++, l++ )
      Akk[l] = _g1h_Integralf ( hole_k, G1_NQUADSQ, jac,
                                0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );

          /* off-diagonal blocks */
  for ( k = 0; k < hole_k; k++ ) {
    aki = &Aki[k*G1_DBDIM*nfunc_a];
    for ( i = 0; i < G1_DBDIM; i++ )
      for ( j = 0; j < nfunc_a; j++ )
        aki[j*G1_DBDIM+i] = _g1h_Integral0f ( &jac[k*G1_NQUADSQ],
                          &lgrc[k*m+i*G1_NQUADSQ], &lgr[j*n+k*G1_NQUADSQ] );
  }

        /* now the right-hand side matrix */
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, k, &fc00, &fc01, &fc10, &fc11,
                                  &fd00, &fd01, &fd10, &fd11 );
        _g1h_TabLaplacianf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                             fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                             &trd[k*G1_NQUADSQ*5], &lgrb[k*G1_NQUADSQ] );
      }
          /* scalar products with the extended basis functions */
    for ( k = 0; k < hole_k; k++ )
      if ( support_b[j] & (0x0001 << k) ) {
        for ( i = 0; i < G1_DBDIM; i++ )
          Bi[(k*G1_DBDIM+i)*nfunc_b+j] = _g1h_Integral0f ( &jac[k*G1_NQUADSQ],
                                      &lgrc[k*m+i*G1_NQUADSQ], &lgrb[k*G1_NQUADSQ] );
      }
      else
        for ( i = 0; i < G1_DBDIM; i++ )
          Bi[(k*G1_DBDIM+i)*nfunc_b+j] = 0.0;

          /* scalar products with the main basis functions */
    for ( i = 0; i < nfunc_a; i++ )
      Bk[i*nfunc_b+j] = _g1h_Integralf ( hole_k, G1_NQUADSQ, jac,
                          0xFFFF, &lgr[i*n], support_b[j], lgrb );
  }

  pkv_SetScratchMemTop ( sp );
  if ( !g1h_DecomposeExtMatrixf ( domain ) )
    return false;

  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeExtFormMatrixf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_SetExtRightSidef ( GHoleDomainf *domain,
                                const float *Bi, const float *Bk,
                                int spdimen, CONST_ float *hole_cp,
                                float *fc00, float *b )
{
  void    *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  unsigned char *bfcpn;
  float  *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  float         *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float  *x, *y, *cp;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c;
  int    i, j, k, l;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  bfcpn = privateG->bfcpn;

  x = pkv_GetScratchMemf ( spdimen*nfunc_b );
  y = pkv_GetScratchMemf ( spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G1_CROSSDEGSUM+6)*2*hole_k*spdimen*sizeof(float) );
  G1GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy ( &x[j*spdimen], cp, spdimen*sizeof(float) );
    for ( i = k = 0;  i < hole_k;  i++, k+= spdimen ) {
      _g1h_GetBFBPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc10, &bc11,  
                                &bd00, &bd01, &bd10, &bd11 );

      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fc10[k*(G1_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fc11[k*(G1_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G1_CROSS11DEG+1)] );

      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],  
                       0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],  
                       0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fd10[k*(G1_CROSS10DEG+1)],  
                       0, y, 0, &fd10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fd11[k*(G1_CROSS11DEG+1)],  
                       0, y, 0, &fd11[k*(G1_CROSS11DEG+1)] );
    }
  }

  if ( b ) {
    pkn_MultMatrixf ( nfunc_c, nfunc_b, nfunc_b, Bi,
                      spdimen, spdimen, x, spdimen, b );
    pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bk,
                      spdimen, spdimen, x, spdimen, &b[hole_k*G1_DBDIM*spdimen] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_SetExtRightSidef*/

boolean _g1h_OutputExtPatchesf ( GHoleDomainf *domain, int spdimen,
                           CONST_ float *x, float *fc00, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  float  *bc00, *bc01, *bd00, *bd01;
  float  *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float  *y, *cp;
  int    hole_k, nfunc_a, nfunc_c, degu, degv;
  int    i, j, k, l, m;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;

  y = pkv_GetScratchMemf ( spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !y )
    goto failure;

  G1GetFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[(nfunc_c+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesf ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );

      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );

      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixf ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
    }
  }  

  for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
    mbs_BezC1CoonsToBezf ( spdimen,
                   G1_CROSS00DEG, &fc00[k*(G1_CROSS00DEG+1)],
                   G1_CROSS01DEG, &fc01[k*(G1_CROSS01DEG+1)],
                   G1_CROSS10DEG, &fc10[k*(G1_CROSS10DEG+1)],
                   G1_CROSS11DEG, &fc11[k*(G1_CROSS11DEG+1)],
                   G1_CROSS00DEG, &fd00[k*(G1_CROSS00DEG+1)],
                   G1_CROSS01DEG, &fd01[k*(G1_CROSS01DEG+1)],
                   G1_CROSS10DEG, &fd10[k*(G1_CROSS10DEG+1)],
                   G1_CROSS11DEG, &fd11[k*(G1_CROSS11DEG+1)],
                   &degu, &degv, y );

    for ( l = j = 0;  l < G1H_FINALDEG-3;  l++ )
      for ( m = 0;  m < G1H_FINALDEG-3;  m++, j++ )
        pkn_SubtractMatrixf ( 1, spdimen, 0, &y[((l+2)*(G1H_FINALDEG+1)+(m+2))*spdimen],
                              0, &x[(i*G1_DBDIM+j)*spdimen],
                              0, &y[((l+2)*(G1H_FINALDEG+1)+(m+2))*spdimen] );

    outpatch ( degu, degv, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_OutputExtPatchesf*/

boolean g1h_ExtFillHolef ( GHoleDomainf *domain,
                           int spdimen, CONST_ float *hole_cp,
                           float *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  float  *fc00, *x;
  int    hole_k, nfunc_a, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->EAmat )
    if ( !g1h_DecomposeExtMatrixf ( domain ) )
      goto failure;

  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ExtFillHolef*/

boolean g1h_ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii, *cmat, *rcmat;
  float  *fc00, *b, *x, *y;
  int    hole_k, nfunc_a, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->ECmat || !privateG1->ERCmat ||
       privateG1->extnconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  cmat = privateG1->ECmat;
  rcmat = privateG1->ERCmat;
  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
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
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixf ( nconstr, nbf, nbf, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, y );
  pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemf ( spdimen*nbf );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ExtFillHoleConstrf*/

boolean g1h_ExtFillHoleAltConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii,*acmat, *arcmat;
  float  *fc00, *b, *x, *y;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->AECmat || !privateG1->AERCmat ||
       privateG1->extnaconstr != naconstr || privateG1->extacdim != spdimen ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  acmat = privateG1->AECmat;
  arcmat = privateG1->AERCmat;
  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
  pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, x );
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
  pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, b );
  pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, Lii, spdimen, spdimen, b );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, b, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ExtFillHoleAltConstrf*/

