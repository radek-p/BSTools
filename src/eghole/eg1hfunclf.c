
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

static boolean g1h_ComputeBBMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  int      hole_k, nfunc_b;
  int      i, j, l, n;
  float    *trd, *jac, *BBmat;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float    *lap;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  if ( !privateG || !privateG1 )
    goto failure;
  if ( !privateG1->Amat && !privateG1->EAmat )
    goto failure;

  if ( !privateG1->BBmat ) {
    nfunc_b   = privateG1->nfunc_b;
    support_b = privateG->support_b;
    BBmat = privateG1->BBmat = malloc ( (nfunc_b*(nfunc_b+1)/2)*sizeof(float) );
    if ( !BBmat )
      goto failure;

    n = hole_k*G1_NQUADSQ;

    tkn = pkv_GetScratchMemf ( 19*G1_NQUAD );
    jac = pkv_GetScratchMemf ( n );
    trd = pkv_GetScratchMemf ( 5*n );
    lap = pkv_GetScratchMemf ( (nfunc_b+1)*n );
    if ( !tkn || !jac || !trd || !lap )
      goto failure;
    hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[6*G1_NQUAD];
    ddhfunc = &dhfunc[6*G1_NQUAD];
    _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
    mbs_TabCubicHFuncDer2f ( 0.0, 1.0, G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc );

    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetDiPatchCurvesf ( domain, i,
                               &c00, &c01, &c10, &c11,
                               &d00, &d01, &d10, &d11 );
      _g1h_TabDiPatchJac2f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*5] );
    }
    for ( j = 0; j < nfunc_b; j++ )
      for ( i = 0; i < hole_k; i++ )
        if ( support_b[j] & (0x0001 << i) ) {
          _g1h_GetBFBPatchCurvesf ( domain, j, i,
                                    &fc00, &fc01, &fc10, &fc11,
                                    &fd00, &fd01, &fd10, &fd11 );
          _g1h_TabLaplacianf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                               fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                               &trd[i*G1_NQUADSQ*5], &lap[(j*hole_k+i)*G1_NQUADSQ] );
        }
    for ( i = l = 0;  i < nfunc_b;  i++ ) {
      for ( j = 0;  j <= i;  j++, l++ )
        BBmat[l] = _g1h_Integralf ( hole_k, G1_NQUADSQ, jac,
                                    support_b[i], &lap[i*n],
                                    support_b[j], &lap[j*n] );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeBBMatrixf*/

float g1h_FunctionalValuef ( GHoleDomainf *domain, int spdimen,
                             const float *hole_cp, const float *acoeff )
{
  void   *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  unsigned char *bfcpn;
  int    nfunc_a, nfunc_b, j;
  double fval;
  float  *Amat, *Bmat, *BBmat, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g1h_ComputeBBMatrixf ( domain ) )
    goto failure;

  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  bfcpn   = privateG->bfcpn;
  Amat    = privateG1->Amat;
  Bmat    = privateG1->Bmat;
  BBmat   = privateG1->BBmat;
  if ( !Amat || !Bmat )
    goto failure;

  a = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  c = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(float) );
  pkn_SymMatrixMultf ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductf ( spdimen*nfunc_b, a, c );

  pkn_SymMatrixMultf ( nfunc_a, Amat, spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bmat,
                    spdimen, spdimen, a, spdimen, b );
  fval += pkn_ScalarProductf ( spdimen*nfunc_a, acoeff, c )
           -2.0*pkn_ScalarProductf ( spdimen*nfunc_a, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return (float)max ( 0.0, 0.25*fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g1h_FunctionalValuef*/

/*
float _g1h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                CONST_ float *hole_cp, CONST_ float *acoeff )
{
  void   *sp;
  G1HolePrivateRecf *privateG1;
  unsigned char *bfcpn;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf, j;
  double fval;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *BBmat, *Lii, *Lki, *Lkk, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g1h_ComputeBBMatrixf ( domain ) )
    goto failure;

  privateG1 = domain->privateG1;
  if ( !_g1h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk,
                                     &Bi, &Bk, &Lii, &Lki, &Lkk ) )
    goto failure;

  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  bfcpn   = privateG1->bfcpn;
  BBmat   = privateG1->BBmat;

  a = pkv_GetScratchMemf ( spdimen*nbf );
  b = pkv_GetScratchMemf ( spdimen*nbf );
  c = pkv_GetScratchMemf ( spdimen*nbf );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(float) );
  pkn_SymMatrixMultf ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductf ( spdimen*nfunc_b, a, c );

  pkn_BSymMatrixMultf ( hole_k, G1_DBDIM, nfunc_a, Aii, Aki, Akk,
                        spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixf ( nfunc_c, nfunc_b, nfunc_b, Bi,
                    spdimen, spdimen, a, spdimen, b );
  pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bk,
                    spdimen, spdimen, a, spdimen, &b[nfunc_c*spdimen] );
  fval += pkn_ScalarProductf ( spdimen*nbf, acoeff, c )
            -2.0*pkn_ScalarProductf ( spdimen*nbf, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return max ( 0.0, fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
}*/ /*g1h_ExtFunctionalValuef*/

/* The function above returns incorrect results. Below is a temporary */
/* replacement, to be used before the function above is fixed.  */

static int _np, _spdimen;
static float *patches;
static void myoutpatch ( int n, int m, const float *cp, void *usrptr )
{
  memcpy ( &patches[_np*_spdimen*(n+1)*(m+1)], cp,
           (n+1)*(m+1)*_spdimen*sizeof(float) );
  _np ++;
} /*myoutpatch*/

float g1h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                CONST_ float *hole_cp, CONST_ float *acoeff )
{
  void     *sp;
  G1HolePrivateRecf *privateG1;
  float   *fc00, *x, *Bi, *Bk, fct;
  float   *tkn, jac, lap, lp;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  point2f  *bdi;
  vector2f di, diu, div, diuu, diuv, divv;
  float   *pi, *piu, *piv, *piuu, *piuv, *pivv;
  float   *fi, *fiu, *fiv, *fiuu, *fiuv, *fivv;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  int      i, j, k, l, dn, dm;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  if ( !privateG1->EAmat )
    if ( !g1h_DecomposeExtMatrixf ( domain ) )
      goto failure;

  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+6)*2*hole_k*spdimen );
  patches = pkv_GetScratchMemf ( hole_k*spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  tkn = pkv_GetScratchMemf ( G1_NQUAD );
  bdi = pkv_GetScratchMem ( (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(point2f) );
  pi  = pkv_GetScratchMemf ( 12*spdimen );
  if ( !x || !fc00 || !patches|| !tkn || !bdi || !pi )
    goto failure;
  piu = &pi[spdimen];    piv = &piu[spdimen];
  piuu = &piv[spdimen];  piuv = &piuu[spdimen];  pivv = &piuv[spdimen];
  fi = &pivv[spdimen];   fiu = &fi[spdimen];     fiv = &fiu[spdimen];
  fiuu = &fiv[spdimen];  fiuv = &fiuu[spdimen];  fivv = &fiuv[spdimen];

  Bi = privateG1->EBmat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  _spdimen = spdimen;
  _np = 0;
  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, acoeff, fc00, NULL,
                                 myoutpatch ) )
    goto failure;

  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
  fct = 0.0;
  for ( k = 0; k < hole_k; k++ ) {
    _g1h_GetDiPatchCurvesf ( domain, k, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    if ( !mbs_BezC1CoonsToBezf ( 2,
               G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
               G1_CROSS10DEG, (float*)c10, G1_CROSS11DEG, (float*)c11,
               G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
               G1_CROSS10DEG, (float*)d10, G1_CROSS11DEG, (float*)d11,
               &dn, &dm, (float*)bdi ) )
      goto failure;
    for ( i = 0; i < G1_NQUAD; i++ )
      for ( j = 0; j < G1_NQUAD; j++ ) {
        mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)bdi,
                             tkn[i], tkn[j],
                             &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
        mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, spdimen,
                             &patches[k*(G1H_FINALDEG+1)*(G1H_FINALDEG+1)*spdimen],
                             tkn[i], tkn[j],
                             pi, piu, piv, piuu, piuv, pivv );
        if ( !pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                        diuv.x, diuv.y, divv.x, divv.y,
                                        spdimen, piu, piv, piuu, piuv, pivv,
                                        fiu, fiv, fiuu, fiuv, fivv ) )
          goto failure;
        jac = diu.x*div.y - diu.y*div.x;
        for ( l = 0, lap = 0.0;  l < spdimen;  l++ ) {
          lp = fiuu[l]+fivv[l];
          lap += lp*lp;
        }
        fct += jac*lap;
      }
  }

  pkv_SetScratchMemTop ( sp );
  return (float)(fct/(4.0*(float)G1_NQUADSQ));

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g1h_ExtFunctionalValuef*/

