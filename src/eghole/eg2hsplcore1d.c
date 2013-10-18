
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_DestroySPrivateDatad ( GHoleDomaind *domain )
{
  G2HoleSPrivateRecd *sprivate;

  if ( (sprivate = domain->SprivateG2) ) {
    if ( sprivate->omcknots ) free ( sprivate->omcknots );
    if ( sprivate->pvknots )  free ( sprivate->pvknots );
    if ( sprivate->pvvknots ) free ( sprivate->pvvknots );
    if ( sprivate->cknots )   free ( sprivate->cknots );
    if ( sprivate->basis_d )  free ( sprivate->basis_d );
    if ( sprivate->fpknots )  free ( sprivate->fpknots );
    if ( sprivate->SAMat )    free ( sprivate->SAMat );
    if ( sprivate->SBMat )    free ( sprivate->SBMat );
    if ( sprivate->SLMat )    free ( sprivate->SLMat );
    if ( sprivate->SCmat )    free ( sprivate->SCmat );   
    if ( sprivate->SRCmat )   free ( sprivate->SRCmat );  
    if ( sprivate->ASCmat )   free ( sprivate->ASCmat );  
    if ( sprivate->ASRCmat )  free ( sprivate->ASRCmat ); 
    free ( sprivate );
  }
  domain->SprivateG2 = NULL;
} /*g2h_DestroySPrivateDatad*/

static boolean g2h_AllocSPrivateDatad ( GHoleDomaind *domain, int nk, int m1, int m2 )
{
  G2HoleSPrivateRecd *sprivate;

  g2h_DestroySPrivateDatad ( domain );
  sprivate = domain->SprivateG2 = malloc ( sizeof(G2HoleSPrivateRecd) );
  if ( !sprivate )
    return false;
  sprivate->nk = nk;
  sprivate->m1 = m1;
  sprivate->m2 = m2;
  sprivate->omcknots = sprivate->pvknots = sprivate->pvvknots =
  sprivate->cknots = sprivate->basis_d = sprivate->fpknots = NULL;
  sprivate->SAMat = sprivate->SBMat = sprivate->SLMat = NULL;
  sprivate->SCmat = sprivate->SRCmat = sprivate->ASCmat = sprivate->ASRCmat = NULL;
  return true;
} /*g2h_AllocSPrivateDatad*/

static boolean _g2h_ConstructSplKnotSequencesd ( GHoleDomaind *domain )
{
  G2HoleSPrivateRecd *sprivate;
  int    nk, m1, m2, m3, m4, m5;
  int    i, j, k0, k1, k2, k3, k4;
  double *knr, *knpv, *knpvv, *knb, *cknots, *fpknots;
  double kn;

  sprivate = domain->SprivateG2;

  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  m3 = m1+G2_AUXDEG0;
  m4 = m1+G2_AUXDEG1;
  m5 = max ( m2, m3+(G2H_FINALDEG-G2_CROSS01DEG) );
  m5 = max ( m5, m4+(G2H_FINALDEG-G2_CROSS02DEG) );
  sprivate->lastomcknot = 2*G2_CROSS00DEG+1+nk*m1;
  sprivate->lastpvknot  = 2*G2_CROSS01DEG+1+nk*m3;
  sprivate->lastpvvknot = 2*G2_CROSS02DEG+1+nk*m4;
  sprivate->lastcknot   = 2*G2H_FINALDEG+1+nk*m2;
  sprivate->lastfpknot  = 2*G2H_FINALDEG+1+nk*m5;
  knr = sprivate->omcknots = malloc ( (sprivate->lastomcknot+1)*sizeof(double) );
  knpv = sprivate->pvknots = malloc ( (sprivate->lastpvknot+1)*sizeof(double) );
  knpvv = sprivate->pvvknots = malloc ( (sprivate->lastpvvknot+1)*sizeof(double) );
  cknots = sprivate->cknots = malloc ( (sprivate->lastcknot+1)*sizeof(double) );
  fpknots = sprivate->fpknots = malloc ( (sprivate->lastfpknot+1)*sizeof(double) );
  if ( !knr || !knpv || !knpvv || !cknots || !fpknots )
    return false;

        /* set the knot pattern for the Bezier polynomials */
  knb = &sprivate->bezknots[0];
  kn = 0.0;
  mbs_SetKnotPatternd ( 0, &kn, G2H_FINALDEG+1, &j, knb );
  kn = 1.0;
  mbs_SetKnotPatternd ( 0, &kn, G2H_FINALDEG+1, &j, &knb[G2H_FINALDEG+1] );
  sprivate->lastbezknot = 2*G2H_FINALDEG+1;
        /* set the knot sequences for the auxiliary basis function patches, */
        /* the basis functions from the blocks C, D, and the final patches. */
  kn = 0.0;
  mbs_SetKnotPatternd ( 0, &kn, (k0 = G2_CROSS00DEG+1), &j, knr );
  mbs_SetKnotPatternd ( 0, &kn, (k1 = G2_CROSS01DEG+1), &j, knpv );
  mbs_SetKnotPatternd ( 0, &kn, (k2 = G2_CROSS02DEG+1), &j, knpvv );
  mbs_SetKnotPatternd ( 0, &kn, (k3 = G2H_FINALDEG+1), &j, cknots );
  mbs_SetKnotPatternd ( 0, &kn, (k4 = G2H_FINALDEG+1), &j, fpknots );
  for ( i = 1; i <= nk; i++ ) {
    kn = (double)i/(double)(nk+1);
    mbs_SetKnotPatternd ( 0, &kn, m1, &j, &knr[k0] );      k0 += m1;
    mbs_SetKnotPatternd ( 0, &kn, m3, &j, &knpv[k1] );     k1 += m3;
    mbs_SetKnotPatternd ( 0, &kn, m4, &j, &knpvv[k2] );    k2 += m4;
    mbs_SetKnotPatternd ( 0, &kn, m2, &j, &cknots[k3] );   k3 += m2;
    mbs_SetKnotPatternd ( 0, &kn, m5, &j, &fpknots[k4] );  k4 += m5;
  }
  kn = 1.0;
  mbs_SetKnotPatternd ( 0, &kn, G2_CROSS00DEG+1, &j, &knr[k0] );
  mbs_SetKnotPatternd ( 0, &kn, G2_CROSS01DEG+1, &j, &knpv[k1] );
  mbs_SetKnotPatternd ( 0, &kn, G2_CROSS02DEG+1, &j, &knpvv[k2] );
  mbs_SetKnotPatternd ( 0, &kn, G2H_FINALDEG+1, &j, &cknots[k3] );
  mbs_SetKnotPatternd ( 0, &kn, G2H_FINALDEG+1, &j, &fpknots[k4] );

  return true;
} /*_g2h_ConstructSplKnotSequencesd*/

static boolean _g2h_FindSplCrossDerDd ( const double *knb,
              const double *b1, const double *c1,
              const double *b2, const double *c2,
              const double *b1b1, const double *twob1c1, const double *c1c1,
              int nzc, int ku, int lknr, const double *knr,
              const double *r0, const double *r0s, const double *r0ss,
              int degpv, int lknpv, const double *knpv, double *pv,
              int degpvv, int lknpvv, const double *knpvv, double *pvv )
{
/* knb - array with 10 zeros followed by 10 ones */
/* b1,c1,b2,c2,b1b1,twob1c1,c1c1 - junction functions, polynomials, */
/* nzc - number of the nonzero function, 0,1 or 2 */
/* ku - number of intervals between knots in the knr sequence */
/* lknr, knr - knot sequence for the auxiliary basis function patch and */
/* its derivatives; r0 is of degree G2H_OMCDEG, r0s and r0u are of */
/* degree G2H_OMCDEG-1, r0ss, r0us, r0uu are of degree G2H_OMCDEG-2, */

  void   *sp;
  double *r0u, *r0uu, *r0us;
  double *aux1, *aux2, *aux3, *kn1, *kn2, *kn3;
  int    size, deg1, lkn1, deg2, lkn2, deg3, lkn3;

  sp = pkv_GetScratchMemTop ();
  r0u = pkv_GetScratchMemd ( 3*(lknr-G2H_OMCDEG)-5 );
  size = lknr+1+ku*(G2H_FINALDEG-G2H_OMCDEG+2);
  aux1 = pkv_GetScratchMemd ( 6*size );
  if ( !r0u || !aux1 )
    goto failure;
  r0uu = &r0u[lknr-G2H_OMCDEG-1];
  r0us = &r0uu[lknr-G2H_OMCDEG-2];
  aux2 = &aux1[size];
  aux3 = &aux2[size];
  kn1 = &aux3[size];
  kn2 = &kn1[size];
  kn3 = &kn2[size];

        /* use the information, which curve is nonzero to avoid */
        /* wasting time on zero milling - more than 3 times speedup */
  switch ( nzc ) {
case 0:   /* nonzero is r0 and its derivatives, r0u, r0uu */
    mbs_FindBSDerivativeC1d ( G2H_OMCDEG, lknr, knr, r0, NULL, NULL, r0u );
    mbs_FindBSDerivativeC1d ( G2H_OMCDEG-1, lknr-2, &knr[1], r0u, NULL, NULL, r0uu );
    mbs_FindBSDerivativeC1d ( G2H_OMCDEG-1, lknr-2, &knr[1], r0s, NULL, NULL, r0us );

      /* compute pv = b1*r0u + c1*r0s = b1*r0u */
    lkn1 = size;
    if ( !mbs_multiMultBSCd ( 1, G2_BF01DEG, 2*G2_BF01DEG+1,
                        &knb[G2H_FINALDEG-G2_BF01DEG],
                        0, b1, 1, 1, G2H_OMCDEG-1, lknr-2, &knr[1], 0, r0u,
                        &deg1, &lkn1, kn1, 0, aux1 ) )
      goto failure;
    mbs_multiAdjustBSCRepd ( 1, 1, deg1, lkn1, kn1, 0, aux1,
                             degpv, lknpv, knpv, 0, pv );

      /* compute pvv = b2*r0u + c2*r0s + b1b1*r0uu + twob1c1*r0us + c1c1*r0ss */
      /* = b2*r0u + b1b1*r0uu */
    lkn1 = size;
    if ( !mbs_multiMultBSCd ( 1, G2_BF02DEG, 2*G2_BF02DEG+1,
                        &knb[G2H_FINALDEG-G2_BF02DEG],
                        0, b2, 1, 1, G2H_OMCDEG-1, lknr-2, &knr[1], 0, r0u,
                        &deg1, &lkn1, kn1, 0, aux1 ) )
      goto failure;
    lkn2 = size;
    if ( !mbs_multiMultBSCd ( 1, 2*G2_BF01DEG, 4*G2_BF01DEG+1,
                        &knb[G2H_FINALDEG-2*G2_BF01DEG],
                        0, b1b1, 1, 1, G2H_OMCDEG-2, lknr-4, &knr[2], 0, r0uu,
                        &deg2, &lkn2, kn2, 0, aux2 ) )
      goto failure;
    if ( !mbs_multiAddBSCurvesd ( 1, 1, deg1, lkn1, kn1, 0, aux1,
                                  deg2, lkn2, kn2, 0, aux2,
                                  &deg3, &lkn3, kn3, 0, aux3 ) )
      goto failure;
    mbs_multiAdjustBSCRepd ( 1, 1, deg3, lkn3, kn3, 0, aux3,
                             degpvv, lknpvv, knpvv, 0, pvv );
    break;

case 1:   /* nonzero is r0s and its derivative, r0us */
    mbs_FindBSDerivativeC1d ( G2H_OMCDEG-1, lknr-2, &knr[1], r0s, NULL, NULL, r0us );

      /* compute pv = b1*r0u + c1*r0s = c1*r0s */
    lkn1 = size;
    if ( !mbs_multiMultBSCd ( 1, G2_CG01DEG, 2*G2_CG01DEG+1,
                        &knb[G2H_FINALDEG-G2_CG01DEG],
                        0, c1, 1, 1, G2H_OMCDEG-1, lknr-2, &knr[1], 0, r0s,
                        &deg1, &lkn1, kn1, 0, aux1 ) )
      goto failure;
    mbs_multiAdjustBSCRepd ( 1, 1, deg1, lkn1, kn1, 0, aux1,
                             degpv, lknpv, knpv, 0, pv );

      /* compute pvv = b2*r0u + c2*r0s + b1b1*r0uu + twob1c1*r0us + c1c1*r0ss */
      /* = c2*r0s + twob1c1*r0us */
    lkn1 = size;
    if ( !mbs_multiMultBSCd ( 1, G2_CG02DEG, 2*G2_CG02DEG+1,
                        &knb[G2H_FINALDEG-G2_CG02DEG],
                        0, c2, 1, 1, G2H_OMCDEG-1, lknr-2, &knr[1], 0, r0s,
                        &deg1, &lkn1, kn1, 0, aux1 ) )
      goto failure;
    lkn2 = size;
    if ( !mbs_multiMultBSCd ( 1, G2_BF01DEG+G2_CG01DEG, 2*(G2_BF01DEG+G2_CG01DEG)+1,
                        &knb[G2H_FINALDEG-G2_BF01DEG-G2_CG01DEG], 0, twob1c1,
                        1, 1, G2H_OMCDEG-2, lknr-4, &knr[2], 0, r0us,
                        &deg2, &lkn2, kn2, 0, aux2 ) )
      goto failure;
    if ( !mbs_multiAddBSCurvesd ( 1, 1, deg1, lkn1, kn1, 0, aux1,
                                  deg2, lkn2, kn2, 0, aux2,
                                  &deg3, &lkn3, kn3, 0, aux3 ) )
      goto failure;
    mbs_multiAdjustBSCRepd ( 1, 1, deg3, lkn3, kn3, 0, aux3,
                             degpvv, lknpvv, knpvv, 0, pvv );
    break;

case 2:   /* nonzero is r0ss */
      /* compute pv = b1*r0u + c1*r0s = 0 */
    memset ( pv, 0, (lknpv-degpv)*sizeof(double) );

      /* compute pvv = b2*r0u + c2*r0s + b1b1*r0uu + twob1c1*r0us + c1c1*r0ss */
      /* = c1c1*r0ss */
    lkn1 = size;
    if ( !mbs_multiMultBSCd ( 1, 2*G2_CG01DEG, 4*G2_CG01DEG+1,
                        &knb[G2H_FINALDEG-2*G2_CG01DEG],
                        0, c1c1, 1, 1, G2H_OMCDEG-2, lknr-4, &knr[2], 0, r0ss,
                        &deg1, &lkn1, kn1, 0, aux1 ) )
      goto failure;
    mbs_multiAdjustBSCRepd ( 1, 1, deg1, lkn1, kn1, 0, aux1,
                             degpvv, lknpvv, knpvv, 0, pvv );
    break;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_FindSplCrossDerDd*/

boolean _g2h_GetSplDBasisAuxpd ( GHoleDomaind *domain, int fn, int cn,
                   int *nzc, double *fcomc, double *fcomcd, double *fcomcdd )
{
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int                nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int                nk, m1, i, j;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  fn -= nfunc_a+nfunc_b+nfunc_c;
  if ( fn < 0 || fn >= nfunc_d )
    return false;

  memset ( fcomc, 0, (G2H_OMCDEG+1+nk*m1)*sizeof(double) );
  memset ( fcomcd, 0, (G2H_OMCDEG+nk*m1)*sizeof(double) );
  memset ( fcomcdd, 0, (G2H_OMCDEG-1+nk*m1)*sizeof(double) );
  if ( cn == fn/(3*nk*m1 ) ) {  /* the auxiliary basis function patches */
                                /* are nonzero for only one curve */
    i = fn % (3*nk*m1);
    j = i / 3;
    switch ( (*nzc = i % 3) ) {
  case 0: fcomc[j+5] = 1.0;    break;
  case 1: fcomcd[j+4] = 1.0;   break;
  case 2: fcomcdd[j+3] = 1.0;  break;
    }
  }

  return true;
} /*_g2h_GetSplDBasisAuxpd*/

boolean _g2h_GetSplDBasisCrossDerd ( GHoleDomaind *domain, int fn, int cn,
                 double *fcomc, double *pv, double *pvv, double *pu, double *puu )
{
/* for the i-th domain division curve and the basis fn-th function from the */
/* block D, get the coefficients representing the value at the curve and the */
/* cross derivatives of two basis function patches adjacent to the curve; */
/* the i-th in the arrays pv, pvv and (i-1) mod hole_k in the arrays pu, puu. */
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int                nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int                nk, m1, i, j, lpkn, lpvkn, lpvvkn;
  double             *basis_d;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  lpkn = sprivate->lastomcknot;
  lpvkn = sprivate->lastpvknot;
  lpvvkn = sprivate->lastpvvknot;

  fn -= nfunc_a+nfunc_b+nfunc_c;
  if ( fn < 0 || fn >= nfunc_d )
    return false;

  memset ( fcomc, 0, (lpkn-G2H_OMCDEG)*sizeof(double) );
  if ( cn == fn/(3*nk*m1) )  {
    i = fn % (3*nk*m1);
    if ( (i % 3) == 0 ) fcomc[i/3+5] = 1.0;
    basis_d = sprivate->basis_d;
    j = fn*2*(lpvkn+lpvvkn-G2_CROSS01DEG-G2_CROSS02DEG);
    memcpy ( pv, &basis_d[j], (lpvkn-G2_CROSS01DEG)*sizeof(double) );
    j += lpvkn-G2_CROSS01DEG;
    memcpy ( pvv, &basis_d[j], (lpvvkn-G2_CROSS02DEG)*sizeof(double) );
    j += lpvvkn-G2_CROSS02DEG;
    memcpy ( pu, &basis_d[j], (lpvkn-G2_CROSS01DEG)*sizeof(double) );
    j += lpvkn-G2_CROSS01DEG;
    memcpy ( puu, &basis_d[j], (lpvvkn-G2_CROSS02DEG)*sizeof(double) );
  }
  else {
    memset ( pv,  0, (lpvkn-G2_CROSS01DEG)*sizeof(double) );
    memset ( pvv, 0, (lpvvkn-G2_CROSS02DEG)*sizeof(double) );
    memset ( pu,  0, (lpvkn-G2_CROSS01DEG)*sizeof(double) );
    memset ( puu, 0, (lpvvkn-G2_CROSS02DEG)*sizeof(double) );
  }

  return true;
} /*_g2h_GetSplDBasisCrossDerd*/

static boolean _g2h_FindSplBasisFunctionsDd ( GHoleDomaind *domain )
{
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int    hole_k, nk, m1, i, j, k, l, fn, fnd, nzc;
  int    nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int    lomckn, lpvkn, lpvvkn;
  double *knpv, *knpvv, *knb, *knr, *fcomc, *fcomcd, *fcomcdd, *basis_d;
  double *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
         *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
         *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01;
  double *pv, *pvv;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  lomckn = sprivate->lastomcknot;
  lpvkn = sprivate->lastpvknot;
  lpvvkn = sprivate->lastpvvknot;
  knb = sprivate->bezknots;
  knr = sprivate->omcknots;
  knpv = sprivate->pvknots;
  knpvv = sprivate->pvvknots;

        /* allocate memory for the basis function patches cross derivatives */
  basis_d = sprivate->basis_d = malloc ( nfunc_d*2*
      (lpvkn+lpvvkn-G2_CROSS01DEG-G2_CROSS02DEG)*sizeof(double) );
  if ( !basis_d )
    goto failure;
        /* allocate memory for the spline auxiliary basis function patches */
  fcomc = pkv_GetScratchMemd ( 3*(G2H_OMCDEG+nk*m1) );
  if ( !fcomc )
    goto failure;
  fcomcd = &fcomc[G2H_OMCDEG+1+nk*m1];
  fcomcdd = &fcomcd[G2H_OMCDEG+nk*m1];

        /* now setup the basis function patches cross derivatives */
  G2GetPolynomialAddresses0 ( privateG2->jfunc,
        b01, c01, f01, g01, b11, c11, f11, g11,
        b02, c02, f02, g02, b12, c12, f12, g12,
        b01b01, twob01c01, c01c01, f01f01, twof01g01, g01g01 );
  for ( i = fnd = 0, fn = nfunc_a+nfunc_b+nfunc_c;  i < hole_k;  i++ ) {
    l = (i+hole_k-1) % hole_k;
    for ( j = 0; j < nk*m1; j++ )
      for ( k = 0;  k < 3;  k++, fnd++, fn++ ) {
        _g2h_GetSplDBasisAuxpd ( domain, fn, i, &nzc, fcomc, fcomcd, fcomcdd );
        pv = &basis_d[fnd*2*(lpvkn+lpvvkn-G2_CROSS01DEG-G2_CROSS02DEG)];
        pvv = &pv[lpvkn-G2_CROSS01DEG];
        if ( !_g2h_FindSplCrossDerDd ( knb,
                   &b01[i*(G2_BF01DEG+1)], &c01[i*(G2_CG01DEG+1)],
                   &b02[i*(G2_BF02DEG+1)], &c02[i*(G2_CG02DEG+1)],
                   &b01b01[i*(2*G2_BF01DEG+1)], &twob01c01[i*(G2_BF01DEG+G2_CG01DEG+1)],
                   &c01c01[i*(2*G2_CG01DEG+1)], nzc, nk+1, lomckn,
                   knr, fcomc, fcomcd, fcomcdd,
                   G2_CROSS01DEG, lpvkn, knpv, pv,
                   G2_CROSS02DEG, lpvvkn, knpvv, pvv ) )
          goto failure;
        pv = &pvv[lpvvkn-G2_CROSS02DEG];
        pvv = &pv[lpvkn-G2_CROSS01DEG];
        if ( !_g2h_FindSplCrossDerDd ( knb,
                   &f01[l*(G2_BF01DEG+1)], &g01[l*(G2_CG01DEG+1)],
                   &f02[l*(G2_BF02DEG+1)], &g02[l*(G2_CG02DEG+1)],
                   &f01f01[l*(2*G2_BF01DEG+1)], &twof01g01[l*(G2_BF01DEG+G2_CG01DEG+1)],
                   &g01g01[l*(2*G2_CG01DEG+1)], nzc, nk+1, lomckn,
                   knr, fcomc, fcomcd, fcomcdd,
                   G2_CROSS01DEG, lpvkn, knpv, pv,
                   G2_CROSS02DEG, lpvvkn, knpvv, pvv ) )
          goto failure;
      }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_FindSplBasisFunctionsDd*/

static boolean _g2h_FindSplBasisFunctionsd ( GHoleDomaind *domain )
{
  void               *sp;
  G2HoleSPrivateRecd *sprivate;
  int                nk, m1, m2, hole_k;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  sprivate = domain->SprivateG2;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
        /* compute the numbers of spline basis functions */
  sprivate->csize = (G2H_FINALDEG-5+nk*m2)*(G2H_FINALDEG-5+nk*m2);
  sprivate->nsfunc_c = hole_k*sprivate->csize;
  sprivate->dsize = 3*nk*m1;
  sprivate->nsfunc_d = hole_k*sprivate->dsize;
        /* construct the knot sequence for the representation of the */
        /* final patches */

  if ( !_g2h_ConstructSplKnotSequencesd ( domain ) )
    goto failure;
        /* find the basis function patches for the block D */
  if ( !_g2h_FindSplBasisFunctionsDd ( domain ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_FindSplBasisFunctionsd*/

boolean g2h_ComputeSplBasisd ( GHoleDomaind *domain, int nk, int m1, int m2 )
{
  if ( nk < 1 || nk > G2H_S_MAX_NK ||
       m1 < 1 || m1 > G2H_S_MAX_M1 ||
       m2 < 1 || m2 > G2H_S_MAX_M2 )
    return false;

  if ( !g2h_AllocSPrivateDatad ( domain, nk, m1, m2 ) )
    return false;

        /* the first step is the construction of the polynomial basis, */
        /* as the junction functions and domain patches will be the same */
  if ( !g2h_ComputeBasisd ( domain ) )
    return false;
  if ( !_g2h_FindSplBasisFunctionsd ( domain ) )
    return false;
  return true;
} /*g2h_ComputeSplBasisd*/

