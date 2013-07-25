
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
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

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DestroySPrivateDatad ( GHoleDomaind *domain )
{
  G1HoleSPrivateRecd *sprivate;

  if ( (sprivate = domain->SprivateG1) ) {
    if ( sprivate->omcknots )  free ( sprivate->omcknots );
    if ( sprivate->pvknots )   free ( sprivate->pvknots );
    if ( sprivate->cknots )    free ( sprivate->cknots );
    if ( sprivate->basis_d )   free ( sprivate->basis_d );
    if ( sprivate->fpknots )   free ( sprivate->fpknots );
    if ( sprivate->SAMat )     free ( sprivate->SAMat );
    if ( sprivate->SBMat )     free ( sprivate->SBMat );
    if ( sprivate->SLMat )     free ( sprivate->SLMat );
    if ( sprivate->Q2SAMat )   free ( sprivate->Q2SAMat );
    if ( sprivate->Q2SBMat )   free ( sprivate->Q2SBMat );
    if ( sprivate->Q2SLMat )   free ( sprivate->Q2SLMat );
    if ( sprivate->SCmat )     free ( sprivate->SCmat );
    if ( sprivate->SRCmat )    free ( sprivate->SRCmat );
    if ( sprivate->ASCmat )    free ( sprivate->ASCmat );
    if ( sprivate->ASRCmat )   free ( sprivate->ASRCmat );
    if ( sprivate->Q2SRCmat )  free ( sprivate->Q2SRCmat );
    if ( sprivate->Q2SARCmat ) free ( sprivate->Q2SARCmat );
    free ( sprivate );
  }
  domain->SprivateG1 = NULL;
} /*g1h_DestroySPrivateDatad*/

static boolean g1h_AllocSPrivateDatad ( GHoleDomaind *domain, int nk, int m1, int m2 )
{
  G1HoleSPrivateRecd *sprivate;

  g1h_DestroySPrivateDatad ( domain );
  sprivate = domain->SprivateG1 = malloc ( sizeof(G1HoleSPrivateRecd) );
  if ( !sprivate )
    return false;
  sprivate->nk = nk;
  sprivate->m1 = m1;
  sprivate->m2 = m2;
  sprivate->omcknots = sprivate->pvknots =
  sprivate->cknots = sprivate->basis_d = sprivate->fpknots = NULL;
  sprivate->SAMat = sprivate->SBMat = sprivate->SLMat = NULL;
  sprivate->Q2SAMat = sprivate->Q2SBMat = sprivate->Q2SLMat = NULL;
  sprivate->SCmat = sprivate->SRCmat = sprivate->ASCmat = sprivate->ASRCmat = NULL;
  sprivate->Q2SRCmat = sprivate->Q2SARCmat = NULL;
  return true;
} /*g1h_AllocSPrivateDatad*/

static boolean _g1h_ConstructSplKnotSequencesd ( GHoleDomaind *domain )
{
  G1HoleSPrivateRecd *sprivate;
  int    nk, m1, m2, m3, m5;
  int    i, j, k0, k1, k3, k4;
  double *knr, *knpv, *knb, *cknots, *fpknots;
  double  kn;

  sprivate = domain->SprivateG1;

  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  m3 = m1+G1_AUXDEG0;
  m5 = max ( m2, m3+(G1H_FINALDEG-G1_CROSS01DEG) );
  sprivate->lastomcknot = 2*G1_CROSS00DEG+1+nk*m1;
  sprivate->lastpvknot  = 2*G1_CROSS01DEG+1+nk*m3;
  sprivate->lastcknot   = 2*G1H_FINALDEG+1+nk*m2;
  sprivate->lastfpknot  = 2*G1H_FINALDEG+1+nk*m5;
  knr = sprivate->omcknots = malloc ( (sprivate->lastomcknot+1)*sizeof(double) );
  knpv = sprivate->pvknots = malloc ( (sprivate->lastpvknot+1)*sizeof(double) );
  cknots = sprivate->cknots = malloc ( (sprivate->lastcknot+1)*sizeof(double) );
  fpknots = sprivate->fpknots = malloc ( (sprivate->lastfpknot+1)*sizeof(double) );
  if ( !knr || !knpv || !cknots || !fpknots )
    return false;

        /* set the knot pattern for the Bezier polynomials */
  knb = &sprivate->bezknots[0];
  kn = 0.0;
  mbs_SetKnotPatternd ( 0, &kn, G1H_FINALDEG+1, &j, knb );
  kn = 1.0;
  mbs_SetKnotPatternd ( 0, &kn, G1H_FINALDEG+1, &j, &knb[G1H_FINALDEG+1] );
  sprivate->lastbezknot = 2*G1H_FINALDEG+1;
        /* set the knot sequences for the auxiliary basis function patches, */
        /* the basis functions from the blocks C, D, and the final patches. */
  kn = 0.0;
  mbs_SetKnotPatternd ( 0, &kn, (k0 = G1_CROSS00DEG+1), &j, knr );
  mbs_SetKnotPatternd ( 0, &kn, (k1 = G1_CROSS01DEG+1), &j, knpv );
  mbs_SetKnotPatternd ( 0, &kn, (k3 = G1H_FINALDEG+1), &j, cknots );
  mbs_SetKnotPatternd ( 0, &kn, (k4 = G1H_FINALDEG+1), &j, fpknots );
  for ( i = 1; i <= nk; i++ ) {
    kn = (double)i/(double)(nk+1);
    mbs_SetKnotPatternd ( 0, &kn, m1, &j, &knr[k0] );      k0 += m1;
    mbs_SetKnotPatternd ( 0, &kn, m3, &j, &knpv[k1] );     k1 += m3;
    mbs_SetKnotPatternd ( 0, &kn, m2, &j, &cknots[k3] );   k3 += m2;
    mbs_SetKnotPatternd ( 0, &kn, m5, &j, &fpknots[k4] );  k4 += m5;
  }
  kn = 1.0;
  mbs_SetKnotPatternd ( 0, &kn, G1_CROSS00DEG+1, &j, &knr[k0] );
  mbs_SetKnotPatternd ( 0, &kn, G1_CROSS01DEG+1, &j, &knpv[k1] );
  mbs_SetKnotPatternd ( 0, &kn, G1H_FINALDEG+1, &j, &cknots[k3] );
  mbs_SetKnotPatternd ( 0, &kn, G1H_FINALDEG+1, &j, &fpknots[k4] );

  return true;
} /*_g1h_ConstructSplKnotSequencesd*/

static boolean _g1h_FindSplCrossDerDd ( const double *knb,
              const double *b1, const double *c1,
              int nzc, int ku, int lknr, const double *knr,
              const double *r0, const double *r0s,
              int degpv, int lknpv, const double *knpv, double *pv )
{
/* knb - array with 5 zeros followed by 5 ones, */
/* b1,c1 - junction functions, polynomials, */
/* nzc - number of the nonzero function, 0 or 1 */
/* ku - number of intervals between knots in the knr sequence */
/* lknr, knr - knot sequence for the auxiliary basis function patch and */
/* its derivatives; r0 is of degree G1H_OMCDEG, and r0s and r0u are of */
/* degree G1H_OMCDEG-1. */

  void   *sp;
  double *r0u;
  double *aux, *kn;
  int    size, deg, lkn;

  sp = pkv_GetScratchMemTop ();
  r0u = pkv_GetScratchMemd ( 2*(lknr-G1H_OMCDEG)-3 );
  size = lknr+1+ku*(G1H_FINALDEG-G1H_OMCDEG+2);
  aux = pkv_GetScratchMemd ( 2*size );
  if ( !r0u || !aux )
    goto failure;
  kn = &aux[size];

        /* use the information, which curve is nonzero to avoid */
        /* wasting time on zero milling */
  switch ( nzc ) {
case 0:   /* nonzero is r0 and its derivative, r0u */
    mbs_FindBSDerivativeC1d ( G1H_OMCDEG, lknr, knr, r0, NULL, NULL, r0u );

      /* compute pv = b1*r0u + c1*r0s = b1*r0u */
    lkn = size;
    mbs_multiMultBSCd ( 1, G1_BF01DEG, 2*G1_BF01DEG+1, &knb[G1H_FINALDEG-G1_BF01DEG],
                        0, b1, 1, 1, G1H_OMCDEG-1, lknr-2, &knr[1], 0, r0u,
                        &deg, &lkn, kn, 0, aux );
    mbs_multiAdjustBSCRepd ( 1, 1, deg, lkn, kn, 0, aux,
                             degpv, lknpv, knpv, 0, pv );
    break;

case 1:   /* nonzero is r0s */
      /* compute pv = b1*r0u + c1*r0s = c1*r0s */
    lkn = size;
    mbs_multiMultBSCd ( 1, G1_CG01DEG, 2*G1_CG01DEG+1, &knb[G1H_FINALDEG-G1_CG01DEG],
                        0, c1, 1, 1, G1H_OMCDEG-1, lknr-2, &knr[1], 0, r0s,
                        &deg, &lkn, kn, 0, aux );
    mbs_multiAdjustBSCRepd ( 1, 1, deg, lkn, kn, 0, aux,
                             degpv, lknpv, knpv, 0, pv );
    break;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_FindSplCrossDerDd*/

boolean _g1h_GetSplDBasisAuxpd ( GHoleDomaind *domain, int fn, int cn,
                                 int *nzc, double *fcomc, double *fcomcd )
{
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int                nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int                nk, m1, i, j;

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  fn -= nfunc_a+nfunc_b+nfunc_c;
  if ( fn < 0 || fn >= nfunc_d )
    return false;

  memset ( fcomc, 0, (G1H_OMCDEG+1+nk*m1)*sizeof(double) );
  memset ( fcomcd, 0, (G1H_OMCDEG+nk*m1)*sizeof(double) );
  if ( cn == fn/(2*nk*m1 ) ) {  /* the auxiliary basis function patches */
                                /* are nonzero for only one curve */
    i = fn % (2*nk*m1);
    j = i / 2;
    switch ( (*nzc = i % 2) ) {
  case 0: fcomc[j+3] = 1.0;    break;
  case 1: fcomcd[j+2] = 1.0;   break;
    }
  }

  return true;
} /*_g1h_GetSplDBasisAuxpd*/

boolean _g1h_GetSplDBasisCrossDerd ( GHoleDomaind *domain, int fn, int cn,
                                     double *fcomc, double *pv, double *pu )
{
/* for the i-th domain division curve and the basis fn-th function from the */
/* block D, get the coefficients representing the value at the curve and the */
/* cross derivatives of two basis function patches adjacent to the curve; */
/* the i-th in the arrays pv, pvv and (i-1) mod hole_k in the arrays pu, puu. */
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int                nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int                nk, m1, i, j, lpkn, lpvkn;
  double             *basis_d;

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  lpkn = sprivate->lastomcknot;
  lpvkn = sprivate->lastpvknot;

  fn -= nfunc_a+nfunc_b+nfunc_c;
  if ( fn < 0 || fn >= nfunc_d )
    return false;

  memset ( fcomc, 0, (lpkn-G1H_OMCDEG)*sizeof(double) );
  if ( cn == fn/(2*nk*m1) )  {
    i = fn % (2*nk*m1);
    if ( (i % 2) == 0 ) fcomc[i/2+3] = 1.0;
    basis_d = sprivate->basis_d;
    j = fn*2*(lpvkn-G1_CROSS01DEG);
    memcpy ( pv, &basis_d[j], (lpvkn-G1_CROSS01DEG)*sizeof(double) );
    j += lpvkn-G1_CROSS01DEG;
    memcpy ( pu, &basis_d[j], (lpvkn-G1_CROSS01DEG)*sizeof(double) );
  }
  else {
    memset ( pv,  0, (lpvkn-G1_CROSS01DEG)*sizeof(double) );
    memset ( pu,  0, (lpvkn-G1_CROSS01DEG)*sizeof(double) );
  }

  return true;
} /*_g1h_GetSplDBasisCrossDerd*/

static boolean _g1h_FindSplBasisFunctionsDd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int    hole_k, nk, m1, i, j, k, l, fn, fnd, nzc;
  int    nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  int    lomckn, lpvkn;
  double *knpv, *knb, *knr, *fcomc, *fcomcd, *basis_d;
  double *b01, *c01, *f01, *g01;
  double *pv;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  lomckn = sprivate->lastomcknot;
  lpvkn = sprivate->lastpvknot;
  knb = sprivate->bezknots;
  knr = sprivate->omcknots;
  knpv = sprivate->pvknots;

        /* allocate memory for the basis function patches cross derivatives */
  basis_d = sprivate->basis_d = malloc ( nfunc_d*2*
                                         (lpvkn-G1_CROSS01DEG)*sizeof(double) );
  if ( !basis_d )
    goto failure;
        /* allocate memory for the spline auxiliary basis function patches */
  fcomc = pkv_GetScratchMemd ( 2*(G1H_OMCDEG+1+nk*m1) );
  if ( !fcomc )
    goto failure;
  fcomcd = &fcomc[G1H_OMCDEG+1+nk*m1];

        /* now setup the basis function patches cross derivatives */
  G1GetPolyAddr0 ( privateG1->jfunc, b01, c01, f01, g01 );
  for ( i = fnd = 0, fn = nfunc_a+nfunc_b+nfunc_c;  i < hole_k;  i++ ) {
    l = (i+hole_k-1) % hole_k;
    for ( j = 0; j < nk*m1; j++ )
      for ( k = 0;  k < 2;  k++, fnd++, fn++ ) {
        _g1h_GetSplDBasisAuxpd ( domain, fn, i, &nzc, fcomc, fcomcd );
        pv = &basis_d[fnd*2*(lpvkn-G1_CROSS01DEG)];
        if ( !_g1h_FindSplCrossDerDd ( knb,
                   &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
                   nzc, nk+1, lomckn,
                   knr, fcomc, fcomcd,
                   G1_CROSS01DEG, lpvkn, knpv, pv ) )
          goto failure;
        pv = &pv[lpvkn-G1_CROSS01DEG];
        if ( !_g1h_FindSplCrossDerDd ( knb,
                   &f01[l*(G1_BF01DEG+1)], &g01[l*(G1_CG01DEG+1)],
                   nzc, nk+1, lomckn,
                   knr, fcomc, fcomcd,
                   G1_CROSS01DEG, lpvkn, knpv, pv ) )
          goto failure;
      }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_FindSplBasisFunctionsDd*/

static boolean _g1h_FindSplBasisFunctionsd ( GHoleDomaind *domain )
{
  void               *sp;
  G1HoleSPrivateRecd *sprivate;
  int                nk, m1, m2, hole_k;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  sprivate = domain->SprivateG1;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
        /* compute the numbers of spline basis functions */
  sprivate->csize = (G1H_FINALDEG-3+nk*m2)*(G1H_FINALDEG-3+nk*m2);
  sprivate->nsfunc_c = hole_k*sprivate->csize;
  sprivate->dsize = 2*nk*m1;
  sprivate->nsfunc_d = hole_k*sprivate->dsize;
        /* construct the knot sequence for the representation of the */
        /* final patches */

  if ( !_g1h_ConstructSplKnotSequencesd ( domain ) )
    goto failure;
        /* find the basis function patches for the block D */
  if ( !_g1h_FindSplBasisFunctionsDd ( domain ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_FindSplBasisFunctionsd*/

boolean g1h_ComputeSplBasisd ( GHoleDomaind *domain, int nk, int m1, int m2 )
{
  if ( nk < 1 || nk > G1H_S_MAX_NK ||
       m1 < 1 || m1 > G1H_S_MAX_M1 ||
       m2 < 1 || m2 > G1H_S_MAX_M2 )
    return false;

  if ( !g1h_AllocSPrivateDatad ( domain, nk, m1, m2 ) )
    return false;

        /* the first step is the construction of the polynomial basis, */
        /* as the junction functions and domain patches will be the same */
  if ( !g1h_ComputeBasisd ( domain ) )
    return false;
  if ( !_g1h_FindSplBasisFunctionsd ( domain ) )
    return false;
  return true;
} /*g1h_ComputeSplBasisd*/

