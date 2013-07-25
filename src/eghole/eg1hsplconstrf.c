
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g1h_SplV0SpaceDimf ( GHoleDomainf *domain )
{
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate; 

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( sprivate )
    return privateG1->nfunc_a + sprivate->nsfunc_c + sprivate->nsfunc_d;
  else
    return -1;
} /*g1h_SplV0SpaceDimf*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_SetSplConstraintMatrixf ( GHoleDomainf *domain,
                                      int nconstr, const float *cmat )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  G1HoleSPrivateRecf *sprivate;
  float *SCmat, *SRCmat, *aa;
  float s, t;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( nconstr <= 0 || nconstr > nfunc_a+nfunc_d )
    goto failure;

  if ( sprivate->SCmat ) free ( sprivate->SCmat );  /* forget previous constraints */
  if ( sprivate->SRCmat ) free ( sprivate->SRCmat );
  if ( sprivate->Q2SRCmat ) {
    free ( sprivate->Q2SRCmat );
    sprivate->Q2SRCmat = NULL;
  }

  if ( !sprivate->SLMat ) {
    if ( !g1h_DecomposeSplMatrixf ( domain ) )
      goto failure;
  }

  SCmat = sprivate->SCmat = malloc ( nbf*nconstr*sizeof(float) );
  SRCmat = sprivate->SRCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(float) );
  if ( !SCmat || !SRCmat )
    goto failure;

  memcpy ( SCmat, cmat, nbf*nconstr*sizeof(float) );
  SCmat = pkv_GetScratchMemf ( (nbf+2)*nconstr );
  if ( !SCmat )
    goto failure;

  aa = &SCmat[nbf*nconstr];
  pkv_TransposeMatrixf ( nconstr, nbf, nbf, sprivate->SCmat, nconstr, SCmat );
  pkn_Block2LowerTrMSolvef ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             sprivate->SLMat, nconstr, nconstr, SCmat );
  if ( !pkn_QRDecomposeMatrixf ( nbf, nconstr, SCmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      SRCmat[pkn_SymMatIndex(i,j)] = SCmat[i*nconstr+j];
  for ( i = 0; i < nconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = SRCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = SRCmat[pkn_SymMatIndex(i,i)];
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  sprivate->splnconstr = nconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( sprivate->SCmat ) free ( sprivate->SCmat );
  if ( sprivate->SRCmat ) free ( sprivate->SRCmat );
  sprivate->SCmat = sprivate->SRCmat = NULL;
  sprivate->splnconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetSplConstraintMatrixf*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_SetSplAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                         int naconstr, const float *acmat )
{
  void *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  float *ASCmat, *ASRCmat, *aa, s, t;
  float *SLmat;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( naconstr <= 0 || naconstr > nfunc_a+nfunc_d )
    goto failure;

  if ( sprivate->ASCmat ) free ( sprivate->ASCmat );  /* forget previous constraints */
  if ( sprivate->ASRCmat ) free ( sprivate->ASRCmat );
  if ( sprivate->Q2SARCmat ) {
    free ( sprivate->Q2SARCmat );
    sprivate->Q2SARCmat = NULL;
  }

  if ( !sprivate->SLMat ) {
    if ( !g1h_DecomposeSplMatrixf ( domain ) )
      goto failure;
  }

  SLmat = sprivate->SLMat;
  ASCmat = sprivate->ASCmat = malloc ( nbf*spdimen*naconstr*sizeof(float) );
  ASRCmat = sprivate->ASRCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(float) );
  if ( !ASCmat || !ASRCmat )
    goto failure;

  memcpy ( ASCmat, acmat, nbf*spdimen*naconstr*sizeof(float) );
  ASCmat = pkv_GetScratchMemf ( (nbf*spdimen+2)*naconstr );
  if ( !ASCmat )
    goto failure;

  aa = &ASCmat[nbf*spdimen*naconstr];
  pkv_TransposeMatrixf ( naconstr, nbf*spdimen, nbf*spdimen,
                         sprivate->ASCmat, naconstr, ASCmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_Block2LowerTrMSolvef ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                               SLmat, naconstr, naconstr, &ASCmat[i*nbf*naconstr] );
  if ( !pkn_QRDecomposeMatrixf ( nbf*spdimen, naconstr, ASCmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < naconstr; i++ )
    for ( j = i; j < naconstr; j++ )
      ASRCmat[pkn_SymMatIndex(i,j)] = ASCmat[i*naconstr+j];
  for ( i = 0; i < naconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = ASRCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = ASRCmat[pkn_SymMatIndex(i,i)];
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  sprivate->splacdim = spdimen;
  sprivate->splnaconstr = naconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( sprivate->ASCmat ) free ( sprivate->ASCmat );
  if ( sprivate->ASRCmat ) free ( sprivate->ASRCmat );
  sprivate->ASCmat = sprivate->ASRCmat = NULL;
  sprivate->splnaconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetSplAltConstraintMatrixf*/

