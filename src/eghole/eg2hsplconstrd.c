
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g2h_SplV0SpaceDimd ( GHoleDomaind *domain )
{
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate; 

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( sprivate )
    return privateG2->nfunc_a + sprivate->nsfunc_c + sprivate->nsfunc_d;
  else
    return -1;
} /*g2h_SplV0SpaceDimd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2h_SetSplConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *SCmat, *SRCmat, *aa;
  double s, t;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( nconstr <= 0 || nconstr > nfunc_a+nfunc_d )
    goto failure;

  if ( sprivate->SCmat ) free ( sprivate->SCmat );  /* forget previous constraints */
  if ( sprivate->SRCmat ) free ( sprivate->SRCmat );

  if ( !sprivate->SLMat ) {
    if ( !g2h_DecomposeSplMatrixd ( domain ) )
      goto failure;
  }

  SCmat = sprivate->SCmat = malloc ( nbf*nconstr*sizeof(double) );
  SRCmat = sprivate->SRCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
  if ( !SCmat || !SRCmat )
    goto failure;

  memcpy ( SCmat, cmat, nbf*nconstr*sizeof(double) );
  SCmat = pkv_GetScratchMemd ( (nbf+2)*nconstr );
  if ( !SCmat )
    goto failure;

  aa = &SCmat[nbf*nconstr];
  pkv_TransposeMatrixd ( nconstr, nbf, nbf, sprivate->SCmat, nconstr, SCmat );
  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             sprivate->SLMat, nconstr, nconstr, SCmat );
  if ( !pkn_QRDecomposeMatrixd ( nbf, nconstr, SCmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
} /*g2h_SetSplConstraintMatrixd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2h_SetSplAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                         int naconstr, const double *acmat )
{
  void *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *ASCmat, *ASRCmat, *aa, s, t;
  double *SLmat;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( naconstr <= 0 || naconstr > nfunc_a+nfunc_d )
    goto failure;

  if ( sprivate->ASCmat ) free ( sprivate->ASCmat );  /* forget previous constraints */
  if ( sprivate->ASRCmat ) free ( sprivate->ASRCmat );

  if ( !sprivate->SLMat ) {
    if ( !g2h_DecomposeSplMatrixd ( domain ) )
      goto failure;
  }

  SLmat = sprivate->SLMat;
  ASCmat = sprivate->ASCmat = malloc ( nbf*spdimen*naconstr*sizeof(double) );
  ASRCmat = sprivate->ASRCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
  if ( !ASCmat || !ASRCmat )
    goto failure;

  memcpy ( ASCmat, acmat, nbf*spdimen*naconstr*sizeof(double) );
  ASCmat = pkv_GetScratchMemd ( (nbf*spdimen+2)*naconstr );
  if ( !ASCmat )
    goto failure;

  aa = &ASCmat[nbf*spdimen*naconstr];
  pkv_TransposeMatrixd ( naconstr, nbf*spdimen, nbf*spdimen,
                         sprivate->ASCmat, naconstr, ASCmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                               SLmat, naconstr, naconstr, &ASCmat[i*nbf*naconstr] );
  if ( !pkn_QRDecomposeMatrixd ( nbf*spdimen, naconstr, ASCmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
} /*g2h_SetSplAltConstraintMatrixd*/

