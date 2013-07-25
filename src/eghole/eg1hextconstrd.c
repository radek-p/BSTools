
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g1h_ExtV0SpaceDimd ( GHoleDomaind *domain )
{
  G1HolePrivateRecd *privateG1;

  privateG1 = domain->privateG1;
  if ( !privateG1 )
    return 0;
  else
    return privateG1->nfunc_a + domain->hole_k*G1_DBDIM;
} /*g1h_ExtV0SpaceDimd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_SetExtConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *ECmat, *ERCmat, *aa, s, t;
  double *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int    hole_k, nfunc_a, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( privateG1->ECmat ) free ( privateG1->ECmat );  /* forget previous constraints */
  if ( privateG1->ERCmat ) free ( privateG1->ERCmat );
  if ( privateG1->Q2ERCMat ) {
    free ( privateG1->Q2ERCMat );
    privateG1->Q2ERCMat = NULL;
  }

  if ( !privateG1->ELmat ) {
    if ( !g1h_DecomposeExtMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG1->nfunc_a;
  nbf = nfunc_a+G1_DBDIM*hole_k;
  if ( nconstr <= 0 || nconstr > nfunc_a )
    goto failure;

  if ( !_g1h_GetExtBlockAddressesd ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  ECmat = privateG1->ECmat = malloc ( nbf*nconstr*sizeof(double) );
  ERCmat = privateG1->ERCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
  if ( !ECmat || !ERCmat )
    goto failure;

  memcpy ( ECmat, cmat, nbf*nconstr*sizeof(double) );
  ECmat = pkv_GetScratchMemd ( (nbf+2)*nconstr );
  if ( !ECmat )
    goto failure;

  aa = &ECmat[nbf*nconstr];
  pkv_TransposeMatrixd ( nconstr, nbf, nbf, privateG1->ECmat, nconstr, ECmat );
  pkn_Block1LowerTrMSolved ( hole_k, G1_DBDIM, nfunc_a, Lii,
                             nconstr, nconstr, ECmat );
  if ( !pkn_QRDecomposeMatrixd ( nbf, nconstr, ECmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      ERCmat[pkn_SymMatIndex(i,j)] = ECmat[i*nconstr+j];
  for ( i = 0; i < nconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = ERCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = ERCmat[pkn_SymMatIndex(i,i)]; 
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG1->extnconstr = nconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG1->ECmat ) free ( privateG1->ECmat );
  if ( privateG1->ERCmat ) free ( privateG1->ERCmat );
  privateG1->ECmat = privateG1->ERCmat = NULL;
  privateG1->extnconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetExtConstraintMatrixd*/

boolean g1h_SetExtAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                         int naconstr, const double *acmat )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *AECmat, *AERCmat, *aa, s, t;
  double *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int    hole_k, nfunc_a, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( privateG1->AECmat ) free ( privateG1->AECmat );  /* forget previous constraints */
  if ( privateG1->AERCmat ) free ( privateG1->AERCmat );
  if ( privateG1->Q2AERCMat ) {
    free ( privateG1->Q2AERCMat );
    privateG1->Q2AERCMat = NULL;
  }

  if ( !privateG1->ELmat ) {
    if ( !g1h_DecomposeExtMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG1->nfunc_a;
  nbf = nfunc_a+G1_DBDIM*hole_k;
  if ( naconstr <= 0 || naconstr > nfunc_a )
    goto failure;

  if ( !_g1h_GetExtBlockAddressesd ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  AECmat = privateG1->AECmat = malloc ( nbf*spdimen*naconstr*sizeof(double) );
  AERCmat = privateG1->AERCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
  if ( !AECmat || !AERCmat )
    goto failure;

  memcpy ( AECmat, acmat, nbf*spdimen*naconstr*sizeof(double) );
  AECmat = pkv_GetScratchMemd ( (nbf*spdimen+2)*naconstr );
  if ( !AECmat )
    goto failure;

  aa = &AECmat[nbf*spdimen*naconstr];
  pkv_TransposeMatrixd ( naconstr, nbf*spdimen, nbf*spdimen,
                         privateG1->AECmat, naconstr, AECmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_Block1LowerTrMSolved ( hole_k, G1_DBDIM, nfunc_a, Lii,
                               naconstr, naconstr, &AECmat[i*nbf*naconstr] );
  if ( !pkn_QRDecomposeMatrixd ( nbf*spdimen, naconstr, AECmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < naconstr; i++ )
    for ( j = i; j < naconstr; j++ )
      AERCmat[pkn_SymMatIndex(i,j)] = AECmat[i*naconstr+j];
  for ( i = 0; i < naconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = AERCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = AERCmat[pkn_SymMatIndex(i,i)]; 
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG1->extacdim = spdimen;
  privateG1->extnaconstr = naconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG1->AECmat ) free ( privateG1->AECmat );
  if ( privateG1->AERCmat ) free ( privateG1->AERCmat );
  privateG1->AECmat = privateG1->AERCmat = NULL;
  privateG1->extnaconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetExtAltConstraintMatrixd*/

