
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
#include "eg2holed.h"

#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g2h_ExtV0SpaceDimd ( GHoleDomaind *domain )
{
  G2HolePrivateRecd *privateG2;

  privateG2 = domain->privateG2;
  if ( !privateG2 )
    return 0;
  else
    return privateG2->nfunc_a + domain->hole_k*G2_DBDIM;
} /*g2h_ExtV0SpaceDimd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2h_SetExtConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double  *ECmat, *ERCmat, *aa, s, t;
  double  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int    hole_k, nfunc_a, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( privateG2->ECmat ) free ( privateG2->ECmat );  /* forget previous constraints */
  if ( privateG2->ERCmat ) free ( privateG2->ERCmat );

  if ( !privateG2->ELmat ) {
    if ( !g2h_DecomposeExtMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG2->nfunc_a;
  nbf = nfunc_a+G2_DBDIM*hole_k;
  if ( nconstr <= 0 || nconstr > nfunc_a )
    goto failure;

  if ( !_g2h_GetExtBlockAddressesd ( domain, &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  ECmat = privateG2->ECmat = malloc ( nbf*nconstr*sizeof(double) );
  ERCmat = privateG2->ERCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
  if ( !ECmat || !ERCmat )
    goto failure;

  memcpy ( ECmat, cmat, nbf*nconstr*sizeof(double) );
  ECmat = pkv_GetScratchMemd ( (nbf+2)*nconstr );
  if ( !ECmat )
    goto failure;

  aa = &ECmat[nbf*nconstr];
  pkv_TransposeMatrixd ( nconstr, nbf, nbf, privateG2->ECmat, nconstr, ECmat );
  pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a, Lii,
                             nconstr, nconstr, ECmat );
  if ( !pkn_QRDecomposeMatrixd ( nbf, nconstr, ECmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG2->extnconstr = nconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG2->ECmat ) free ( privateG2->ECmat );
  if ( privateG2->ERCmat ) free ( privateG2->ERCmat );
  privateG2->ECmat = privateG2->ERCmat = NULL;
  privateG2->extnconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SetExtConstraintMatrixd*/

boolean g2h_SetExtAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                         int naconstr, const double *acmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double  *AECmat, *AERCmat, *aa, s, t;
  double  *Aii, *Aki, *Akk, *Bi, *Bk, *Lii;
  int    hole_k, nfunc_a, nbf, i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( privateG2->AECmat ) free ( privateG2->AECmat );  /* forget previous constraints */
  if ( privateG2->AERCmat ) free ( privateG2->AERCmat );

  if ( !privateG2->ELmat ) {
    if ( !g2h_DecomposeExtMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG2->nfunc_a;
  nbf = nfunc_a+G2_DBDIM*hole_k;
  if ( naconstr <= 0 || naconstr > nfunc_a )
    goto failure;

  if ( !_g2h_GetExtBlockAddressesd ( domain,
                          &Aii, &Aki, &Akk, &Bi, &Bk, &Lii ) )
    goto failure;

  AECmat = privateG2->AECmat = malloc ( nbf*spdimen*naconstr*sizeof(double) );
  AERCmat = privateG2->AERCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
  if ( !AECmat || !AERCmat )
    goto failure;

  memcpy ( AECmat, acmat, nbf*spdimen*naconstr*sizeof(double) );
  AECmat = pkv_GetScratchMemd ( (nbf*spdimen+2)*naconstr );
  if ( !AECmat )
    goto failure;

  aa = &AECmat[nbf*spdimen*naconstr];
  pkv_TransposeMatrixd ( naconstr, nbf*spdimen, nbf*spdimen,
                         privateG2->AECmat, naconstr, AECmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a, Lii,
                               naconstr, naconstr, &AECmat[i*nbf*naconstr] );
  if ( !pkn_QRDecomposeMatrixd ( nbf*spdimen, naconstr, AECmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG2->extacdim = spdimen;
  privateG2->extnaconstr = naconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG2->AECmat ) free ( privateG2->AECmat );
  if ( privateG2->AERCmat ) free ( privateG2->AERCmat );
  privateG2->AECmat = privateG2->AERCmat = NULL;
  privateG2->extnconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SetExtAltConstraintMatrixd*/

