
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_Q2DrawMatricesf ( GHoleDomainf *domain,
                           void (*drawmatrix)(int nfa, int nfb,
                                              float *amat, float *bmat) )
{
  void                *sp;
  float               *amat, *bmat;
  int                 na, nb;
  G1HolePrivateRecf   *privateG1;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1 )
    goto wayout;
  if ( !privateG1->Q2AMat || !privateG1->Q2BMat )
    goto wayout;

  na = privateG1->nfunc_a;
  na = na*(na+1)/2;
  nb = privateG1->nfunc_a*privateG1->nfunc_b;
  amat = pkv_GetScratchMemf ( na );
  bmat = pkv_GetScratchMemf ( nb );
  if ( !amat || !bmat )
    goto wayout;

  memcpy ( amat, privateG1->Q2AMat, na*sizeof(float) );
  memcpy ( bmat, privateG1->Q2BMat, nb*sizeof(float) );
  drawmatrix ( privateG1->nfunc_a, privateG1->nfunc_b, amat, bmat );
 
wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawMatricesf*/

void g1h_Q2DrawExtMatricesf ( GHoleDomainf *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 float *Aii, float *Bi) )
{
  void  *sp;
  G1HolePrivateRecf   *privateG1;
  float *amat, *bmat;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, Asize, Bsize;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1 )
    goto wayout;
  if ( !privateG1->Q2EAMat || !privateG1->Q2EBMat )
    goto wayout;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = G1_DBDIM*hole_k;
  Asize = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
  Bsize = (nfunc_a+nfunc_c)*nfunc_b;
  amat = pkv_GetScratchMemf ( Asize );
  bmat = pkv_GetScratchMemf ( Bsize );
  if ( !amat || !bmat )
    goto wayout;
  memcpy ( amat, privateG1->Q2EAMat, Asize*sizeof(float) );
  memcpy ( bmat, privateG1->Q2EBMat, Bsize*sizeof(float) );
  drawmatrix ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawExtMatricesf*/

void g1h_Q2DrawSplMatricesf ( GHoleDomainf *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 float *Aii, float *Bi) )
{
  void *sp;
  G1HolePrivateRecf   *privateG1;
  G1HoleSPrivateRecf  *privateS;
  float *amat, *bmat;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, asize, bsize;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  privateS  = domain->SprivateG1;
  if ( !privateG1 || !privateS )
    goto wayout;
  if ( !privateS->Q2SAMat || !privateS->Q2SBMat )
    goto wayout;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = privateS->nsfunc_c;
  nfunc_d = privateS->nsfunc_d;
  asize = pkn_Block3ArraySize ( hole_k-1, nfunc_c/hole_k,
                  nfunc_c/hole_k+nfunc_d+nfunc_a );
  bsize = nfunc_b*(nfunc_a+nfunc_c+nfunc_d);
  amat = pkv_GetScratchMemf ( asize );
  bmat = pkv_GetScratchMemf ( bsize );
  if ( !amat || !bmat )
    goto wayout;
  memcpy ( amat, privateS->Q2SAMat, asize*sizeof(float) );
  memcpy ( bmat, privateS->Q2SBMat, bsize*sizeof(float) );
  drawmatrix ( hole_k-1, nfunc_c/hole_k, nfunc_c/hole_k+nfunc_d+nfunc_a,
               amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawSplMatricesf*/

