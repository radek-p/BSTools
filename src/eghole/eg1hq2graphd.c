
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
void g1h_Q2DrawMatricesd ( GHoleDomaind *domain,
                           void (*drawmatrix)(int nfa, int nfb,
                                              double *amat, double *bmat) )
{
  void                *sp;
  double              *amat, *bmat;
  int                 na, nb;
  G1HolePrivateRecd   *privateG1;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1 )
    goto wayout;
  if ( !privateG1->Q2AMat || !privateG1->Q2BMat )
    goto wayout;

  na = privateG1->nfunc_a;
  na = na*(na+1)/2;
  nb = privateG1->nfunc_a*privateG1->nfunc_b;
  amat = pkv_GetScratchMemd ( na );
  bmat = pkv_GetScratchMemd ( nb );
  if ( !amat || !bmat )
    goto wayout;

  memcpy ( amat, privateG1->Q2AMat, na*sizeof(double) );
  memcpy ( bmat, privateG1->Q2BMat, nb*sizeof(double) );
  drawmatrix ( privateG1->nfunc_a, privateG1->nfunc_b, amat, bmat );
 
wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawMatricesd*/

void g1h_Q2DrawExtMatricesd ( GHoleDomaind *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 double *Aii, double *Bi) )
{
  void   *sp;
  G1HolePrivateRecd   *privateG1;
  double *amat, *bmat;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, Asize, Bsize;

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
  amat = pkv_GetScratchMemd ( Asize );
  bmat = pkv_GetScratchMemd ( Bsize );
  if ( !amat || !bmat )
    goto wayout;
  memcpy ( amat, privateG1->Q2EAMat, Asize*sizeof(double) );
  memcpy ( bmat, privateG1->Q2EBMat, Bsize*sizeof(double) );
  drawmatrix ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawExtMatricesd*/

void g1h_Q2DrawSplMatricesd ( GHoleDomaind *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 double *Aii, double *Bi) )
{
  void *sp;
  G1HolePrivateRecd   *privateG1;
  G1HoleSPrivateRecd  *privateS;
  double *amat, *bmat;
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
  amat = pkv_GetScratchMemd ( asize );
  bmat = pkv_GetScratchMemd ( bsize );
  if ( !amat || !bmat )
    goto wayout;
  memcpy ( amat, privateS->Q2SAMat, asize*sizeof(double) );
  memcpy ( bmat, privateS->Q2SBMat, bsize*sizeof(double) );
  drawmatrix ( hole_k-1, nfunc_c/hole_k, nfunc_c/hole_k+nfunc_d+nfunc_a,
               amat, bmat );

wayout:
  pkv_SetScratchMemTop ( sp );
} /*g1h_Q2DrawSplMatricesd*/

