
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    float *patchmatrix;
    int   i;
  } _g1h_auxstr;

static void _g1h_outsympatchmatrixf ( int n, int m, const float *cp, void *usrptr )
{
  _g1h_auxstr *auxs;
  float       *patchmatrix;
  int         nrows;

  auxs = (_g1h_auxstr*)usrptr;
  nrows = (n+1)*(m+1);  /* (9+1)*(9+1) == 100 */
  if ( auxs->i ) {
    patchmatrix = &auxs->patchmatrix[nrows*(6*auxs->i+1)];
    pkv_TransposeMatrixf ( nrows, 6, 7, &cp[1], nrows, patchmatrix );
  }
  else
    pkv_TransposeMatrixf ( nrows, 7, 7, cp, nrows, auxs->patchmatrix );
  auxs->i ++;
} /*_g1h_outsympatchmatrixf*/

boolean g1h_GetSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix )
{
  void        *sp;
  G1HolePrivateRecf *privateG1;
  int         hole_k;
  float       *cp;
  _g1h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Lmat )
    if ( !g1h_DecomposeMatrixf ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemf ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(float) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g1h_FillHolef ( domain, 7, cp, NULL, &auxs,
                        _g1h_outsympatchmatrixf ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetSymPatchMatrixf*/

boolean g1h_GetExtSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix )
{
  void        *sp;
  G1HolePrivateRecf *privateG1;
  int         hole_k;
  float       *cp;
  _g1h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->ELmat )
    if ( !g1h_DecomposeExtMatrixf ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemf ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(float) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g1h_ExtFillHolef ( domain, 7, cp, NULL, &auxs,
                           _g1h_outsympatchmatrixf ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetExtSymPatchMatrixf*/

boolean g1h_Q2GetSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix )
{
  void        *sp;
  G1HolePrivateRecf *privateG1;
  int         hole_k;
  float       *cp;
  _g1h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixf ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemf ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(float) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g1h_Q2FillHolef ( domain, 7, cp, NULL, &auxs,
                          _g1h_outsympatchmatrixf ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2GetSymPatchMatrixf*/

boolean g1h_Q2GetExtSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix )
{
  void        *sp;
  G1HolePrivateRecf *privateG1;
  int         hole_k;
  float       *cp;
  _g1h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2ELMat )
    if ( !g1h_Q2ExtDecomposeMatrixf ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemf ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(float) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g1h_Q2ExtFillHolef ( domain, 7, cp, NULL, &auxs,
                             _g1h_outsympatchmatrixf ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2GetExtSymPatchMatrixf*/

boolean g1h_MatrixFillSymHolef ( int hole_k, const float *patchmatrix,
                                 int spdimen, const float *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const float *cp,
                                                    void *usrptr ) )
{
  void  *sp;
  int   i, j, l, m, n, ncp;
  float *cp;

  sp = pkv_GetScratchMemTop ();
  ncp = (G1H_FINALDEG+1)*(G1H_FINALDEG+1);
  cp = pkv_GetScratchMemf ( spdimen*ncp );
  if ( !cp )
    goto failure;
  for ( i = 0; i < hole_k; i++ ) {
    memset ( cp, 0, spdimen*ncp*sizeof(float) );
    pkn_MultMatrixf ( ncp, 1, 1, patchmatrix, spdimen, spdimen, hole_cp,
                      spdimen, cp );
    for ( j = 0; j < hole_k; j++ ) {
      l = (i-j+hole_k) % hole_k;
      m = 12*l+1;
      n = 6*j;
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+1)*ncp],
                           spdimen, spdimen, &hole_cp[m*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+2)*ncp],
                           spdimen, spdimen, &hole_cp[(m+1)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+3)*ncp],
                           spdimen, spdimen, &hole_cp[(m+3)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+4)*ncp],
                           spdimen, spdimen, &hole_cp[(m+4)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+5)*ncp],
                           spdimen, spdimen, &hole_cp[(m+6)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddf ( ncp, 1, 1, &patchmatrix[(n+6)*ncp],
                           spdimen, spdimen, &hole_cp[(m+7)*spdimen],
                           spdimen, cp );
    }
    outpatch ( G1H_FINALDEG, G1H_FINALDEG, cp, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_MatrixFillSymHolef*/

