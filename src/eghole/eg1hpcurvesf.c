
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


boolean g1h_GetFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                   CONST_ float *hole_cp, float *acoeff,
                                   void (*outcurve) ( int n, const float *cp ) )
{
  void          *sp;
  int           hole_k;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  int           nfunc_a, nfunc_b;
  float         *bcp, *cp, *y;
  float         *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn;
  bcp = pkv_GetScratchMemf ( hole_k*spdimen*(G1_CROSS00DEG+1) );
  y = pkv_GetScratchMemf ( spdimen*(G1_CROSS00DEG+1) );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*(G1_CROSS00DEG+1)*sizeof(float) );
  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFBPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc10, &bc11,
                                &bd00, &bd01, &bd10, &bd11 );
      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &bcp[k*(G1_CROSS00DEG+1)] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesf ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );
      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &bcp[k*(G1_CROSS00DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k+= spdimen )
    outcurve ( G1_CROSS00DEG, &bcp[k*(G1_CROSS00DEG+1)] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetFinalPatchCurvesf*/

boolean g1h_GetExtFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, const float *cp ) )
{
  int hole_k, nfunc_c;

  hole_k = domain->hole_k;
  nfunc_c = hole_k*G1_DBDIM;
  return g1h_GetFinalPatchCurvesf ( domain, spdimen, hole_cp,
                                    &acoeff[spdimen*nfunc_c], outcurve );
} /*g1h_GetExtFinalPatchCurvesf*/

boolean g1h_GetABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc )
{
  G1HolePrivateRecf *privateG1;
  int   hole_k, nfunc_a;
  float *bc00, *bc01, *bd00, *bd01;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    return false;
  _g1h_GetBFAPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bd00, &bd01 );
  memcpy ( bfpc, bc00, (G1_CROSS00DEG+1)*sizeof(float) );
  return true;
} /*g1h_GetABasisFPatchCurvef*/

boolean g1h_GetBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc )
{
  int   hole_k, nfunc_b;
  float *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    return false;
  _g1h_GetBFBPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bc10, &bc11,
                            &bd00, &bd01, &bd10, &bd11 );
  memcpy ( bfpc, bc00, (G1_CROSS00DEG+1)*sizeof(float) );
  return true;
} /*g1h_GetBBasisFPatchCurvef*/

