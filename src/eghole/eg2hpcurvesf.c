
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


boolean g2h_GetFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                   CONST_ float *hole_cp, float *acoeff,
                                   void (*outcurve) ( int n, const float *cp ) )
{
  void          *sp;
  int           hole_k;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  int           nfunc_a, nfunc_b;
  float         *bcp, *cp, *y;
  float         *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
                *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn;
  bcp = pkv_GetScratchMemf ( hole_k*spdimen*(G2_CROSS00DEG+1) );
  y = pkv_GetScratchMemf ( spdimen*(G2_CROSS00DEG+1) );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*(G2_CROSS00DEG+1)*sizeof(float) );
  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12, 
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &bcp[k*(G2_CROSS00DEG+1)] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &bcp[k*(G2_CROSS00DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k+= spdimen )
    outcurve ( G2_CROSS00DEG, &bcp[k*(G2_CROSS00DEG+1)] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetFinalPatchCurvesf*/

boolean g2h_GetExtFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, const float *cp ) )
{
  int hole_k, nfunc_c;

  hole_k = domain->hole_k;
  nfunc_c = hole_k*G2_DBDIM;
  return g2h_GetFinalPatchCurvesf ( domain, spdimen, hole_cp,
                                    &acoeff[spdimen*nfunc_c], outcurve );
} /*g2h_GetExtFinalPatchCurvesf*/

boolean g2h_GetABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc )
{
  G2HolePrivateRecf *privateG2;
  int   hole_k, nfunc_a;
  float *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    return false;
  _g2h_GetBFAPatchCurvesf ( domain, fn, i, &bc00, &bc02, &bc01,
                            &bd00, &bd01, &bd02 );
  memcpy ( bfpc, bc00, (G2_CROSS00DEG+1)*sizeof(float) );
  return true;
} /*g2h_GetABasisFPatchCurvef*/

boolean g2h_GetBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc )
{
  int    hole_k, nfunc_b;
  float *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
        *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    return false;
  _g2h_GetBFBPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                            &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
  memcpy ( bfpc, bc00, (G2_CROSS00DEG+1)*sizeof(float) );
  return true;
} /*g2h_GetBBasisFPatchCurvef*/

