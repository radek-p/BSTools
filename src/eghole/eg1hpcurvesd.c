
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


boolean g1h_GetFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                   CONST_ double *hole_cp, double *acoeff,
                                   void (*outcurve) ( int n, const double *cp ) )
{
  void          *sp;
  int           hole_k;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  int           nfunc_a, nfunc_b;
  double        *bcp, *cp, *y;
  double        *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn;
  bcp = pkv_GetScratchMemd ( hole_k*spdimen*(G1_CROSS00DEG+1) );
  y = pkv_GetScratchMemd ( spdimen*(G1_CROSS00DEG+1) );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*(G1_CROSS00DEG+1)*sizeof(double) );
  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFBPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc10, &bc11,
                                &bd00, &bd01, &bd10, &bd11 );
      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &bcp[k*(G1_CROSS00DEG+1)] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );
      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*(G1_CROSS00DEG+1)],
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
} /*g1h_GetFinalPatchCurvesd*/

boolean g1h_GetExtFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, const double *cp ) )
{
  int hole_k, nfunc_c;

  hole_k = domain->hole_k;
  nfunc_c = hole_k*G1_DBDIM;
  return g1h_GetFinalPatchCurvesd ( domain, spdimen, hole_cp,
                                    &acoeff[spdimen*nfunc_c], outcurve );
} /*g1h_GetExtFinalPatchCurvesd*/

boolean g1h_GetABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc )
{
  G1HolePrivateRecd *privateG1;
  int   hole_k, nfunc_a;
  double *bc00, *bc01, *bd00, *bd01;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    return false;
  _g1h_GetBFAPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bd00, &bd01 );
  memcpy ( bfpc, bc00, (G1_CROSS00DEG+1)*sizeof(double) );
  return true;
} /*g1h_GetABasisFPatchCurved*/

boolean g1h_GetBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc )
{
  int   hole_k, nfunc_b;
  double *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    return false;
  _g1h_GetBFBPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bc10, &bc11,
                            &bd00, &bd01, &bd10, &bd11 );
  memcpy ( bfpc, bc00, (G1_CROSS00DEG+1)*sizeof(double) );
  return true;
} /*g1h_GetBBasisFPatchCurved*/

