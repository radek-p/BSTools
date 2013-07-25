
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
typedef struct {
    double *patchmatrix;
    int    i;
  } _g2h_auxstr;

static void _g2h_outsympatchmatrixd ( int n, int m, const double *cp, void *usrptr )
{
  _g2h_auxstr *auxs;
  double      *patchmatrix;
  int         nrows;

  auxs = (_g2h_auxstr*)usrptr;
  nrows = (n+1)*(m+1);  /* (9+1)*(9+1) == 100 */
  if ( auxs->i ) {
    patchmatrix = &auxs->patchmatrix[nrows*(6*auxs->i+1)];
    pkv_TransposeMatrixd ( nrows, 6, 7, &cp[1], nrows, patchmatrix );
  }
  else
    pkv_TransposeMatrixd ( nrows, 7, 7, cp, nrows, auxs->patchmatrix );
  auxs->i ++;
} /*_g2h_outsympatchmatrixd*/

boolean g2h_GetSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix )
{
  void        *sp;
  G2HolePrivateRecd *privateG2;
  int         hole_k;
  double       *cp;
  _g2h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->Lmat )
    if ( !g2h_DecomposeMatrixd ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemd ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(double) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g2h_FillHoled ( domain, 7, cp, NULL, &auxs,
                        _g2h_outsympatchmatrixd ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSymPatchMatrixd*/

boolean g2h_GetExtSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix )
{
  void        *sp;
  G2HolePrivateRecd *privateG2;
  int         hole_k;
  double       *cp;
  _g2h_auxstr auxs;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  if ( !privateG2->ELmat )
    if ( !g2h_DecomposeExtMatrixd ( domain ) )
      goto failure;
  cp = pkv_GetScratchMemd ( 7*(12*hole_k+1) );
  if ( !cp )
    goto failure;
  memset ( cp, 0, 7*(12*hole_k+1)*sizeof(double) );
  cp[0] = cp[8] = cp[16] = cp[31] = cp[39] = cp[54] = cp[62] = 1.0;
  auxs.patchmatrix = patchmatrix;
  auxs.i = 0;
  if ( !g2h_ExtFillHoled ( domain, 7, cp, NULL, &auxs,
                           _g2h_outsympatchmatrixd ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetExtSymPatchMatrixd*/

boolean g2h_MatrixFillSymHoled ( int hole_k, const double *patchmatrix,
                                 int spdimen, const double *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const double *cp,
                                                    void *usrptr ) )
{
  void  *sp;
  int   i, j, l, m, n, ncp;
  double *cp;

  sp = pkv_GetScratchMemTop ();
  ncp = (G2H_FINALDEG+1)*(G2H_FINALDEG+1);
  cp = pkv_GetScratchMemd ( spdimen*ncp );
  if ( !cp )
    goto failure;
  for ( i = 0; i < hole_k; i++ ) {
    memset ( cp, 0, spdimen*ncp*sizeof(double) );
    pkn_MultMatrixd ( ncp, 1, 1, patchmatrix, spdimen, spdimen, hole_cp,
                      spdimen, cp );
    for ( j = 0; j < hole_k; j++ ) {
      l = (i-j+hole_k) % hole_k;
      m = 12*l+1;
      n = 6*j;
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+1)*ncp],
                           spdimen, spdimen, &hole_cp[m*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+2)*ncp],
                           spdimen, spdimen, &hole_cp[(m+1)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+3)*ncp],
                           spdimen, spdimen, &hole_cp[(m+3)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+4)*ncp],
                           spdimen, spdimen, &hole_cp[(m+4)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+5)*ncp],
                           spdimen, spdimen, &hole_cp[(m+6)*spdimen],
                           spdimen, cp );
      pkn_MultMatrixAddd ( ncp, 1, 1, &patchmatrix[(n+6)*ncp],
                           spdimen, spdimen, &hole_cp[(m+7)*spdimen],
                           spdimen, cp );
    }
    outpatch ( G2H_FINALDEG, G2H_FINALDEG, cp, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_MatrixFillSymHoled*/

