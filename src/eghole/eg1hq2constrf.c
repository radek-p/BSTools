
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
boolean g1h_Q2FillHoleConstrf ( GHoleDomainf *domain,
                                int spdimen, CONST_ float *hole_cp,
                                int nconstr, CONST_ float *constr,
                                float *acoeff, void *usrptr,
                                void (*outpatch) ( int n, int m, const float *cp,
                                                   void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecf *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  float  *lmat, *cmat, *rcmat, *aa, *b, *x, *y;
  float  *fc00;
  float  s, t;
  int    i, j;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Cmat || privateG1->nconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixf ( domain ) )
      goto failure;
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  if ( !(rcmat = privateG1->Q2RCMat) ) {
    rcmat = privateG1->Q2RCMat =
              (float*)malloc ( (nconstr*(nconstr+1))/2*sizeof(float) );
    if ( !rcmat )
      goto failure;
    cmat = pkv_GetScratchMemf ( (nfunc_a+2)*nconstr );
    if ( !cmat )
      goto failure;
    aa = &cmat[nfunc_a*nconstr];
    pkv_TransposeMatrixf ( nconstr, nfunc_a, nfunc_a, privateG1->Cmat,
                           nconstr, cmat );
    pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, nconstr, nconstr, cmat,
                             nconstr, cmat  );
    if ( !pkn_QRDecomposeMatrixf ( nfunc_a, nconstr, cmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < nconstr; i++ )
      for ( j = i; j < nconstr; j++ )
        rcmat[pkn_SymMatIndex(i,j)] = cmat[i*nconstr+j];
    for ( i = 0; i < nconstr; i++ ) {
      for ( j = 0, s = 0.0;  j < i; j++ ) {
        t = rcmat[pkn_SymMatIndex(i,j)];
        s += t*t;
      }
      t = rcmat[pkn_SymMatIndex(i,i)];
      if ( t*t < 1.0e-10*s ) {
        domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
        goto failure;
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  cmat = privateG1->Cmat;

  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSidef ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixf ( nconstr, nfunc_a, nfunc_a, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixf ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemf ( spdimen*nfunc_a );
  if ( !y ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixf ( nconstr, nfunc_a, nfunc_a, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemf ( spdimen*nfunc_a );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g1h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHoleConstrf*/

boolean g1h_Q2FillHoleAltConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int naconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecf *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  float  *lmat, *acmat, *arcmat, *aa, *b, *x, *y;
  float  *fc00;
  float  s, t;
  int    i, j;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->ACmat ||
       privateG1->naconstr != naconstr || privateG1->acdim != spdimen ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixf ( domain ) )
      goto failure;
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  if ( !(arcmat = privateG1->Q2ARCMat) ) {
    arcmat = privateG1->Q2ARCMat =
              (float*)malloc ( (naconstr*(naconstr+1))/2*sizeof(float) );
    if ( !arcmat )
      goto failure;
    acmat = pkv_GetScratchMemf ( (nfunc_a*spdimen+2)*naconstr );
    if ( !acmat )
      goto failure;
    aa = &acmat[nfunc_a*spdimen*naconstr];
    pkv_TransposeMatrixf ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a,
                           privateG1->ACmat, naconstr, acmat );
    for ( i = 0; i < spdimen; i++ )
      pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, naconstr,
                                naconstr, &acmat[i*nfunc_a*naconstr],
                                naconstr, &acmat[i*nfunc_a*naconstr] );
    if ( !pkn_QRDecomposeMatrixf ( nfunc_a*spdimen, naconstr, acmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < naconstr; i++ )
      for ( j = i; j < naconstr; j++ )
        arcmat[pkn_SymMatIndex(i,j)] = acmat[i*naconstr+j];
    for ( i = 0; i < naconstr; i++ ) {
      for ( j = i, s = 0.0;  j < naconstr;  j++  ) {
        t = arcmat[pkn_SymMatIndex(i,j)];
        s += t*t;
      }
      t = arcmat[pkn_SymMatIndex(i,i)];
      if ( t*t < 1.0e-10*s ) {
        domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
        goto failure;
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  acmat = privateG1->ACmat;

  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  y = pkv_GetScratchMemf ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSidef ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixf ( nfunc_a, spdimen, spdimen, x, nfunc_a, y );
  pkn_MultMatrixf ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixf ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixf ( naconstr, nfunc_a*spdimen, nfunc_a*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixf ( spdimen, nfunc_a, nfunc_a, y, spdimen, b );
  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, y );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g1h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHoleAltConstrf*/

