
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecf *privateG1;
  float  *lmat, *cmat, *rcmat, *aa;
  float  *fc00, *b, *x, *y;
  float  s, t;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  int    i, j;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->ECmat || privateG1->extnconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  if ( !privateG1->Q2EAMat )
    if ( !g1h_Q2ExtComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat )
    if ( !g1h_Q2ExtDecomposeMatrixf ( domain ) )
      goto failure;
  lmat = privateG1->Q2ELMat;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  if ( !(rcmat = privateG1->Q2ERCMat) ) {
    rcmat = privateG1->Q2ERCMat =
              (float*)malloc ( (nconstr*(nconstr+1))/2*sizeof(float) );
    if ( !rcmat )
      goto failure;
    cmat = pkv_GetScratchMemf ( (nbf+2)*nconstr );
    if ( !cmat )
      goto failure;
    aa = &cmat[nbf*nconstr];
    pkv_TransposeMatrixf ( nconstr, nbf, nbf, privateG1->ECmat,
                           nconstr, cmat );
    pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                               nconstr, nconstr, cmat );
    if ( !pkn_QRDecomposeMatrixf ( nbf, nconstr, cmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < nconstr; i++ )
      for ( j = i;  j < nconstr; j++ )
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
  cmat = privateG1->ECmat;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetExtRightSidef ( domain, privateG1->Q2EBMat,
                                &privateG1->Q2EBMat[nfunc_c*nfunc_b],
                                spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, x );
  pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixf ( nconstr, nbf, nbf, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixf ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolvef ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemf ( spdimen*nbf );
  if ( !y ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixf ( nconstr, nbf, nbf, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, y );
  pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, y );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemf ( spdimen*nbf );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtFillHoleConstrf*/

boolean g1h_Q2ExtFillHoleAltConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) )
{
  void    *sp;
  G1HolePrivateRecf *privateG1;
  float  *lmat, *acmat, *arcmat, *aa;
  float  *fc00, *b, *x, *y;
  float  s, t;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  int    i, j;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->AECmat ||
       privateG1->extnaconstr != naconstr || privateG1->extacdim != spdimen ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  if ( privateG1->Q2AMat )
    if ( !g1h_Q2ExtComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2ELMat )
    if ( !g1h_Q2ExtDecomposeMatrixf ( domain ) )
      goto failure;
  lmat = privateG1->Q2ELMat;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  acmat = privateG1->AECmat;
  if ( !(arcmat = privateG1->Q2AERCMat) ) {
    arcmat = privateG1->Q2ARCMat =
               (float*)malloc ( (naconstr*(naconstr+1))/2*sizeof(float) );
    if ( !arcmat )
      goto failure;
    acmat = pkv_GetScratchMemf ( (nbf*spdimen+2)*naconstr );
    if ( !acmat )
      goto failure;
    aa = &acmat[nbf*spdimen*naconstr];
    pkv_TransposeMatrixf ( naconstr, spdimen*nbf, spdimen*nbf,
                           privateG1->AECmat, naconstr, acmat );
    for ( i = 0; i < spdimen; i++ )
      pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                                 naconstr, naconstr, &acmat[i*nbf*naconstr] );
    if ( !pkn_QRDecomposeMatrixf ( nbf*spdimen, naconstr, acmat, aa ) ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
    for ( i = 0; i < naconstr; i++ )
      for ( j = i; j < naconstr; j++ )
        arcmat[pkn_SymMatIndex(i,j)] = acmat[i*naconstr+j];
    for ( i = 0; i < naconstr; i++ ) {
      for ( j = 0, s = 0.0;  j < i;  j++ ) {
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
  acmat = privateG1->AECmat;

  b = pkv_GetScratchMemf ( spdimen*nbf );
  x = pkv_GetScratchMemf ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetExtRightSidef ( domain, privateG1->Q2EBMat,
                &privateG1->Q2EBMat[nfunc_c*nfunc_b], spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, x );
  pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixf ( nbf, spdimen, spdimen, x, nbf, y );
  pkn_MultMatrixf ( naconstr, spdimen*nbf, spdimen*nbf, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixf ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolvef ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixf ( naconstr, nbf*spdimen, nbf*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixf ( spdimen, nbf, nbf, y, spdimen, b );
  pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, b );
  pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, lmat,
                             spdimen, spdimen, b );
  pkn_SubtractMatrixf ( 1, spdimen*nbf, 0, x, 0, b, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(float) );

  if ( !_g1h_OutputExtPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ExtFillHoleAltConstrf*/

