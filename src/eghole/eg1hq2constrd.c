
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
boolean g1h_Q2FillHoleConstrd ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                int nconstr, CONST_ double *constr,
                                double *acoeff, void *usrptr,
                                void (*outpatch) ( int n, int m, const double *cp,
                                                   void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double  *lmat, *cmat, *rcmat, *aa, *b, *x, *y;
  double  *fc00;
  double  s, t;
  int    i, j;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Cmat || privateG1->nconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixd ( domain ) )
      goto failure;
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  if ( !(rcmat = privateG1->Q2RCMat) ) {
    rcmat = privateG1->Q2RCMat =
              (double*)malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
    if ( !rcmat )
      goto failure;
    cmat = pkv_GetScratchMemd ( (nfunc_a+2)*nconstr );
    if ( !cmat )
      goto failure;
    aa = &cmat[nfunc_a*nconstr];
    pkv_TransposeMatrixd ( nconstr, nfunc_a, nfunc_a, privateG1->Cmat,
                           nconstr, cmat );
    pkn_LowerTrMatrixSolved ( nfunc_a, lmat, nconstr, nconstr, cmat,
                             nconstr, cmat  );
    if ( !pkn_QRDecomposeMatrixd ( nfunc_a, nconstr, cmat, aa ) ) {
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

  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSided ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixd ( nconstr, nfunc_a, nfunc_a, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixd ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemd ( spdimen*nfunc_a );
  if ( !y ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixd ( nconstr, nfunc_a, nfunc_a, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemd ( spdimen*nfunc_a );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHoleConstrd*/

boolean g1h_Q2FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double  *lmat, *acmat, *arcmat, *aa, *b, *x, *y;
  double  *fc00;
  double  s, t;
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
    if ( !g1h_Q2ComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixd ( domain ) )
      goto failure;
    
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  if ( !(arcmat = privateG1->Q2ARCMat) ) {
    arcmat = privateG1->Q2ARCMat =
              (double*)malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
    if ( !arcmat )
      goto failure;
    acmat = pkv_GetScratchMemd ( (nfunc_a*spdimen+2)*naconstr );
    if ( !acmat )
      goto failure;
    aa = &acmat[nfunc_a*spdimen*naconstr];
    pkv_TransposeMatrixd ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a,
                           privateG1->ACmat, naconstr, acmat );
    for ( i = 0; i < spdimen; i++ )
      pkn_LowerTrMatrixSolved ( nfunc_a, lmat, naconstr,
                                naconstr, &acmat[i*nfunc_a*naconstr],
                                naconstr, &acmat[i*nfunc_a*naconstr] );
    if ( !pkn_QRDecomposeMatrixd ( nfunc_a*spdimen, naconstr, acmat, aa ) ) {
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

  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetRightSided ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
              /* now x is the solution without the constraints and it is time */
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixd ( nfunc_a, spdimen, spdimen, x, nfunc_a, y );
  pkn_MultMatrixd ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixd ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixd ( naconstr, nfunc_a*spdimen, nfunc_a*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixd ( spdimen, nfunc_a, nfunc_a, y, spdimen, b );
  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, y );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, y, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nfunc_a, 0, x, 0, y, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHoleAltConstrd*/

