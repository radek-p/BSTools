
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g1h_V0SpaceDimd ( GHoleDomaind *domain )
{
  G1HolePrivateRecd *privateG1;

  privateG1 = domain->privateG1;
  if ( !privateG1 )
    return 0;
  else
    return privateG1->nfunc_a;
} /*g1h_V0SpaceDimd*/

/* ////////////////////////////////////////////////////////////////////////// */
/* Get derivatives of order 0 ... 2 of the basis function patches             */
/* in the direction of the boundar curves cno, for all basis functions.       */
/* This procedure may be useful for setting up the constraints.               */
/* In more general cases mixed partial derivatives might be needed, but this  */
/* possibility is not serviced at the moment.                                 */
boolean g1h_GetBPDerivativesd ( GHoleDomaind *domain, int cno, double *val )
{
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a;
  double *bbr0, *cr0, *v;
  int    i;

  hole_k = domain->hole_k;
  if ( cno < 0 || cno >= hole_k )
    goto failure;

  privateG1 = domain->privateG1;
  if ( !privateG1 )
    goto failure;
  bbr0 = privateG1->basis_a;
  if ( !bbr0 )
    goto failure;

  nfunc_a = privateG1->nfunc_a;
  for ( i = 0; i < nfunc_a; i++ ) {
    cr0 = &bbr0[(i*hole_k+cno)*(G1_CROSS00DEG+1)];
    v = &val[3*i];
    if ( !mbs_BCHornerDer2C1d ( G1_CROSS00DEG, cr0, 0.0, &v[0], &v[1], &v[2] ) )
      goto failure;
  }
  return true;

failure:
  return false;
} /*g1h_GetBPDerivativesd*/

boolean g1h_GetBFuncPatchd ( GHoleDomaind *domain, int fn, int pn, double *bp )
{
  G1HolePrivateRecd *privateG1;
  int    nfunc_a, nfunc_b, degu, degv;
  double *bc00, *bc01, *bc10, *bc11,
         *bd00, *bd01, *bd10, *bd11;
  double zero[] = { 0.0, 0.0 };

  privateG1 = domain->privateG1;
  if ( !privateG1 )
    return false;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  
  if ( fn < 0 )
    return false;
  else if ( fn < nfunc_a ) {
    _g1h_GetBFAPatchCurvesd ( domain, fn, pn, &bc00, &bc01, &bd00, &bd01 );
    mbs_BezC1CoonsToBezd ( 1,
               G1_CROSS00DEG, bc00, G1_CROSS01DEG, bc01, 1, zero, 1, zero,
               G1_CROSS00DEG, bd00, G1_CROSS01DEG, bd01, 1, zero, 1, zero,
               &degu, &degv, bp );
  }
  else if ( fn < nfunc_a+nfunc_b ) {
    _g1h_GetBFBPatchCurvesd ( domain, fn-nfunc_a, pn,
                              &bc00, &bc01, &bc10, &bc11,
                              &bd00, &bd01, &bd10, &bd11 );
    mbs_BezC1CoonsToBezd ( 1,
               G1_CROSS00DEG, bc00, G1_CROSS01DEG, bc01,
               G1_CROSS10DEG, bc10, G1_CROSS11DEG, bc11,
               G1_CROSS00DEG, bd00, G1_CROSS01DEG, bd01,
               G1_CROSS10DEG, bd10, G1_CROSS11DEG, bd11, &degu, &degv, bp );
  }
  else
    return false;

  return true;
} /*g1h_GetBFuncPatchd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g1h_SetConstraintMatrixd ( GHoleDomaind *domain,
                                   int nconstr, const double *cmat )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *Cmat, *RCmat, *aa, s, t;
  int    nfunc_a, i, j;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( privateG1->Cmat ) free ( privateG1->Cmat );  /* forget previous constraints */
  if ( privateG1->RCmat ) free ( privateG1->RCmat );
  if ( privateG1->Q2RCMat ) {
    free ( privateG1->Q2RCMat );
    privateG1->Q2RCMat = NULL;
  }

  if ( !privateG1->Lmat ) {
    if ( !g1h_DecomposeMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG1->nfunc_a;
  if ( nconstr <= 0 || nconstr > nfunc_a )
    goto failure;

  Cmat = privateG1->Cmat = malloc ( nfunc_a*nconstr*sizeof(double) );
  RCmat = privateG1->RCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
  if ( !Cmat || !RCmat )
    goto failure;

  memcpy ( Cmat, cmat, nfunc_a*nconstr*sizeof(double) );
  Cmat = pkv_GetScratchMemd ( (nfunc_a+2)*nconstr );
  if ( !Cmat )
    goto failure;

  aa = &Cmat[nfunc_a*nconstr];
  pkv_TransposeMatrixd ( nconstr, nfunc_a, nfunc_a, privateG1->Cmat, nconstr, Cmat );
  pkn_LowerTrMatrixSolved ( nfunc_a, privateG1->Lmat,
                            nconstr, nconstr, Cmat, nconstr, Cmat );
  if ( !pkn_QRDecomposeMatrixd ( nfunc_a, nconstr, Cmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      RCmat[pkn_SymMatIndex(i,j)] = Cmat[i*nconstr+j];
  for ( i = 0; i < nconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = RCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = RCmat[pkn_SymMatIndex(i,i)];
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG1->nconstr = nconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG1->Cmat ) free ( privateG1->Cmat );
  if ( privateG1->RCmat ) free ( privateG1->RCmat );
  privateG1->Cmat = privateG1->RCmat = NULL;
  privateG1->nconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetConstraintMatrixd*/

boolean g1h_SetAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int naconstr, const double *acmat )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *ACmat, *ARCmat, *aa, s, t;
  int    nfunc_a, i, j;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( privateG1->ACmat ) free ( privateG1->ACmat );  /* forget previous constraints */
  if ( privateG1->ARCmat ) free ( privateG1->ARCmat );
  if ( privateG1->Q2ARCMat ) {
    free ( privateG1->Q2ARCMat );
    privateG1->Q2ARCMat = NULL;
  }

  if ( !privateG1->Lmat ) {
    if ( !g1h_DecomposeMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG1->nfunc_a;
  if ( naconstr <= 0 || naconstr > nfunc_a )
    goto failure;

  ACmat = privateG1->ACmat = malloc ( nfunc_a*spdimen*naconstr*sizeof(double) );
  ARCmat = privateG1->ARCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
  if ( !ACmat || !ARCmat )
    goto failure;

  memcpy ( ACmat, acmat, nfunc_a*spdimen*naconstr*sizeof(double) );
  ACmat = pkv_GetScratchMemd ( (nfunc_a*spdimen+2)*naconstr );
  if ( !ACmat )
    goto failure;

  aa = &ACmat[nfunc_a*spdimen*naconstr];
  pkv_TransposeMatrixd ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a,
                         privateG1->ACmat, naconstr, ACmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_LowerTrMatrixSolved ( nfunc_a, privateG1->Lmat,
                              naconstr, naconstr, &ACmat[i*nfunc_a*naconstr],
                              naconstr, &ACmat[i*nfunc_a*naconstr] );
  if ( !pkn_QRDecomposeMatrixd ( nfunc_a*spdimen, naconstr, ACmat, aa ) ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  for ( i = 0; i < naconstr; i++ )
    for ( j = i; j < naconstr; j++ )
      ARCmat[pkn_SymMatIndex(i,j)] = ACmat[i*naconstr+j];
  for ( i = 0; i < naconstr; i++ ) {
    for ( j = 0, s = 0.0;  j < i;  j++ ) {
      t = ARCmat[pkn_SymMatIndex(i,j)];
      s += t*t;
    }
    t = ARCmat[pkn_SymMatIndex(i,i)];
    if ( t*t < 1.0e-10*s ) {
      domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG1->acdim = spdimen;
  privateG1->naconstr = naconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG1->ACmat ) free ( privateG1->ACmat );
  if ( privateG1->ARCmat ) free ( privateG1->ARCmat );
  privateG1->ACmat = privateG1->ARCmat = NULL;
  privateG1->naconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SetAltConstraintMatrixd*/

