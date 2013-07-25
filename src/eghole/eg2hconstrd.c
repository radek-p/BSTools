
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
int g2h_V0SpaceDimd ( GHoleDomaind *domain )
{
  G2HolePrivateRecd *privateG2;

  privateG2 = domain->privateG2;
  if ( !privateG2 )
    return 0;
  else
    return privateG2->nfunc_a;
} /*g2h_V0SpaceDimd*/

/* ////////////////////////////////////////////////////////////////////////// */
/* Get derivatives of order 0 ... 4 of the basis function patches             */
/* in the direction of the boundar curves cno, for all basis functions.       */
/* This procedure may be useful for setting up the constraints.               */
/* In more general cases mixed partial derivatives might be needed, but this  */
/* possibility is not serviced at the moment.                                 */
boolean g2h_GetBPDerivativesd ( GHoleDomaind *domain, int cno, double *val )
{
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a;
  double *bbr0, *cr0, *v;
  int    i, j, k, l;

  hole_k = domain->hole_k;
  if ( cno < 0 || cno >= hole_k )
    goto failure;

  privateG2 = domain->privateG2;
  if ( !privateG2 )
    goto failure;
  bbr0 = privateG2->basis_a;
  if ( !bbr0 )
    goto failure;

  nfunc_a = privateG2->nfunc_a;
  for ( i = 0; i < nfunc_a; i++ ) {
    cr0 = &bbr0[(i*hole_k+cno)*(G2_CROSS00DEG+1)];
    v = &val[5*i];
        /* currently there is no procedure computing derivatives */
        /* of Bezier curves of order 4, therefore they are computed "by hand" */
    memcpy ( v, cr0, 5*sizeof(double) );
    for ( j = 1; j <= 4; j++ )
      for ( k = 4; k >= j; k-- )
        v[k] -= v[k-1];
    for ( j = 1, k = l = G2_CROSS00DEG; j <= 4; j++, l *= --k )
      v[j] *= (double)l;
  }
  return true;

failure:
  return false;
} /*g2h_GetBPDerivativesd*/

boolean g2h_GetBFuncPatchd ( GHoleDomaind *domain, int fn, int pn, double *bp )
{
  G2HolePrivateRecd *privateG2;
  int    nfunc_a, nfunc_b, degu, degv;
  double *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  double zero[] = { 0.0, 0.0 };

  privateG2 = domain->privateG2;
  if ( !privateG2 )
    return false;

  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  
  if ( fn < 0 )
    return false;
  else if ( fn < nfunc_a ) {
    _g2h_GetBFAPatchCurvesd ( domain, fn, pn,
                              &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
    mbs_BezC2CoonsToBezd ( 1,
               G2_CROSS00DEG, bc00, G2_CROSS01DEG, bc01, G2_CROSS02DEG, bc02,
               1, zero, 1, zero, 1, zero,
               G2_CROSS00DEG, bd00, G2_CROSS01DEG, bd01, G2_CROSS02DEG, bd02,
               1, zero, 1, zero, 1, zero, &degu, &degv, bp );
  }
  else if ( fn < nfunc_a+nfunc_b ) {
    _g2h_GetBFBPatchCurvesd ( domain, fn-nfunc_a, pn,
                              &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                              &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
    mbs_BezC2CoonsToBezd ( 1,
               G2_CROSS00DEG, bc00, G2_CROSS01DEG, bc01, G2_CROSS02DEG, bc02,
               G2_CROSS10DEG, bc10, G2_CROSS11DEG, bc11, G2_CROSS12DEG, bc12,
               G2_CROSS00DEG, bd00, G2_CROSS01DEG, bd01, G2_CROSS02DEG, bd02,
               G2_CROSS10DEG, bd10, G2_CROSS11DEG, bd11, G2_CROSS02DEG, bd12,
               &degu, &degv, bp );
  }
  else
    return false;

  return true;
} /*g2h_GetBFuncPatchd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2h_SetConstraintMatrixd ( GHoleDomaind *domain,
                                   int nconstr, const double *cmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double *Cmat, *RCmat, *aa, s, t;
  int    nfunc_a, i, j;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( privateG2->Cmat ) free ( privateG2->Cmat );  /* forget previous constraints */
  if ( privateG2->RCmat ) free ( privateG2->RCmat );

  if ( !privateG2->Lmat ) {
    if ( !g2h_DecomposeMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG2->nfunc_a;
  if ( nconstr <= 0 || nconstr > nfunc_a )
    goto failure;

  Cmat = privateG2->Cmat = malloc ( nfunc_a*nconstr*sizeof(double) );
  RCmat = privateG2->RCmat = malloc ( (nconstr*(nconstr+1))/2*sizeof(double) );
  if ( !Cmat || !RCmat )
    goto failure;

  memcpy ( Cmat, cmat, nfunc_a*nconstr*sizeof(double) );
  Cmat = pkv_GetScratchMemd ( (nfunc_a+2)*nconstr );
  if ( !Cmat )
    goto failure;

  aa = &Cmat[nfunc_a*nconstr];
  pkv_TransposeMatrixd ( nconstr, nfunc_a, nfunc_a, privateG2->Cmat, nconstr, Cmat );
  pkn_LowerTrMatrixSolved ( nfunc_a, privateG2->Lmat,
                            nconstr, nconstr, Cmat, nconstr, Cmat );
  if ( !pkn_QRDecomposeMatrixd ( nfunc_a, nconstr, Cmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG2->nconstr = nconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG2->Cmat ) free ( privateG2->Cmat );
  if ( privateG2->RCmat ) free ( privateG2->RCmat );
  privateG2->Cmat = privateG2->RCmat = NULL;
  privateG2->nconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SetConstraintMatrixd*/

boolean g2h_SetAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int naconstr, const double *acmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  double *ACmat, *ARCmat, *aa, s, t;
  int    nfunc_a, i, j;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  if ( privateG2->ACmat ) free ( privateG2->ACmat );  /* forget previous constraints */
  if ( privateG2->ARCmat ) free ( privateG2->ARCmat );

  if ( !privateG2->Lmat ) {
    if ( !g2h_DecomposeMatrixd ( domain ) )
      goto failure;
  }

  nfunc_a = privateG2->nfunc_a;
  if ( naconstr <= 0 || naconstr > nfunc_a )
    goto failure;

  ACmat = privateG2->ACmat = malloc ( nfunc_a*spdimen*naconstr*sizeof(double) );
  ARCmat = privateG2->ARCmat = malloc ( (naconstr*(naconstr+1))/2*sizeof(double) );
  if ( !ACmat || !ARCmat )
    goto failure;

  memcpy ( ACmat, acmat, nfunc_a*spdimen*naconstr*sizeof(double) );
  ACmat = pkv_GetScratchMemd ( (nfunc_a*spdimen+2)*naconstr );
  if ( !ACmat )
    goto failure;

  aa = &ACmat[nfunc_a*spdimen*naconstr];
  pkv_TransposeMatrixd ( naconstr, spdimen*nfunc_a, spdimen*nfunc_a,
                         privateG2->ACmat, naconstr, ACmat );
  for ( i = 0; i < spdimen; i++ )
    pkn_LowerTrMatrixSolved ( nfunc_a, privateG2->Lmat,
                              naconstr, naconstr, &ACmat[i*nfunc_a*naconstr],
                              naconstr, &ACmat[i*nfunc_a*naconstr] );
  if ( !pkn_QRDecomposeMatrixd ( nfunc_a*spdimen, naconstr, ACmat, aa ) ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
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
      domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
      goto failure;
    }
  }

  privateG2->acdim = spdimen;
  privateG2->naconstr = naconstr;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( privateG2->ACmat ) free ( privateG2->ACmat );
  if ( privateG2->ARCmat ) free ( privateG2->ARCmat );
  privateG2->ACmat = privateG2->ARCmat = NULL;
  privateG2->naconstr = 0;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SetAltConstraintMatrixd*/

