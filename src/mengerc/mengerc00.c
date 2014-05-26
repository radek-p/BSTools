
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "mengerc.h"

#include "mengercprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_TabBasisFunctions ( int n, int nqkn, mengerc_data *md )
{
  void   *sp;
  double *kn, *bsf, *bsf1, *bsf2, *dbsf, *dbsf1, *ddbsf;
  int    i, j, fnz, nnz;

  sp = pkv_GetScratchMemTop ();
  memset ( md, 0, sizeof(mengerc_data) );
  md->nqkn = nqkn;
  md->n = n;
  md->qkn = malloc ( nqkn*sizeof(double) );
  md->qc = malloc ( nqkn*sizeof(double) );
  md->bsf = malloc ( nqkn*(n+1)*sizeof(double) );
  md->dbsf = malloc ( nqkn*(n+1)*sizeof(double) );
  md->ddbsf = malloc ( nqkn*(n+1)*sizeof(double) );
  md->bsf1 = malloc ( nqkn*n*sizeof(double) );
  md->dbsf1 = malloc ( nqkn*n*sizeof(double) );
  if ( !md->qkn || !md->qc || !md->bsf || !md->ddbsf || !md->bsf1 || !md->dbsf1 )
    return false;
    /* tworzenie wezlow i wspolczynnikow kwadratury */
  switch ( nqkn ) {
case 1:
    if ( !pkn_QuadRectanglesd ( 0.0, 1.0, 1, md->qkn, md->qc ) )
      return false;
    break;
case 2:
    if ( !pkn_QuadGaussLegendre4d ( 0.0, 1.0, 2, md->qkn, md->qc ) )
      return false;
    break;
case 4:
    if ( !pkn_QuadGaussLegendre8d ( 0.0, 1.0, 4, md->qkn, md->qc ) )
      return false;
    break;
case 6:
    if ( !pkn_QuadGaussLegendre12d ( 0.0, 1.0, 6, md->qkn, md->qc ) )
      return false;
    break;
case 8:
    if ( !pkn_QuadGaussLegendre16d ( 0.0, 1.0, 8, md->qkn, md->qc ) )
      return false;
    break;
case 10:
    if ( !pkn_QuadGaussLegendre20d ( 0.0, 1.0, 10, md->qkn, md->qc ) )
      return false;
    break;
default:
    return false;
  }

  kn = pkv_GetScratchMemd ( 3*n+2 );
  if ( !kn )
    return false;
  for ( i = 0; i < 2*n+3; i++ )
    kn[i] = (double)i;
  bsf2 = &kn[2*n+3];
  for ( i = 0; i < nqkn; i++ ) {
    bsf = &md->bsf[i*(n+1)];
    bsf1 = &md->bsf1[i*n];
    dbsf = &md->dbsf[i*(n+1)];
    dbsf1 = &md->dbsf1[i*n];
    ddbsf = &md->ddbsf[i*(n+1)];
        /* oblicz wartosci funkcji bazowych stopnia n */
    mbs_deBoorBasisd ( n, 2*n+2, kn, (double)n+md->qkn[i],
                       &fnz, &nnz, bsf );
        /* oblicz wartosci funkcji bazowych stopnia n-1 */
    mbs_deBoorBasisd ( n-1, 2*n, kn, (double)(n-1)+md->qkn[i],
                       &fnz, &nnz, bsf1 );
        /* oblicz wartosci funkcji bazowych stopnia n-2 */
    mbs_deBoorBasisd ( n-2, 2*n-2, kn, (double)(n-2)+md->qkn[i],
                       &fnz, &nnz, bsf2 );
        /* oblicz pochodne funkcji bazowych stopnia n */
    dbsf[0] = -bsf1[0];
    for ( j = 1; j < n; j++ )
      dbsf[j] = bsf1[j-1]-bsf1[j];
    dbsf[n] = bsf1[n-1];
        /* oblicz pochodne funkcji bazowych stopnia n-1 */
    dbsf1[0] = -bsf2[0];
    for ( j = 1; j < n-1; j++ )
      dbsf1[j] = bsf2[j-1]-bsf2[j];
    dbsf1[n-1] = bsf2[n-2];
        /* oblicz pochodne drugiego rzedu funkcji bazowych stopnia n */
    ddbsf[0] = -dbsf1[0];
    for ( j = 1; j < n; j++ )
      ddbsf[j] = dbsf1[j-1]-dbsf1[j];
    ddbsf[n] = dbsf1[n-1];
  }

  pkv_SetScratchMemTop ( sp );
  return true;
} /*mengerc_TabBasisFunctions*/

boolean mengerc_BindACurve ( mengerc_data *md,
                             int deg, int lkn, double *knots, point3d *cpoints,
                             int nqkn, double w, double *kara, boolean alt_scale )
{
  int      ncp, i;
  double   acp;
  point3d  sc;
  int      n, nn;

  memset ( md, 0, sizeof(mengerc_data) );
  if ( !mengerc_TabBasisFunctions ( deg, nqkn, md ) )
    goto failure;
  md->lkn = lkn;
  md->knots = knots;
  md->cpoints = cpoints;
  md->L = 1.0;
  md->w = w;
  md->alt_scale = alt_scale;
  if ( kara )
    memcpy ( md->kara, kara, MENGERC_NPPARAM*sizeof(double) );
  else
    for ( i = 0; i < MENGERC_NPPARAM; i++ )
      md->kara[i] = 1.0;
  if ( !mengerc_intD ( md, lkn, knots, cpoints, &md->L, &acp ) )
    goto failure;
  ncp = lkn-deg;
  mengerc_GravityCentre ( ncp-deg, cpoints, &sc );
  md->mdi = mengerc_FindRemotestPoint ( ncp-deg, cpoints, &sc );
  n = 3*(ncp-deg);
  nn = (n*(n+1))/2;
  md->fx = malloc ( (15*n+6*nn)*sizeof(double) );
  md->mcpoints = malloc ( 3*(lkn-2*deg)*sizeof(double) );
  if ( !md->fx || !md->mcpoints )
    goto failure;
  md->gx = &md->fx[n];
  md->hx = &md->gx[n];
  md->ggkM = &md->hx[n];  md->ggR1 = &md->ggkM[n];  md->ggR2 = &md->ggR1[n];
  md->ggR3 = &md->ggR2[n];  md->ggR4 = &md->ggR3[n];  md->ggR5 = &md->ggR4[n];
  md->hgkM = &md->ggR5[n];  md->hgR1 = &md->hgkM[n];  md->hgR2 = &md->hgR1[n];
  md->hgR3 = &md->hgR2[n];  md->hgR4 = &md->hgR3[n];  md->hgR5 = &md->hgR4[n];
  md->hhkM = &md->hgR5[n];  md->hhR1 = &md->hhkM[nn];  md->hhR2 = &md->hhR1[nn];
  md->hhR3 = &md->hhR2[nn];  md->hhR4 = &md->hhR3[nn];  md->hhR5 = &md->hhR4[nn];
  memset ( md->gx, 0, 3*n*sizeof(double) );
  return true;

failure:
  mengerc_UntieTheCurve ( md );
  return false;
} /*mengerc_BindACurve*/

void mengerc_UntieTheCurve ( mengerc_data *md )
{
  if ( md->qkn )      free ( md->qkn );
  if ( md->qc )       free ( md->qc );
  if ( md->bsf )      free ( md->bsf );
  if ( md->dbsf )     free ( md->dbsf );
  if ( md->ddbsf )    free ( md->ddbsf );
  if ( md->bsf1 )     free ( md->bsf1 );
  if ( md->dbsf1 )    free ( md->dbsf1 );
  if ( md->fx )       free ( md->fx );
  if ( md->mcpoints ) free ( md->mcpoints );
} /*mengerc_UntieTheCurve*/

