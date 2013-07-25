
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"

#define EPS 1.0e-7

/* ////////////////////////////////////////////////////////////////////////// */
void _g2mbl_Testfghd ( void *data, char *fn )
{
  void             *sp;
  mesh_lmt_optdata *d;
  int              nv;
  point3d          *mvcp, *auxmvcp;
  int              ndomelems, *domelcpind;
  int              nvcp, nvars, *nncpi, *vncpi;
  meshdom_elem     *domelem;
  double           *ftab, *gtab, *htab;
  int              hsize, *hprof;
  double           **hrows, **lhrows;
  double           *coeff, *dcoeff, *grad;
  double           func, f1, f2, sca;
  int              i, j;
  FILE             *f;

  sp = pkv_GetScratchMemTop ();
        /*extract data to local variables */
  d = data;
  nv = d->nv;
  mvcp = d->mvcp;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  ftab = d->ftab;
  gtab = d->gtab;
  htab = d->htab;
  nvcp = d->nvcp;
  nvars = d->nvars;
  nncpi = d->nncpi;
  vncpi = d->vncpi;
  hsize = d->hsize;
  hprof = d->hprof;
  hrows = d->hrows;
  lhrows = d->lhrows;
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvars );
  if ( !coeff )
    goto failure;
  dcoeff = &coeff[nvars];
  grad = &dcoeff[nvars];
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp )
    goto failure;
  memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );

  if ( !g2mbl_UFuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                          nv, mvcp, ndomelems, domelem, domelcpind, nncpi,
                          true, ftab, gtab, htab,
                          &func, nvars, grad, hsize, hprof, hrows ) )
      goto failure;

  f = fopen ( fn, "w+" );
  if ( !f )
    goto failure;

  fprintf ( f, "func = %17.10e\ngrad:\n", func );
  for ( i = 0; i < nvcp; i++ ) {
    sca = auxmvcp[vncpi[i]].x;
    auxmvcp[vncpi[i]].x += EPS;
    f1 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    auxmvcp[vncpi[i]].x = sca-EPS;
    f2 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    coeff[3*i] = (f1-f2)/(2.0*EPS);
    auxmvcp[vncpi[i]].x = sca;

    sca = auxmvcp[vncpi[i]].y;
    auxmvcp[vncpi[i]].y += EPS;
    f1 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    auxmvcp[vncpi[i]].y = sca-EPS;
    f2 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    coeff[3*i+1] = (f1-f2)/(2.0*EPS);
    auxmvcp[vncpi[i]].y = sca;

    sca = auxmvcp[vncpi[i]].z;
    auxmvcp[vncpi[i]].z += EPS;
    f1 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    auxmvcp[vncpi[i]].z = sca-EPS;
    f2 = g2mbl_UFuncd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                        ndomelems, domelem, domelcpind, true, ftab );
    coeff[3*i+2] = (f1-f2)/(2.0*EPS);
    auxmvcp[vncpi[i]].z = sca;
  }
  for ( i = 0; i < nvars; i++ )
    fprintf ( f, "%3d: %17.10e %17.10e\n", i, grad[i], coeff[i] );

  fprintf ( f, "hessian:\n" );
  for ( i = 0; i < nvcp; i++ ) {
    sca = auxmvcp[vncpi[i]].x;
    auxmvcp[vncpi[i]].x += EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, coeff );
    auxmvcp[vncpi[i]].x = sca-EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, grad );
    auxmvcp[vncpi[i]].x = sca;
    for ( j = hprof[3*i]; j <= 3*i; j++ )
      lhrows[3*i][j] = (coeff[j]-grad[j])/(2.0*EPS);

    sca = auxmvcp[vncpi[i]].y;
    auxmvcp[vncpi[i]].y += EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, coeff );
    auxmvcp[vncpi[i]].y = sca-EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, grad );
    auxmvcp[vncpi[i]].y = sca;
    for ( j = hprof[3*i+1]; j <= 3*i+1; j++ )
      lhrows[3*i+1][j] = (coeff[j]-grad[j])/(2.0*EPS);

    sca = auxmvcp[vncpi[i]].z;
    auxmvcp[vncpi[i]].z += EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, coeff );
    auxmvcp[vncpi[i]].z = sca-EPS;
    g2mbl_UFuncGradd ( d->nkn1, d->aqcoeff, d->aNitabs, d->aJac, nv, auxmvcp,
                       ndomelems, domelem, domelcpind, nncpi, true, ftab, gtab,
                       &f1, nvars, grad );
    auxmvcp[vncpi[i]].z = sca;
    for ( j = hprof[3*i+2]; j <= 3*i+2; j++ )
      lhrows[3*i+2][j] = (coeff[j]-grad[j])/(2.0*EPS);
  }
  for ( i = 0; i < nvars; i++ ) {
    for ( j = hprof[i]; j <= i; j++ )
      fprintf ( f, "%3d %3d: %17.10e %17.10e\n", i, j, hrows[i][j], lhrows[i][j] );
    fprintf ( f, "\n" );
  }

  fclose ( f );

failure:
  pkv_SetScratchMemTop ( sp );
  return;
} /*_g2mbl_Testfghd*/

void _g2mbl_GetHessianProfiled ( void *data, int *nvars, int **hprof )
{
  mesh_lmt_optdata *d;

  d = data;
  *nvars = d->nvars;
  *hprof = d->hprof;
} /*_g2mbl_GetHessianProfiled*/

void _g2mbl_GetVariablePointOrder ( void *data, int *nvcp, int **vncpi )
{
  mesh_lmt_optdata *d;

  d = data;
  *nvcp = d->nvcp;
  *vncpi = d->vncpi;
} /*_g2mbl_GetVariablePointOrder*/

void _g2mbl_GetHessiand ( void *data, int *nvars, int **hprof, int *hsize,
                          double ***hrows, double ***lhrows )
{
  void             *sp;
  mesh_lmt_optdata *d;
  int              nv;
  point3d          *mvcp;
  int              ndomelems, *domelcpind;
  int              *nncpi;
  meshdom_elem     *domelem;
  double           *ftab, *gtab, *htab;
  double           func, *grad;

  sp = pkv_GetScratchMemTop ();
        /*extract data to local variables */
  d = data;
  nv = d->nv;
  mvcp = d->mvcp;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  ftab = d->ftab;
  gtab = d->gtab;
  htab = d->htab;
  *nvars = d->nvars;
  nncpi = d->nncpi;
  *hsize = d->hsize;
  *hprof = d->hprof;
  *hrows = d->hrows;
  *lhrows = d->lhrows;
  grad = pkv_GetScratchMemd ( *nvars );
  if ( !grad )
    goto failure;
  if ( !g2mbl_UFuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                          nv, mvcp, ndomelems, domelem, domelcpind, nncpi,
                          true, ftab, gtab, htab,
                          &func, *nvars, grad, *hsize, *hprof, *hrows ) ) {
failure:
    *nvars = *hsize = 0;
    *hprof = NULL;
    *hrows = *lhrows = NULL;
  }
  pkv_SetScratchMemTop ( sp );
} /*_g2mbl_GetHessiand*/

