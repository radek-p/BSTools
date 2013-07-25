
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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

/* ////////////////////////////////////////////////////////////////////////// */
#define EPS0      1.0e-3
#define EPS1      1.0e-6
#define EPS2      1.0e-7
#define EPS3      1.0e-4
#define EPS4      1.0e-16
#define EPS6      1.0e-5
#define EPS7      1.0e-8
#define ALPHA     0.05
#define MAXNTN   20
#define MAXCG  1000
#define MAXGTN   10
#define MAXBTN   10
#define MAXCTN   20
#define THR1     0.9
#define THR2     0.9
#define THR3     0.5/*0.25*/
#define THR4     0.5
#define GRTHR    0.02

/* ///////////////////////////////////////////////////////////////////////// */
double _g2mbl_AuxAltBNuFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int nvcp, int *vncpi,
                           int ndomelems, int *domelind, meshdom_elem *domelem,
                           int *domelcpind, double *ftab,
                           int nvars, int hsize, int *hprof,
                           double *Hessian, double *LHessian, double **lhrows,
                           double *grad, double *dcoeff, point3d *auxmvcp,
                           double nu )
{
  if ( _g2bl_ShiftDecompHessiand ( nvars, hsize, hprof, Hessian,
                                   LHessian, lhrows, nu ) ) {
#ifdef _DEBUG
printf ( "F" );
#endif
    pkn_NRBLowerTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, dcoeff, 1, dcoeff );
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
    return g2mbl_AFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                          ndomelems, domelind, domelem, domelcpind, true, ftab );
  }
  else
    return MYINFINITY;
} /*_g2mbl_AuxAltBNuFuncd*/

boolean _g2mbl_CMPMultRTHR3x3d ( int nrowsa, int nnza, index3 *ai, double *ac,
                                 int ncolsb, int nnzb, index3 *bi, double *bc,
                                 int nnz1,
                                 int hsize, int *hprof, double **hrows );

static boolean _g2mbl_CMPSetupCoarseTerm ( mesh_lmt_optdata *d )
{
  void          *sp;
  double        *Hbl, *rmnzc;
  int           nHbl, *cHbl;
  nzHbl_rowdesc *iHbl;
  int           brmnnz, bHnnz;
  index3        *brmnzi, *Hblnzi;
  int           *hprof, hsize, nvcp, nwcp;
  double        **hrows;
  int           i, j, k, l, j0, j1;

  sp = pkv_GetScratchMemTop ();
  Hbl = d->Hbl;
  rmnzc = d->rmnzc;
  nHbl = d->Hblsize;
  iHbl = d->iHbl;
  cHbl = d->cHbl;
  nvcp = d->nvcp;
  nwcp = d->nwcp;
  brmnnz = d->pmnnz;
  brmnzi = d->pmnzi;
  hprof = d->phprof;
  hsize = d->phsize;
  hrows = d->phrows;
        /* find the sparse representation of the Hessian matrix */
  bHnnz = 2*nHbl-nvcp;
  Hblnzi = pkv_GetScratchMem ( bHnnz*sizeof(index3) );
  if ( !Hblnzi )
    goto failure;
   for ( i = l = 0;  i < nvcp;  i++ ) {
     j0 = iHbl[i].firsthbl;
     j1 = j0 + iHbl[i].nhbl - 1;
     for ( j = j0; j < j1; j++ ) {
       k = cHbl[j];
       Hblnzi[l  ].i = i;
       Hblnzi[l  ].j = k;
       Hblnzi[l++].k = 9*j;
       Hblnzi[l  ].i = k;
       Hblnzi[l  ].j = i;
       Hblnzi[l++].k = 9*j;
     }
     Hblnzi[l].i = Hblnzi[l].j = i;
     Hblnzi[l++].k = 9*j1;
   }
        /* multiply the block Hessian by the refinement matrix */
  if ( !_g2mbl_CMPMultRTHR3x3d ( nvcp, bHnnz, Hblnzi, Hbl,
                                 nwcp, brmnnz, brmnzi, rmnzc, d->nnz1,
                                 hsize, hprof, hrows ) )
    goto failure;


  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPSetupCoarseTerm*/

boolean _g2mbl_MultQixAltd ( int n, void *usrdata, const double *x, double *Qix )
{
  void             *sp;
  mesh_lmt_optdata *md;
  block_desc       *block;
  int              *pprof, *nncpi, *vncpi;
  double           **plrows, *aux;
  int              i, j, k, l, nb, nbp, nbl, maxbl;
  int              nvars1;
  int              *phprof, nvcp, nwcp, rmnnz;
  double           **phrows, *rmnzc, *px, *qx;
  index3           *rmnzi;

  sp = pkv_GetScratchMemTop ();
  md = (mesh_lmt_optdata*)usrdata;
  nbl = md->nbl;
  block = md->block;
        /* find the maximal block size */
  maxbl = block[0].nvars;
  for ( i = 1; i < nbl; i++ )
    maxbl = max ( maxbl, block[i].nvars );
  aux = pkv_GetScratchMemd ( maxbl );
  if ( !aux )
    goto failure;

  nncpi = md->nncpi;
  memset ( Qix, 0, n*sizeof(double) );
  for ( i = 0; i < nbl; i ++ ) {
    vncpi = block[i].vncpi;
    nbp = block[i].ncp;
    nb = block[i].nvars;
    pprof = block[i].hprof;
    plrows = block[i].lhrows;
        /* permute the right-hand side */
    for ( j = 0; j < nbp; j++ )
      memcpy ( &aux[3*j], &x[3*nncpi[vncpi[j]]], 3*sizeof(double) );
        /* solve the triangular systems */
    pkn_NRBLowerTrSolved ( nb, pprof, plrows[0], plrows,
                           1, 1, aux, 1, aux );
    pkn_NRBUpperTrSolved ( nb, pprof, plrows[0], plrows,
                           1, 1, aux, 1, aux );
        /* permute the variables */
    for ( j = 0; j < nbp; j++ ) {
      k = 3*nncpi[vncpi[j]];
      l = 3*j;
      Qix[k]   += aux[l];
      Qix[k+1] += aux[l+1];
      Qix[k+2] += aux[l+2];
    }
  }
  pkv_SetScratchMemTop ( sp );
        /* use the coarse mesh preconditioner term, if present */
  if ( md->phrows ) {

    phprof = md->phprof;
    phrows = md->phrows;
    nvcp   = md->nvcp;     /* n == 3*nvcp */
    nwcp   = md->nwcp;
    nvars1 = 3*nwcp;
    rmnnz = md->pmnnz;
    rmnzi = md->pmnzi;
    rmnzc = md->rmnzc;
    px = pkv_GetScratchMemd ( n+nvars1 );
    if ( !px )
      goto way_out;
    qx = &px[nvars1];
    if ( !pkn_MultSPsubMTVectord ( n, nvars1, rmnnz, rmnzi, rmnzc,
                                   3, x, px ) )
      goto way_out;
    pkn_NRBLowerTrSolved ( nvars1, phprof, phrows[0], phrows, 1, 1, px, 1, px );
    pkn_NRBUpperTrSolved ( nvars1, phprof, phrows[0], phrows, 1, 1, px, 1, px );
    if ( !pkn_MultSPsubMVectord ( nvcp, nwcp, rmnnz, rmnzi, rmnzc, 3, px, qx ) )
      goto way_out;
    pkn_AddMatrixd ( 1, n, 0, Qix, 0, qx, 0, Qix );
  }
way_out:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MultQixAltd*/

boolean g2mbl_IterBlSurfaceOptAltBLMTd ( void *data, boolean *finished )
{
  void             *sp;
  mesh_lmt_optdata *d;
  block_desc       *bd;
  int              nv;
  point3d          *mvcp, *auxmvcp;
  int              ndomelems, bdomelems, *domelcpind, *domelind;
  nzHbl_rowdesc    *iHbl;
  int              *cHbl, Hblsize;
  double           *Hbl;
  int              nbl, nvcp, nbcp, nvars, *nncpi, *vncpi, nkn;
  meshdom_elem     *domelem;
  double           *ftab, *gtab, *htab, **Nitabs, **Jac, *qcoeff;
  int              hsize, *hprof;
  double           *Hessian, **hrows, *LHessian, **lhrows;
  double           *coeff, *dcoeff, *grad;
  double           func, gn, lco, ldco;
  boolean          positive;
  int              i, k, itm;
  double           lmin, lmax, ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge;

#ifdef _DEBUG
#define SUCCESS_CG \
  { \
printf ( "func = %12.7e\n", d->func ); \
    *finished = true; \
    goto next_iter; \
  }
#else
#define SUCCESS_CG \
  { \
    *finished = true; \
    goto next_iter; \
  }
#endif

#define MLFUNC(nu) \
_g2mbl_AuxAltBNuFuncd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, nbcp, vncpi, \
        bdomelems, domelind, domelem, domelcpind, ftab, nvars, hsize, hprof, \
        Hessian, LHessian, lhrows, grad, dcoeff, auxmvcp, nu )

#define RECORD_MIN(g,fnu) \
  { if ( fnu < fge ) { \
      memcpy ( coeff, dcoeff, nvars*sizeof(double) ); \
      bd->func = fge = fnu; \
      bd->nu = ge = g; \
      d->newpoint = true; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
  *finished = false;
        /* extract data to local variables */
  d = data;
#ifdef _DEBUG
printf ( "%d:", d->itn );
#endif
  fga = fgb = fgc = d->func;
  nv = d->nv;
  mvcp = d->mvcp;
  nvcp = d->nvcp;
  ftab = d->ftab;
  gtab = d->gtab;
  htab = d->htab;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  nbl = d->nbl;

        /* try a Newton iteration for the entire system of equations */
  if ( d->itn < d->next_entire || d->ibl > 0 )
    goto try_with_one_block;
  for ( i = 0; i < nbl; i++ )
    if ( d->block[i].hblstat < 3 )
      goto try_with_one_block;

        /* assign the other data to local variables */
  nvars = d->nvars;
  nncpi = d->nncpi;
  vncpi = d->vncpi;
  nkn    = d->nkn2;
  Nitabs = d->bNitabs;
  Jac    = d->bJac;
  qcoeff = d->bqcoeff;
  Hblsize = d->Hblsize;
  iHbl = d->iHbl;
  cHbl = d->cHbl;
  Hbl = d->Hbl;

  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  grad = pkv_GetScratchMemd ( 3*nvars );
  if ( !auxmvcp || !grad )
    goto failure;
  coeff = &grad[nvars];
  dcoeff = &grad[nvars];

  memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
  lco = 0.0;
  for ( i = 0; i < nvcp; i++ )
    lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
  lco = sqrt ( lco );

  d->last_ntn = true;
  if ( d->use_cg ) {
#ifdef _DEBUG
printf ( "H*" );
#endif
#ifdef COUNT
  d->nH ++; 
#endif
    if ( !g2mbl_B3x3FuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                          nv, mvcp, ndomelems, domelem, domelcpind, nncpi,
                          true, ftab, gtab, htab,
                          &func, nvars, grad, Hblsize, iHbl, cHbl, Hbl ) )
      goto failure;
        /* prepare the preconditioner */
#ifdef _DEBUG
printf ( "(cg) " );
#endif
    for ( i = 0; i < nbl; i++ ) {
      bd = &d->block[i];
      if ( !g2mbl_AFuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp,
                          bd->ndomelems, bd->domelind, domelem, domelcpind,
                          nncpi, bd->nncpi, false,
                          ftab, gtab, htab, &fga, bd->nvars, dcoeff,
                          bd->hsize, bd->hprof, bd->hrows ) )
        goto failure;
      bd->hblstat = 1;
    }
    for ( i = 0; i < nbl; i++ ) {
      bd = &d->block[i];
      memcpy ( bd->lhrows[0], bd->hrows[0], bd->hsize*sizeof(double) );
      if ( !pkn_NRBSymCholeskyDecompd ( bd->nvars, bd->hprof,
                                        bd->lhrows[0], bd->lhrows, NULL ) ) {
#ifdef _DEBUG
printf ( "!" );
#endif
        d->next_entire = d->itn + /*2**/nbl;
        goto try_with_one_block;
      }
      else {
        bd->hblstat = 3;
#ifdef _DEBUG
printf ( "+" );
#endif
      }
    }
        /* setup and decompose the coarse mesh preconditioner term */
    if ( d->phrows ) {
      if ( !_g2mbl_CMPSetupCoarseTerm ( d ) )
        goto failure;
#ifdef _DEBUG
printf ( "(" );
#endif
      if ( !pkn_NRBSymCholeskyDecompd ( 3*d->nwcp, d->phprof,
                                        d->phrows[0], d->phrows, NULL ) ) {
#ifdef _DEBUG
printf ( "!) " );
#endif
        goto try_with_one_block;
      }
#ifdef _DEBUG
printf ( "+)" );
#endif
    }

    if ( nkn != d->nkn1 ) {
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
      d->nG ++;
#endif
      g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                         ndomelems, domelem, domelcpind, nncpi,
                         true, ftab, gtab, &func, nvars, grad );
    }
    d->func = func;
        /* now the actual Newton method */
    gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
    for ( i = 0; i < MAXNTN; i++ ) {
#ifdef _DEBUG
printf ( ";" );
#endif
      memset ( dcoeff, 0, nvars*sizeof(double) );
      pkn_PCGd ( nvars, d, grad, dcoeff, _g2mbl_MultAxd, _g2mbl_MultQixAltd,
                 min ( nvars, MAXCG ), EPS3*gn*gn, EPS4*gn*gn, &itm );
#ifdef _DEBUG
printf ( " %d ", itm );
#endif
      ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
      _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
      g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                         ndomelems, domelem, domelcpind, nncpi,
                         true, ftab, gtab, &func, nvars, grad );
      if ( func < d->func ) {
        memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
        if ( ldco <= EPS6*lco && d->func-func < EPS7*func )
          SUCCESS_CG;
        d->func = func;
        if ( ldco <= EPS1*lco )
          SUCCESS_CG;
        ga = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        if ( ga > THR4*gn ) {
#ifdef _DEBUG
printf ( ";" );
#endif
          memset ( dcoeff, 0, nvars*sizeof(double) );
          pkn_PCGd ( nvars, d, grad, dcoeff, _g2mbl_MultAxd, _g2mbl_MultQixAltd,
                     min ( nvars, MAXCG ), EPS3*gn*gn, EPS4*gn*gn, &itm );
#ifdef _DEBUG
printf ( " %d ", itm );
#endif
          ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
          _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "F" );
#endif
          func = g2mbl_UFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                                ndomelems, domelem, domelcpind, true, ftab );
          if ( func < d->func ) {
            memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
            d->func = func;
            if ( ldco <= EPS1*lco )
              SUCCESS_CG;
          }
          goto next_iter;
        }
        gn = ga;
      }
      else if ( i == 0 ) {
        d->next_entire = d->itn + /*2**/nbl;
        goto try_with_one_block;
      }
      else
        break;
    }
  }
  else {
    hsize = d->hsize;
    hprof = d->hprof;
    Hessian = d->Hessian;
    hrows = d->hrows;
#ifdef _DEBUG
printf ( "H*" );
#endif
#ifdef COUNT
d->nH ++; 
#endif
    if ( !g2mbl_UFuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                          nv, mvcp, ndomelems, domelem, domelcpind, nncpi,
                          true, ftab, gtab, htab,
                          &func, nvars, grad, hsize, hprof, hrows ) )
      goto failure;
        /* no need to copy the Hessian */
    if ( !pkn_NRBSymCholeskyDecompd ( nvars, hprof, Hessian, hrows, NULL ) ) {
      for ( i = 0; i < nbl; i++ ) {
        bd = &d->block[i];
        if ( !g2mbl_AFuncGradHessiand ( d->nkn1, d->aqcoeff,
                      d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp,
                      bd->ndomelems, bd->domelind, domelem, domelcpind,
                      nncpi, bd->nncpi, false,
                      ftab, gtab, htab, &fga, bd->nvars, dcoeff,
                      bd->hsize, bd->hprof, bd->hrows ) )
          goto failure;
        bd->hblstat = 1;
      }
      d->next_entire = d->itn + /*2**/nbl;
      goto try_with_one_block;
    }
    if ( nkn != d->nkn1 ) {
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
      d->nG ++;
#endif
      g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                         ndomelems, domelem, domelcpind, nncpi,
                         true, ftab, gtab, &func, nvars, grad );
    }
    d->func = func;
        /* now the actual Newton method */
    gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
    for ( i = 0; i < MAXNTN; i++ ) {
#ifdef _DEBUG
printf ( "," );
#endif
      pkn_NRBLowerTrSolved ( nvars, hprof, Hessian, hrows, 1, 1, grad, 1, dcoeff );
      pkn_NRBUpperTrSolved ( nvars, hprof, Hessian, hrows, 1, 1, dcoeff, 1, dcoeff );
      ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
      _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
      g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                         ndomelems, domelem, domelcpind, nncpi,
                         true, ftab, gtab, &func, nvars, grad );
      if ( func < d->func ) {
        memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
        if ( ldco <= EPS6*lco && d->func-func < EPS7*func )
          SUCCESS_CG;
        d->func = func;
        if ( ldco <= EPS1*lco )
          SUCCESS_CG;
        ga = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        if ( ga > THR4*gn ) {
#ifdef _DEBUG
printf ( "," );
#endif
          pkn_NRBLowerTrSolved ( nvars, hprof, Hessian, hrows, 1, 1, grad, 1, dcoeff );
          pkn_NRBUpperTrSolved ( nvars, hprof, Hessian, hrows, 1, 1, dcoeff, 1, dcoeff );
          ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
          _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "F" );
#endif
          func = g2mbl_UFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                                ndomelems, domelem, domelcpind, true, ftab );
          if ( func < d->func ) {
            memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
            d->func = func;
            if ( ldco <= EPS1*lco )
              SUCCESS_CG;
          }
          goto next_iter;
        }
        gn = ga;
      }
      else if ( i == 0 ) {
        for ( i = 0; i < nbl; i++ ) {
          bd = &d->block[i];
          if ( !g2mbl_AFuncGradHessiand ( d->nkn1, d->aqcoeff,
                        d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp,
                        bd->ndomelems, bd->domelind, domelem, domelcpind,
                        nncpi, bd->nncpi, false,
                        ftab, gtab, htab, &fga, bd->nvars, dcoeff,
                        bd->hsize, bd->hprof, bd->hrows ) )
            goto failure;
          bd->hblstat = 1;
        }
        d->next_entire = d->itn + /*2**/nbl;
        goto try_with_one_block;
      }
      else
        break;
    }
  }
  goto next_iter;

try_with_one_block:
  pkv_SetScratchMemTop ( sp );
  d->last_ntn = false;
          /* determine the block for the current iteration */
  if ( d->itn > 0 ) {
            /* the idea commented out below does not work very well */
/*
    i = d->nextblock;
    if ( d->block[i].hblstat == 3)
*/
      i = d->nextblock = (d->nextblock+1) % nbl;
  }
  else
    i = d->nextblock = 0;
#ifdef _DEBUG
printf ( "(block %d):", d->nextblock );
#endif
  bd = &d->block[i];
  bdomelems = bd->ndomelems;
  domelind = bd->domelind;
  nbcp = bd->ncp;
  nvars = bd->nvars;
  nncpi = bd->nncpi;
  vncpi = bd->vncpi;
  hsize = bd->hsize;
  hprof = bd->hprof;
  hrows = bd->hrows;
  Hessian = hrows[0];
  lhrows = bd->lhrows;
  LHessian = lhrows[0];
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvars );
  if ( !coeff )
    goto failure;
  dcoeff = &coeff[nvars];
  grad = &dcoeff[nvars];
  if ( d->accurate ) {
    nkn    = d->nkn2;
    Nitabs = d->bNitabs;
    Jac    = d->bJac;
    qcoeff = d->bqcoeff;
  }
  else {
    nkn    = d->nkn1;
    Nitabs = d->aNitabs;
    Jac    = d->aJac;
    qcoeff = d->aqcoeff;
  }
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp )
    goto failure;
  memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
  lco = 0.0;
  for ( i = 0; i < nbcp; i++ )
    lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
  lco = sqrt ( lco );

  if ( bd->hblstat == 0 ) {
#ifdef _DEBUG
printf ( "H" );
#endif
#ifdef COUNT
  d->nH ++;
#endif
    if ( !g2mbl_AFuncGradHessiand ( d->nkn1, d->aqcoeff,
                         d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                         nv, mvcp, bdomelems, domelind, domelem, domelcpind,
                         d->nncpi, nncpi, true, ftab, gtab, htab,
                         &func, nvars, grad, hsize, hprof, hrows ) )
      goto failure;
    bd->hblstat = 1;
    d->newpoint = false;
    if ( nkn != d->nkn1 ) {
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
      d->nG ++;
#endif
      if ( !g2mbl_AFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                         bdomelems, domelind, domelem, domelcpind, bd->nncpi,
                         true, ftab, gtab, &func, nvars, grad ) )
        goto failure;
    }
  }
  else {  /* try to reuse the Hessian block from the previous iteration */
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
    d->nG ++;
#endif
    if ( !g2mbl_AFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                       bdomelems, domelind, domelem, domelcpind, bd->nncpi,
                       true, ftab, gtab, &func, nvars, grad ) )
      goto failure;
  }
  bd->func = func;
  gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
  if ( d->itn < nbl )
    bd->gn0 = gn;
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
        /* first, try the Newton method */
  if ( bd->hblstat >= 2 )
    positive = true;
  else {
    memcpy ( LHessian, Hessian, hsize*sizeof(double) );
    positive = pkn_NRBSymCholeskyDecompd ( nvars, hprof, LHessian, lhrows, NULL );
    if ( positive )
      bd->hblstat = 2;
#ifdef _DEBUG
if ( positive ) printf ( "+" );
#endif
  }
  if ( positive ) {
    pkn_NRBLowerTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, dcoeff, 1, dcoeff );
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    _g2mbl_AddCPIncrement ( nbcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );

#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
    d->nG ++;
#endif
    if ( !g2mbl_AFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                       bdomelems, domelind, domelem, domelcpind, bd->nncpi,
                       true, ftab, gtab, &fga, nvars, coeff ) )
      goto failure;
    ga = sqrt ( pkn_ScalarProductd ( nvars, coeff, coeff ) );
    if ( fga >= bd->func )
      goto switch_to_lmt;
    bd->func = fga;
    memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
    if ( ga < THR3*gn )
      bd->hblstat = 3;
    else
      bd->hblstat = 0;
#ifdef _DEBUG
printf ( "," );
#endif
    pkn_NRBLowerTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, coeff, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, dcoeff, 1, dcoeff );
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    _g2mbl_AddCPIncrement ( nbcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );

#ifdef _DEBUG
printf ( "F" );
#endif
#ifdef COUNT
    d->nF ++;
#endif
    func = g2mbl_AFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                          bdomelems, domelind, domelem, domelcpind, true, ftab );
    if ( func >= bd->func ) {
      bd->hblstat = 0;
      goto next_iter;
    }
#ifdef COUNT
    d->nN ++;
#endif
    memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
    bd->func = func;
    d->newpoint = true;
    if ( !bd->accurate ) {  /* switch accuracy */
#ifdef _DEBUG
printf ( "*" );
#endif
      bd->accurate = true;
      d->ibl --;
      if ( d->ibl == 0 ) {
        for ( k = 0; k < nbl; k++ )       /* this is to prevent effects of changing */
          d->block[k].func = MYINFINITY;  /* quadrature approximation error */
        d->accurate = true;
#ifdef _DEBUG
printf ( "*" );
#endif
      }
    }
    goto next_iter;
  }
#ifdef _DEBUG
else
  printf ( "!" );
#endif
switch_to_lmt:
#ifdef COUNT
  d->nLM ++;
#endif
  bd->hblstat = 0;

        /* minimization along the Levenberg-Marquardt trajectory */
  ga = 0.0;  fga = MYINFINITY;
  gb = ge = MYINFINITY;
  fgb = fge = func;
  memset ( coeff, 0, nvars*sizeof(double) );
          /* get the initial nu for bracketing */
  if ( bd->nu <= 0.0 ) {
    if ( !pkn_NRBSymFindEigenvalueIntervald ( nvars, hprof, Hessian, hrows,
                                              &lmin, &lmax ) )
      goto failure;
    bd->nu = gc = ALPHA*(lmax-lmin);
  }
  else
    gc = bd->nu;
  fgc = MLFUNC ( gc );
  RECORD_MIN ( gc, fgc );
        /* bracketing */
  for ( i = 0; fgc >= fga; i++ ) {
    ga = gc;  fga = fgc;
    gc *= (double)(i+2);  /* trying the Fibonacci sequence */
    fgc = MLFUNC ( gc );
    RECORD_MIN ( gc, fgc );
  }
  gd = gc*(double)(i+2);
  fgd = MLFUNC ( gd );
  RECORD_MIN ( gd, fgd );
  if ( fgd < fgc ) {
    do {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = gc*(double)(i+2);  i++;
      fgd = MLFUNC ( gd );
      RECORD_MIN ( gd, fgd );
    } while ( fgd < fgc );
    gb = gd;  fgb = fgd;
  }
  else {
    gb = gd;  fgb = fgd;
    for ( i = 0; ; i++ ) {
      if ( ga > 0.0 )
        gd = sqrt ( ga*gc );
      else if ( gc > 1.0e4 )
        gd = sqrt ( gc );
      else
        gd = 0.01*gc;
      fgd = MLFUNC ( gd );
      RECORD_MIN ( gd, fgd );
      if ( i >= MAXGTN && fge < func )
        goto finish_it;
      if ( fgd >= fga ) {
        ga = gd;  fga = fgd;
      }
      else if ( fgd < fgc ) {
        gb = gc;  fgb = fgc;  gc = gd;  fgc = fgd;
      }
      else {
        ga = gd;  fga = fgd;
        break;
      }
    }
  }
        /* find second point in the bracket to start the golden ratio division */
  if ( fga < fgb ) {
    gd = gc;  fgd = fgc;
    gc = GRDIV ( ga, gd );
    fgc = MLFUNC ( gc );
    RECORD_MIN ( gc, fgc );
  }
  else {
    gd = GRDIV ( gb, gc );
    fgd = MLFUNC ( gd );
    RECORD_MIN ( gd, fgd );
  }
        /* golden ratio division */
#ifdef _DEBUG
printf ( "g" );
#endif
  for ( i = 0; i < MAXGTN; i++ ) {
/*
    if ( fgc < fgd ) {
      gb = gd;  fgb = fgd;
      gd = gc;  fgd = fgc;
      gc = GRDIV ( ga, gd );
      fgc = MLFUNC ( gc );
      RECORD_MIN ( gc, fgc )
    }
    else {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = GRDIV ( gb, gc );
      fgd = MLFUNC ( gd );
      RECORD_MIN ( gd, fgd )
    }
*/
    if ( _g2mbl_DivideIntervald ( &ga, &gc, &gd, &gb,
                                  &fga, &fgc, &fgd, &fgb ) ) {
      fgc = MLFUNC ( gc );
      RECORD_MIN ( gc, fgc )
    }
    else {
      fgd = MLFUNC ( gd );
      RECORD_MIN ( gd, fgd )
    }
    if ( fga - fge < GRTHR*(func-fge) ||
         fga - fge < EPS2*fge )
      break;
  }
finish_it:
  _g2mbl_AddCPIncrement ( nbcp, vncpi, mvcp, (vector3d*)coeff, mvcp  );

next_iter:
  d->itn ++;
#ifdef _DEBUG
printf ( "\n" );
#endif
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_IterBlSurfaceOptAltBLMTd*/

boolean g2mbl_FindBlSurfaceAltBLMTd ( int nv, BSMvertex *mv, int *mvhei, 
                                     point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                     int nfac, BSMfacet *mfac, int *mfhei,
                                     byte *mkcp,
                                     double C, double dO, double dM,
                                     int maxit, int nkn1, int nkn2, int nbl )
{
  void             *sp;
  mesh_lmt_optdata *data;
  int              itn;
  boolean          finished;

  sp = pkv_GetScratchMemTop ();
  data = NULL;
  if ( !g2mbl_InitBlSurfaceOptAltBLMTd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                      nfac, mfac, mfhei, mkcp,
                                      C, dO, dM, nkn1, nkn2, nbl,
                                      (void**)&data ) )
    goto failure;

  finished = false;
  for ( itn = 0; itn < maxit; itn ++ ) {
    if ( !g2mbl_IterBlSurfaceOptAltBLMTd ( (void*)data, &finished ) )
      goto failure;
    if ( finished )
      break;
  }

#ifdef COUNT
printf ( "nF = %d, nG = %d, nH = %d, nN = %d, nLM = %d\n",
         data->nF, data->nG, data->nH, data->nN, data->nLM );
#endif

  g2mbl_OptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_OptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_FindBlSurfaceAltBLMTd*/

