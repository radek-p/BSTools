
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
#include "msgpool.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define EPS0    5.0e-10
#define EPS1    1.0e-8 
#define EPS2    1.0e-7
#define MAXNTN  16
#define MAXGTN  16
#define MAXSTN  30
#define MAXBTN  10
#define MAXCTN  20
#define THR1    0.9
#define THR2    0.75
#define THR3    0.25
#define GRTHR   0.02

/* ///////////////////////////////////////////////////////////////////////// */
double _g2mbl_AuxNuFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int nvcp, int *vncpi,
                           int ndomelems, meshdom_elem *domelem, int *domelcpind,
                           double *ftab, int nvars, int hsize, int *hprof,
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
    return g2mbl_UFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                          ndomelems, domelem, domelcpind, true, ftab );
  }
  else
    return MYINFINITY;
} /*_g2mbl_AuxNuFuncd*/

boolean g2mbl_IterBlSurfaceOptLMTd ( void *data, boolean *finished )
{
  void             *sp;
  mesh_lmt_optdata *d;
  int              nv;
  point3d          *mvcp, *auxmvcp;
  int              ndomelems, *domelcpind;
  int              nvcp, nvars, *nncpi, *vncpi, nkn;
  meshdom_elem     *domelem;
  double           *ftab, *gtab, *htab, **Nitabs, **Jac, *qcoeff;
  int              hsize, *hprof;
  double           *Hessian, **hrows, *LHessian, **lhrows;
  double           *coeff, *dcoeff, *grad;
  double           func, gn, gn1, lco, ldco, df, dq, rk;
  boolean          all, positive;
  int              i, ktn;
  double           lmin, lmax, ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge;

#ifdef _DEBUG
#define SUCCESS \
  { \
printf ( "\nfunc = %12.7e, gn = %12.7e", func, gn ); \
    *finished = true; \
    goto next_iter; \
  }
#else
#define SUCCESS \
  { \
    *finished = true; \
    goto next_iter; \
  }
#endif

#ifdef _DEBUG
#define SWITCH_ACCURACY \
  { \
printf ( "*" ); \
    d->accurate = d->newpoint = true; \
    d->func = MYINFINITY;  /* this is to prevent effects of changing */ \
                           /* quadrature approximation error */ \
    goto next_iter; \
  }
#else
#define SWITCH_ACCURACY \
  { \
    d->accurate = d->newpoint = true; \
    d->func = MYINFINITY;  /* this is to prevent effects of changing */ \
                           /* quadrature approximation error */ \
    goto next_iter; \
  }
#endif

#define MLFUNC(nu) \
_g2mbl_AuxNuFuncd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, nvcp, vncpi, \
        ndomelems, domelem, domelcpind, ftab, nvars, hsize, hprof, \
        Hessian, LHessian, lhrows, grad, dcoeff, auxmvcp, nu )

#define RECORD_MIN(g,fnu) \
  { if ( fnu < fge ) { \
      memcpy ( coeff, dcoeff, nvars*sizeof(double) ); \
      d->func = fge = fnu; \
      d->nu[0] = d->nu[1] = ge = g; \
      d->newpoint = true; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
  *finished = false;
        /*extract data to local variables */
  d = data;
#ifdef _DEBUG
printf ( "%2d: ", d->itn );
#endif
  fga = fgc = func = d->func;
  gn = d->gn0;
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
  Hessian = d->Hessian;
  hrows = d->hrows;
  LHessian = d->LHessian;
  lhrows = d->lhrows;
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvars );
  if ( !coeff )
    goto failure;
  dcoeff = &coeff[nvars];
  grad = &dcoeff[nvars];
  if ( d->accurate ) {
    nkn = d->nkn2;
    Nitabs = d->bNitabs;
    Jac    = d->bJac;
    qcoeff = d->bqcoeff;
  }
  else {
    nkn = d->nkn1;
    Nitabs = d->aNitabs;
    Jac    = d->aJac;
    qcoeff = d->aqcoeff;
  }
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp )
    goto failure;
  memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
  lco = 0.0;
  for ( i = 0; i < nvcp; i++ )
    lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
  lco = sqrt ( lco );

  all = true;
  if ( d->newpoint ) {
#ifdef _DEBUG
printf ( "H" );
#endif
#ifdef COUNT
  d->nH ++;
#endif
    if ( !g2mbl_UFuncGradHessiand ( d->nkn1, d->aqcoeff,
                          d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac,
                          nv, mvcp, ndomelems, domelem, domelcpind, nncpi,
                          all, ftab, gtab, htab,
                          &func, nvars, grad, hsize, hprof, hrows ) )
      goto failure;
    d->newpoint = false;
  }
  else {
    if ( d->accurate ) {
      SUCCESS
    }
    else
      SWITCH_ACCURACY
  }
  if ( nkn != d->nkn1 )
    g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                       ndomelems, domelem, domelcpind, nncpi,
                       true, ftab, gtab, &func, nvars, grad );

  d->func = min ( func, d->func );
  gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
  if ( !d->itn )
    d->gn0 = gn;
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
        /* first, try the Newton method */
  memcpy ( LHessian, Hessian, hsize*sizeof(double) );
  positive = pkn_NRBSymCholeskyDecompd ( nvars, hprof, LHessian, lhrows, NULL );
  if ( positive ) {
#ifdef _DEBUG
printf ( "+" );
#endif
    pkn_NRBLowerTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, dcoeff, 1, dcoeff );
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
    if ( !d->accurate ) {
      if ( ldco <= 5.0*EPS1*lco )
        SWITCH_ACCURACY
    }
    else if ( ldco <= 5.0*EPS1*lco )
      SUCCESS

#ifdef _DEBUG
printf ( "F" );
#endif
#ifdef COUNT
    d->nF ++;
#endif
    fga = g2mbl_UFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                          ndomelems, domelem, domelcpind, true, ftab );
    if ( fga >= d->func )
      goto switch_to_lmt;
#ifdef COUNT
    d->nN ++;
#endif
    memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
    df = fga - d->func;
    d->func = fga;
    d->newpoint = true;
        /* test if the functional value decrement is sufficient */
    if ( !_g2bl_ComputeDeltaQd ( nvars, hprof, hrows, grad, dcoeff, &dq ) )
      goto failure;
    rk = df/dq;
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
    d->nG ++;
#endif
    if ( !g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                             ndomelems, domelem, domelcpind, nncpi,
                             true, ftab, gtab, &fga, nvars, coeff ) )
      goto failure;
    gn1 = sqrt ( pkn_ScalarProductd ( nvars, coeff, coeff ) );
    if ( gn1 > gn )
      goto switch_to_lmt;
    if ( !d->accurate ) {
      if ( gn < 5.0*EPS0*d->gn0 )
        SWITCH_ACCURACY
    }
    else if ( gn <= EPS0*d->gn0 )
      SUCCESS

    memcpy ( grad, coeff, nvars*sizeof(double) );
    if ( rk > THR1 ) {
        /* make additional iterations of the Newton method with */
        /* the same Hessian matrix */
      for ( ktn = 0; ktn < MAXNTN; ktn++ ) {
#ifdef _DEBUG
printf ( "," );
#endif
        gn1 = gn;
        pkn_NRBLowerTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, grad, 1, dcoeff );
        pkn_NRBUpperTrSolved ( nvars, hprof, LHessian, lhrows, 1, 1, dcoeff, 1, dcoeff );
        _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
        d->nG ++;
#endif
        if ( !g2mbl_UFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                                 ndomelems, domelem, domelcpind, nncpi,
                                 true, ftab, gtab, &fga, nvars, grad ) )
          goto failure;
        gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        if ( fga < d->func ) {
          d->func = fga;
          memcpy ( mvcp, auxmvcp, nv*sizeof(point3d) );
          if ( !d->accurate ) {
            if ( gn <= 5.0*EPS0*d->gn0 )
              SWITCH_ACCURACY
          }
          else if ( gn <= EPS0*d->gn0 )
            SUCCESS
          ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
          if ( !d->accurate ) {
            if ( ldco <= 5.0*EPS1*lco )
              SWITCH_ACCURACY
          }
          else if ( ldco <= EPS1*lco )
            SUCCESS
        }
        else {
            /* there may be fga > d->func due to rounding errors, */
            /* hence a way out is necessary */
          if ( !d->accurate ) {
            if ( gn <= 5.0*EPS1*d->gn0 )
              SWITCH_ACCURACY
          }
          else if ( gn <= EPS0*d->gn0 )
            SUCCESS
        }
        if ( gn > THR3*gn1 )
          goto next_iter;
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
        /* minimization along the Levenberg-Marquardt trajectory */
  ga = 0.0;  fga = MYINFINITY;
  gb = ge = MYINFINITY;
  fgb = fge = func;
  memset ( coeff, 0, nvars*sizeof(double) );
          /* get the initial nu for bracketing */
  if ( d->nu[0] <= 0.0 ) {
    if ( !pkn_NRBSymFindEigenvalueIntervald ( nvars, hprof, hrows[0], hrows,
                                              &lmin, &lmax ) ) {
printf ( "%s\n", ERRMSG_23 );
      goto failure;
    }
    d->nu[0] = d->nu[1] = gc = 0.01*(lmax-lmin);
  }
  else
    gc = sqrt ( d->nu[0]*d->nu[1] );
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
      if ( i >=  MAXGTN && fge < func )
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
  _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)coeff, mvcp );

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
} /*g2mbl_IterBlSurfaceOptLMTd*/

boolean g2mbl_FindBlSurfaceLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                  point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  byte *mkcp,
                                  double C, double dO, double dM,
                                  int maxit, int nkn1, int nkn2 )
{
  void             *sp;
  mesh_lmt_optdata *data;
  int              itn;
  boolean          finished;

  sp = pkv_GetScratchMemTop ();
  data = NULL;
  if ( !g2mbl_InitBlSurfaceOptLMTd ( nv, mv, mvhei, mvcp, nhe, mhe,
                                     nfac, mfac, mfhei, mkcp,
                                     C, dO, dM, nkn1, nkn2, (void**)&data ) )
    goto failure;

  finished = false;
  for ( itn = 0; itn < maxit; itn ++ ) {
    if ( !g2mbl_IterBlSurfaceOptLMTd ( (void*)data, &finished ) )
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
} /*g2mbl_FindBlSurfaceLMTd*/

