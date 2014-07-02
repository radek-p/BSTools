
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
#include "raybez.h"
#include "mengerc.h"

#include "mengercprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
static boolean MHEig ( int N, void *usrdata, double *x, double *f )
{
  void         *sp;
  mengerc_data *md;
  int          n, nn, i;
  double       *h, *eigval, emin, emax, ef;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;
  n = 3*(md->lkn - 2*md->deg);
  nn = (n*(n+1))/2;
  eigval = pkv_GetScratchMemd ( n+nn );
  if ( !eigval )
    goto failure;
  h = &eigval[n];
  if ( !md->heigok ) {
    memcpy ( h, md->hhkMe, nn*sizeof(double) );
    if ( pkn_SymMatFindEigenvaluesd ( n, h, eigval ) ) {
      emin = emax = eigval[0];
      for ( i = 1; i < n; i++ )
        if ( eigval[i] < emin ) emin = eigval[i];
        else if ( eigval[i] > emax ) emax = eigval[i];
      md->hmin = emin;
      md->hmax = emax;
    }  
    md->heigok = true;
  }
  pkn_AddMatrixMd ( 1, nn, 0, md->hhkMe, 0, md->hhR1, x[0], 0, h );
  pkn_AddMatrixMd ( 1, nn, 0, h, 0, md->hhR2, x[1], 0, h );
  pkn_AddMatrixMd ( 1, nn, 0, h, 0, md->hhR3, x[2], 0, h );
  pkn_AddMatrixMd ( 1, nn, 0, h, 0, md->hhR4, x[3], 0, h );
  pkn_AddMatrixMd ( 1, nn, 0, h, 0, md->hhR5, x[4], 0, h );
  if ( pkn_SymMatFindEigenvaluesd ( n, h, eigval ) ) {
    emin = emax = eigval[0];
    for ( i = 1; i < n; i++ )
      if ( eigval[i] < emin ) emin = eigval[i];
      else if ( eigval[i] > emax ) emax = eigval[i];
/*    *f = ef = 1.0e-3*emax - emin; */
/*    *f = ef = (1.0e-9*emax + 1.0e-3)*emax - emin; */
    md->_emin = emin;
    md->_emax = emax;

    ef = (1.0e-9*emax + 1.0e-3)*emax - (0.01*emin*emin+1.0)*emin;
    if ( emax > md->hmax )
      ef += 1.0e-9*(emax-md->hmax)*(emax-md->hmax)*(emax-md->hmax);
    if ( ef < md->ef ) {
      md->emin = emin;
      md->emax = emax;
      md->ef = ef;
    }
    *f = md->f = ef;
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  *f = 1.0e308; /* infinity */
  pkv_SetScratchMemTop ( sp );
  return false;
} /*MHEig*/

#define EPS 1.0e-3
static boolean MHEigGrad ( int N, void *usrdata,
                           double *x, double *f, double *grad )
{
  void   *sp;
  int    i;
  double _f, h, s;

  sp = pkv_GetScratchMemTop ();
  if ( !MHEig ( 5, usrdata, x, &_f ) )
    goto failure;
  *f = _f;
  for ( i = 0; i < 5; i++ ) {
    s = x[i];
    h = EPS*s;
    x[i] += h;
    if ( !MHEig ( 5, usrdata, x, &grad[i] ) )
      goto failure;
    x[i] = s-h;
    if ( !MHEig ( 5, usrdata, x, &_f ) )
      goto failure;
    grad[i] = (grad[i]-_f)/(h+h);
    x[i] = s;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*MHEigGrad*/

static boolean MHEigGH ( int N, void *usrdata,
                         double *x, double *f, double *g, double *h )
{
  int    i, j;
  double ht[5], fa[5], fb[5], _f, s, t, d, e, ha, hb;

  if ( !MHEig ( 5, usrdata, x, &_f ) )
    goto failure;
  *f = _f;
  for ( i = 0; i < 5; i++ ) {
    ht[i] = d = EPS*x[i];
    s = x[i];
    x[i] = s+d;
    if ( !MHEig ( 5, usrdata, x, &fa[i] ) )
      goto failure;
    x[i] = s-d;
    if ( !MHEig ( 5, usrdata, x, &fb[i] ) )
      goto failure;
    g[i] = (fa[i]-fb[i])/(d+d);
    h[pkn_LowerTrMatIndex(i,i)] = (fa[i]+fb[i]-_f-_f)/(d*d);
    for ( j = 0; j < i; j++ ) {
      e = ht[j];
      t = x[j];
      x[i] = s+d;
      x[j] = t+e;
      if ( !MHEig ( 5, usrdata, x, &ha ) )
        goto failure;
      x[i] = s-d;
      x[j] = t-e;
      if ( !MHEig ( 5, usrdata, x, &hb ) )
        goto failure;
      h[pkn_LowerTrMatIndex(i,j)] = (ha+hb+_f+_f-fa[i]-fb[i]-fa[j]-fb[j])/(2.0*d*e);
      x[j] = t;
    }
    x[i] = s;
  }
  return true;

failure:
  return false;
} /*MHEigGH*/

static boolean MHEigTunnel ( int N, void *usrdata,
                             double *x0, double *x1, boolean *went_out )
{
  *went_out = false;
  if ( x1[0] < 0.01   ) *went_out = true;
  if ( x1[1] < 1.0    )  *went_out = true;
  if ( x1[2] < 1.0e+2 )  *went_out = true;
  if ( x1[3] < 1.0    )  *went_out = true;
  if ( x1[4] < 1.0    )  *went_out = true;
  return true;
} /*MHEigTunnel*/

static int Qq ( char *format, double hmin, double hmax )
{
  return printf ( format, hmin, hmax );
} /*Qq*/

boolean mengerc_OptPenaltyParams1 ( mengerc_data *md, boolean wide )
{
/*#define MAXITER 20*/
#define MAXITER 10
  void      *sp;
  int       n, nn;
  double    *penalty_param, nu;
  double    x[5], y[5], z[5], emin, emax, func, f, *g, *h;
  double    fct1[7] = { 1.0, 0.5, 2.0, 0.25, 4.0, 0.125, 8.0 };
  int       limits1[5] = {7,7,7,7,7};
  double    fct2[5] = {1.0, 0.25, 4.0, 0.0625, 16.0};
  int       limits2[5] = {5,5,5,5,5};
  double    *fct;
  int       *limits;
  int       cnt[5], cntmin[5], i;
  boolean   went_out, progress;

  sp = pkv_GetScratchMemTop ();
  progress = false;
  n = 3*(md->lkn - 2*md->deg);
  nn = (n*(n+1))/2;
  g = pkv_GetScratchMemd ( n+nn );
  if ( !g )
    goto failure;
  h = &g[n];
  penalty_param = md->penalty_param;
  if ( !mengerc_IntegralMengerfgh ( n, (void*)md, &md->cpoints[0].x,
                                    &f, g, h ) )
    goto failure;
  md->heigok = false;

/* 1st stage - brute force search for an initial point */
  if ( wide ) {
    fct = fct1;
    limits = limits1;
  }
  else {
    fct = fct2;
    limits = limits2;
  }
  memcpy ( y, penalty_param, 5*sizeof(double) );
  emin = -1.0e-300;
  emax = func = 1.0e+308;
  memset ( cnt, 0, 5*sizeof(int) );
  do {
    for ( i = 0; i < 5; i++ )
      x[i] = y[i]*fct[cnt[i]];
    MHEigTunnel ( 5, NULL, penalty_param, x, &went_out );
    if ( went_out )
      continue;
    if ( !MHEig ( 5, (void*)md, x, &f ) )
      goto second_st;
    if ( f < func ) {
      func = f;
      emin = md->emin = md->_emin;
      emax = md->emax = md->_emax;
      memcpy ( z, x, 5*sizeof(double) );
      memcpy ( cntmin, cnt, 5*sizeof(int) );
    }
  } while ( pkv_IncMultiCounter ( 5, limits, cnt ) );
  memcpy ( penalty_param, z, 5*sizeof(double) );
  for ( i = 0; i < 5; i++ )
    if ( cntmin[i] ) {
      progress = true;
      break;
    }

  Qq ( "hmin = %10.5g, hmax = %10.5g\n", md->hmin, md->hmax );
  printf ( "emin = %10.5g, emax = %10.5g, ef = %10.5g\n", emin, emax, func );
  for ( i = 0; i < 5; i++ )
    printf ( "%10.5g, %2d\n", z[i], cntmin[i] );

/* 2nd stage - minimization along Levenberg-Marquardt trajectories */
second_st:
  md->ef = 1.0e+308;
  nu = -1.0;
  for ( i = 0; i < MAXITER; i++ ) {
    switch ( pkn_NLMIterd ( 5, (void*)md, penalty_param,
                            MHEig, MHEigGrad, MHEigGH, NULL, MHEigTunnel,
                            -1.0e10, 1.0e-5, 1.0e-5, &nu ) ) {
  case PKN_LMT_ERROR:
      goto failure;
  case PKN_LMT_CONTINUE_N:
  case PKN_LMT_CONTINUE_LM:
  case PKN_LMT_CONTINUE_LM_P:
      break;
  case PKN_LMT_FOUND_MINIMUM:
  case PKN_LMT_FOUND_ZEROGRAD:
  case PKN_LMT_FOUND_ZEROGRAD_P:
  case PKN_LMT_NO_PROGRESS:
  case PKN_LMT_CROSSED_LIMIT:
      goto way_out;
  default:
      goto failure;
    }
  }
way_out:
  printf ( ": emin = %10.5g, emax = %10.5g, ef = %10.5g\n",
           md->emin, md->emax, md->ef );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return progress;
} /*mengerc_OptPenaltyParams1*/

/* ////////////////////////////////////////////////////////////////////////// */
static boolean SDFunc ( int N, void *usrdata, double *x, double *f )
{
#define GFACT 1.0e-3
  void         *sp;
  mengerc_data *md;
  double       _f, gn, *grad, kk[5], nu, *scp;
  int          ncp, n;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;
  ncp = md->lkn-md->deg;
  n = 3*(ncp-md->deg);
  memcpy ( kk, md->penalty_param, 5*sizeof(double) );
  scp = pkv_GetScratchMemd ( N+3*ncp );
  if ( !scp )
    goto failure;
  grad = &scp[3*ncp];
  memcpy ( scp, md->cpoints, ncp*sizeof(point3d) );
  memcpy ( md->penalty_param, x, 5*sizeof(double) );

  nu = -1.0;
  switch ( pkn_NLMIterd ( n, usrdata, &md->cpoints[0].x, mengerc_IntegralMengerf,
                          mengerc_IntegralMengerfg, mengerc_IntegralMengerfgh,
                          mengerc_IntegralMengerTransC, mengerc_HomotopyTest,
                          0.0, 1.0e-4, 1.0e-7, &nu ) ) {
case PKN_LMT_CONTINUE_N:
case PKN_LMT_CONTINUE_LM:
case PKN_LMT_CONTINUE_LM_P:
case PKN_LMT_FOUND_MINIMUM:
case PKN_LMT_FOUND_ZEROGRAD:
case PKN_LMT_FOUND_ZEROGRAD_P:
    if ( !mengerc_IntegralMengerfg ( N, usrdata, &md->cpoints[0].x,
                                     &_f, grad ) )
      goto failure;
    gn = sqrt ( pkn_ScalarProductd ( N, grad, grad ) );
    *f = _f += GFACT*gn;
    break;

case PKN_LMT_FOUND_BARRIER:
    *f = _f = 1.0e+308;
    break;

case PKN_LMT_ERROR:
case PKN_LMT_CROSSED_LIMIT:
case PKN_LMT_NO_PROGRESS:
default:
    goto failure;
  }
  if ( _f < md->fmin ) {
    memcpy ( md->mcpoints, md->cpoints, n*sizeof(double) );
    md->fmin = _f;
  }

  memcpy ( md->penalty_param, kk, 5*sizeof(double) );
  memcpy ( md->cpoints, scp, ncp*sizeof(point3d) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( scp )
    memcpy ( md->cpoints, scp, ncp*sizeof(point3d) );
  memcpy ( md->penalty_param, kk, 5*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return false;
#undef GFACT
} /*SDFunc*/

static boolean SDFuncGrad ( int N, void *usrdata, double *x, double *f, double *g )
{
  void   *sp;
  int    i;
  double fa[5], fb[5], _f, s, d;

  sp = pkv_GetScratchMemTop ();

  if ( !SDFunc ( 5, usrdata, x, &_f ) )
    goto failure;
  *f = _f;
  for ( i = 0; i < 5; i++ ) {
    d = EPS*x[i];
    s = x[i];
    x[i] = s+d;
    if ( !SDFunc ( 5, usrdata, x, &fa[i] ) )
      goto failure;
    x[i] = s-d;
    if ( !SDFunc ( 5, usrdata, x, &fb[i] ) )
      goto failure;
    g[i] = (fa[i]-fb[i])/(d+d);
    x[i] = s;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*SDFuncGrad*/

boolean mengerc_OptPenaltyParams2 ( mengerc_data *md )
{
/*#define MAXIT 20*/
#define MAXIT 20
  void   *sp;
  int    n, i;
  double *penalty_param, nu;
  double f, *g, _f;

  sp = pkv_GetScratchMemTop ();
  n = 3*(md->lkn - 2*md->deg);
  g = pkv_GetScratchMemd ( n );
  if ( !g )
    goto failure;
  penalty_param = md->penalty_param;
  if ( !mengerc_IntegralMengerfg ( n, (void*)md, &md->cpoints[0].x, &f, g ) )
    goto failure;
  md->heigok = false;
  memset ( md->mcpoints, 0, n*sizeof(double) );
  md->fmin = 1.0e+308;

  md->ef = 1.0e+308;
  nu = -1.0;
  for ( i = 0; i < MAXIT; i++ ) {
    switch ( pkn_SDIterd ( 5, (void*)md, penalty_param,
                            SDFunc, SDFuncGrad, NULL, MHEigTunnel,
                            -1.0e10, 1.0e-5, 1.0e-5, &nu ) ) {
  case PKN_SD_CONTINUE:
      break;
  case PKN_SD_FOUND_BARRIER:
  case PKN_SD_FOUND_ZEROGRAD:
  case PKN_SD_NO_PROGRESS:
  case PKN_SD_CROSSED_LIMIT:
      goto way_out;
  case PKN_SD_ERROR:
  default:
      goto failure;
    }
  }

way_out:
  MHEig ( 5, (void*)md, penalty_param, &_f  );
  printf ( "# emin = %10.5g, emax = %10.5g, ef = %10.5g\n",
           md->_emin, md->_emax, md->f );
  if ( md->fmin < _f ) {
    memcpy ( md->cpoints, md->mcpoints, n*sizeof(double) );
    memcpy ( &md->cpoints[md->lkn-2*md->deg], md->cpoints, 3*md->deg );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_OptPenaltyParams2*/

/* ////////////////////////////////////////////////////////////////////////// */
static double _mengerc_opt3 ( void *usrptr, double lf )
{
  mengerc_data *md;
  double       f, penalty_param[5];

  md = (mengerc_data*)usrptr;
  pkn_MultMatrixNumd ( 1, 5, 0, md->penalty_param, exp(lf), 0, penalty_param );
  MHEig ( 5, usrptr, penalty_param, &f );
  return f;
} /*_mengerc_opt3*/

boolean mengerc_OptPenaltyParams3 ( mengerc_data *md )
{
  void    *sp;
  int     n, nn;
  double  lf, f, *g, *h;
  boolean error;

  sp = pkv_GetScratchMemTop ();
  n = 3*(md->lkn - 2*md->deg);
  nn = (n*(n+1))/2;
  g = pkv_GetScratchMemd ( n+nn );
  if ( !g )
    goto failure;
  h = &g[n];
  if ( !mengerc_IntegralMengerfgh ( n, (void*)md, &md->cpoints[0].x,
                                    &f, g, h ) )
    goto failure;
  md->heigok = false;
  lf = pkn_GoldenRatd ( _mengerc_opt3, (void*)md,
                        log (0.01), log (100.0), 1.0e-3, &error );
  if ( !error ) {
    pkn_MultMatrixNumd ( 1, 5, 0, md->penalty_param, exp(lf), 0, md->penalty_param );
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  else {
failure:
    pkv_SetScratchMemTop ( sp );
    return false;
  }
} /*mengerc_OptPenaltyParams3*/

