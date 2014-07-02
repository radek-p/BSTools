
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
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pkgeom.h"
#include "pknum.h"
#include "multibs.h"
#include "mengerc.h"

#include "mengercprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean mengerc_InitMCOptimization ( int deg, int lkn, double *knots,
                                     point3d *cpoints, double w,
                                     double penalty_param[MENGERC_NPPARAM],
                                     int nqkn, int npthr, int opt,
                                     mengerc_data *md )
{
  void   *sp;
  double f;
  int    nvars;

  sp = pkv_GetScratchMemTop ();
  nvars = 3*(lkn-2*deg);
  if ( !mengerc_BindACurve ( md, deg, lkn, knots, cpoints, nqkn,
                             w, penalty_param, w >= 4.0 ) )
    goto failure;
  md->npthr = npthr;
  md->pp_opt = opt;

  memcpy ( md->x, cpoints, nvars*sizeof(double) );
  md->pretransf = true;
  if ( !mengerc_IntegralMengerTransC ( nvars, (void*)md, md->x ) )
    goto failure;

  md->ppopt = true;
  switch ( opt ) {
case MENGERC_OPT_FULL1:
    mengerc_OptPenaltyParams1 ( md, false/*true*/ );
    break;
case MENGERC_OPT_FULL2:
    mengerc_OptPenaltyParams2 ( md );
    break;
case MENGERC_OPT_PART:
    mengerc_OptPenaltyParams3 ( md );
    break;
case MENGERC_OPT_NONE:
default:
    md->ppopt = false;
    break;
  }
  mengerc_IntegralMengerfg ( nvars, (void*)md, md->x, &f, md->g );
  md->lastf = f;
  md->gn = sqrt ( pkn_ScalarProductd ( nvars, md->g, md->g ) );
  md->ppopt = false;
  md->nu = -1.0;
  md->pretransf = true;
  md->fcnt = md->itc = 0;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  mengerc_UntieTheCurve ( md );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_InitMCOptimization*/

boolean mengerc_IterMCOptimization ( mengerc_data *md, boolean *finished )
{
  int     nvars;
  double  f, gn, lastf, lastgn;
  boolean opt_res;

  nvars = md->nvars;
  lastf = md->lastf;
  lastgn = md->gn;
  md->itres = pkn_NLMIterd ( nvars, (void*)md, md->x,
                             mengerc_IntegralMengerf,
                             mengerc_IntegralMengerfg,
                             mengerc_IntegralMengerfgh,
                             mengerc_IntegralMengerTransC,
                             mengerc_HomotopyTest,
                             0.0, 1.0e-4, 1.0e-7, &md->nu );

  mengerc_IntegralMengerfg ( nvars, (void*)md, md->x, &f, md->g );
  gn = sqrt ( pkn_ScalarProductd ( nvars, md->g, md->g ) );
  md->ppopt = false;

  switch ( md->itres ) {
case PKN_LMT_ERROR:
case PKN_LMT_CROSSED_LIMIT:
    goto failure;

case PKN_LMT_CONTINUE_N:
    if ( md->fcnt > 0 )
      md->fcnt --;
    break;

case PKN_LMT_CONTINUE_LM:
case PKN_LMT_CONTINUE_LM_P:
    md->fcnt ++;
    break;

case PKN_LMT_FOUND_MINIMUM:
case PKN_LMT_FOUND_ZEROGRAD:
case PKN_LMT_FOUND_ZEROGRAD_P:
case PKN_LMT_FOUND_BARRIER:
    *finished = true;
    goto next_iter;

case PKN_LMT_NO_PROGRESS:
    if ( md->pp_opt != MENGERC_OPT_NONE && md->fcnt > 0 )
      goto opt_param;
    else {
      *finished = true;
      goto next_iter;
    }

default:
    goto failure;
  }

  md->mdi = mengerc_ModifyRemotestPoint ( md->lkn-md->deg, md->cpoints,
                                          &md->sc, md->mdi );
  if ( md->fcnt >= 5+(int)(sqrt((double)(md->itc))) ) {
opt_param:
    if ( f <= 0.9*lastf || gn <= 0.9*lastgn )
      goto next_iter;
    opt_res = false;
    md->ppopt = true;
    switch ( md->pp_opt ) {
  case MENGERC_OPT_FULL1:
      opt_res = mengerc_OptPenaltyParams1 ( md, false );
      break;
  case MENGERC_OPT_FULL2:
      opt_res |= mengerc_OptPenaltyParams2 ( md );
      break;
  case MENGERC_OPT_PART:
      opt_res |= mengerc_OptPenaltyParams3 ( md );
      break;
  case MENGERC_OPT_NONE:
      md->ppopt = false;
  default:
      break;
    }
    if ( opt_res )
      md->fcnt = 0;
    else
      md->fcnt /= 2;
  }
next_iter:
  md->lastf = f;
  md->gn = gn;
  md->itc ++;
  return true;

failure:
  return false;
} /*mengerc_IterMCOptimization*/

boolean mengerc_OptimizeMengerCurvature (
                      int deg, int lkn, double *knots, point3d *cpoints,
                      double w, double penalty_param[MENGERC_NPPARAM],
                      int nqkn, int npthr, int opt, int maxit,
                      void (*outiter)(void *usrdata,
                                      boolean ppopt, int mdi,
                                      int itres, int it, double f, double g),
                      void *usrdata )
{
  mengerc_data md;
  int          i;
  boolean      finished;

  if ( !mengerc_InitMCOptimization ( deg, lkn, knots, cpoints, w,
                                     penalty_param, nqkn, npthr, opt, &md ) )
    return false;
  if ( outiter )
    outiter ( usrdata, md.ppopt, md.mdi, 0, -2, md.lastf, md.gn );
  finished = false;
  for ( i = 0;  i < maxit;  i++ ) {
    if ( !mengerc_IterMCOptimization ( &md, &finished ) )
      goto failure;
    if ( finished )
      break;
    if ( outiter )
      outiter ( usrdata, md.ppopt, md.mdi, i+1, md.itres, md.lastf, md.gn );
  }
  mengerc_UntieTheCurve ( &md );
  return true;

failure:
  mengerc_UntieTheCurve ( &md );
  return false;
} /*mengerc_OptimizeMengerCurvature*/

