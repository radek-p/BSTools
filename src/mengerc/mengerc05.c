
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

#define MAXLKNOT 1000
#define MAXITER  200

/* ////////////////////////////////////////////////////////////////////////// */
boolean mengerc_OptimizeMengerCurvature (
                      int deg, int lkn, double *knots, point3d *cpoints,
                      double w, double penalty_param[5],
                      int nqkn, int nthr, int opt,
                      void (*outiter)(void *usrdata,
                                      int itres, int it, double f, double g),
                      void *usrdata )
{
  void         *sp;
  mengerc_data md;
  int          n, i, fcnt;
  double       *x, *g, f, gn;
  double       nu, lastf, lastgn;
  boolean      opt_res;
  int          itres;

  sp = pkv_GetScratchMemTop ();
  if ( !mengerc_BindACurve ( &md, deg, lkn, knots, cpoints, nqkn,
                             w, penalty_param, w >= 4.0 ) )
    goto failure;

  n = 3*(lkn-2*deg);
  x = pkv_GetScratchMemd ( n );
  g = pkv_GetScratchMemd ( n );
  if ( !x || !g )
    goto failure;

  md.npthr = nthr;
  if ( !mengerc_IntegralMengerfg ( n, (void*)&md, &cpoints[0].x, &f, g ) )
    goto failure;

  memcpy ( x, cpoints, n*sizeof(double) );
  md.pretransf = true;
  if ( !mengerc_IntegralMengerTransC ( n, (void*)&md, x ) )
    goto failure;

  switch ( opt ) {
case MENGERC_OPT_FULL1:
    mengerc_OptPenaltyParams1 ( &md, false/*true*/ );
    break;
case MENGERC_OPT_FULL2:
    mengerc_OptPenaltyParams2 ( &md );
    break;
case MENGERC_OPT_PART:
    mengerc_OptPenaltyParams3 ( &md );
    break;
case MENGERC_OPT_NONE:
default:
    break;
  }
  mengerc_IntegralMengerfg ( n, (void*)&md, x, &f, g );
  lastf = f;
  gn = sqrt ( pkn_ScalarProductd ( n, g, g ) );
    if ( outiter )
      outiter ( usrdata, 0, -1, f, gn );

  nu = -1.0;
  md.pretransf = true;
  fcnt = 0;
  for ( i = 0;  i < MAXITER;  i++ ) {
    lastf = f;
    lastgn = gn;
    itres = pkn_NLMIterd ( n, (void*)&md, x,
                           mengerc_IntegralMengerf,
                           mengerc_IntegralMengerfg,
                           mengerc_IntegralMengerfgh,
                           mengerc_IntegralMengerTransC,
                           mengerc_HomotopyTest, 0.0, 1.0e-4, 1.0e-7, &nu );

    mengerc_IntegralMengerfg ( n, (void*)&md, x, &f, g );
    gn = sqrt ( pkn_ScalarProductd ( n, g, g ) );
    if ( outiter )
      outiter ( usrdata, i+1, itres, f, gn );

    switch ( itres ) {
  case PKN_LMT_ERROR:
  case PKN_LMT_CROSSED_LIMIT:
      goto failure;

  case PKN_LMT_CONTINUE_N:
      if ( fcnt > 0 )
        fcnt --;
      break;

  case PKN_LMT_CONTINUE_LM:
  case PKN_LMT_CONTINUE_LM_P:
      fcnt ++;
      break;

  case PKN_LMT_FOUND_MINIMUM:
  case PKN_LMT_FOUND_ZEROGRAD:
  case PKN_LMT_FOUND_ZEROGRAD_P:
  case PKN_LMT_FOUND_BARRIER:
      goto finish;

  case PKN_LMT_NO_PROGRESS:
      if ( opt != MENGERC_OPT_NONE && fcnt > 0 )
        goto opt_param;
      else
        goto finish;

  default:
      goto failure;
    }

    md.mdi = mengerc_ModifyRemotestPoint ( lkn-deg, cpoints, &md.sc, md.mdi );
    if ( fcnt >= 5+(int)(sqrt(i)) ) {
opt_param:
      if ( f <= 0.9*lastf || gn <= 0.9*lastgn )
        goto next_iter1;
      opt_res = false;
      switch ( opt ) {
    case MENGERC_OPT_FULL1:
        opt_res = mengerc_OptPenaltyParams1 ( &md, false );
        break;
    case MENGERC_OPT_FULL2:
        opt_res |= mengerc_OptPenaltyParams2 ( &md );
        break;
    case MENGERC_OPT_PART:
        opt_res |= mengerc_OptPenaltyParams3 ( &md );
        break;
    case MENGERC_OPT_NONE:
    default:
        break;
      }
      if ( opt_res )
        fcnt = 0;
      else
        fcnt /= 2;
    }
next_iter1:
    ;
  }

finish:
  memcpy ( penalty_param, md.penalty_param, 5*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_OptimizeMengerCurvature*/

