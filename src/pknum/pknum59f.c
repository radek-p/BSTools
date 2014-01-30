
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

#define MYINFINITY 1.0e38
#define MAXGITER   16
#define GRTHR       0.01
#define EPS2        1.0e-7

#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

/* ///////////////////////////////////////////////////////////////////////// */
static float _pkn_SD_NuFuncf ( int n, void *usrdata,
                               pkn_NLMTevalfuncf funcf, pkn_NLMTtunnelfuncf tunnel,
                               float nu, float *x, float *grad,
                               float *incr, float *auxx,
                               boolean *went_out, float *incn )
{
  float fnu; 

  pkn_MultMatrixNumf ( 1, n, 0, grad, 1.0/nu, 0, incr );
  *incn = sqrt ( pkn_ScalarProductf ( n, incr, incr ) );
  pkn_SubtractMatrixf ( 1, n, 0, x, 0, incr, 0, auxx ); 
  if ( tunnel ) {
    if ( tunnel ( n, usrdata, x, auxx, went_out ) ) {
      if ( *went_out )
        return MYINFINITY;
    }
    else
      return -MYINFINITY;  /* error */
  }
  else
    *went_out = false;
  if ( funcf ( n, usrdata, auxx, &fnu ) )
    return fnu;
  else
    return -MYINFINITY;  /* error */
} /*_pkn_SD_NuFuncf*/

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_SDIterf ( int n, void *usrdata, float *x,
                  pkn_NLMTevalfuncf funcf, pkn_NLMTevalfuncgf funcfg,
                  pkn_NLMTtransfuncf trans, pkn_NLMTtunnelfuncf tunnel,
                  float lowerbound, float eps, float delta,
                  float *nu )
{
  void    *sp;
  float   f, *grad, *incr, *minx, *auxx;
  float   incn, incne;
  float   ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  int     i, j;
  boolean progress, went_out;
  int     result;

#define SDFUNC(nu) \
  _pkn_SD_NuFuncf ( n, usrdata, funcf, tunnel, \
                    nu, x, grad, incr, auxx, &went_out, &incn )

/* this macrodefinition must always be preceded by the previous one */
#define RECORD_MIN(nu,fnu) \
  { if ( went_out && incn < delta ) { \
      result = PKN_SD_FOUND_BARRIER; \
      goto finish_sd; \
    } \
    if ( fnu < lowerbound ) { \
      result = PKN_SD_CROSSED_LIMIT; \
      goto way_out; \
    } \
    if ( fnu < fge ) { \
      memcpy ( minx, auxx, n*sizeof(float) ); \
      ge = nu;  fge = fnu;  incne = incn; \
      progress = true; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
  result = PKN_SD_CONTINUE;
  grad = pkv_GetScratchMemf ( 4*n );
  if ( !grad ) {
    result = PKN_SD_ERROR;
    goto way_out;
  }
  incr = &grad[n];
  minx = &incr[n];
  auxx = &minx[n];

        /* pre-transformation */
  if ( trans )
    if ( !trans ( n, usrdata, x ) ) {
      result = PKN_SD_ERROR;
      goto way_out;
    }
        /* compute the function and gradient */
  progress = false;
  if ( !funcfg ( n, usrdata, x, &f, grad ) ) {
    result = PKN_SD_ERROR;
    goto way_out;
  }
        /* minimization along the steepest descent direction */
  ga = 0.0;              fga = MYINFINITY;
  gb = ge = MYINFINITY;  fgb = fge = f;
  memcpy ( minx, x, n*sizeof(float) );
  incne = MYINFINITY;
        /* get the initial nu for bracketing */
  if ( *nu <= 0.0 )
    gc = 1.0;
  else
    gc = *nu;
  fgc = SDFUNC ( gc );
  RECORD_MIN ( gc, fgc );
        /* bracketing */
  for ( i = 0; fgc >= fga; i++ ) {
    ga = gc;  fga = fgc;
    gc *= (float)(i+2);  /* trying the Fibonacci sequence */
    fgc = SDFUNC ( gc );
    RECORD_MIN ( gc, fgc );
    if ( i >= MAXGITER ) {
      if ( progress ) {
        *nu = ge;
        goto finish_sd;
      }
      else {
        result = PKN_SD_NO_PROGRESS;
        goto way_out;
      }
    }
  }
  gd = gc*(float)(i+2);
  fgd = SDFUNC ( gd );
  RECORD_MIN ( gd, fgd );
  if ( fgd < fgc ) {
    for ( j = 0; j < MAXGITER; j++ ) {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = gc*(float)(i+2);  i++;
      fgd = SDFUNC ( gd );
      RECORD_MIN ( gd, fgd );
      if ( fgd < fgc )
        goto cont1;
    }
    if ( progress ) {
      *nu = ge;
      goto finish_sd;
    }
    else {
      result = PKN_SD_NO_PROGRESS;
      goto way_out;
    }
cont1:
    gb = gd;  fgb = fgd;
  }
  else {
    gb = gd;  fgb = fgd;
    for ( i = 0; ; i++ ) {
      if ( ga > 0.0 )        gd = sqrt ( ga*gc );
      else if ( gc > 1.0+4 ) gd = sqrt ( gc );
      else                   gd = 0.01*gc;
      fgd = SDFUNC ( gd );
      RECORD_MIN ( gd, fgd );
      if ( i >= MAXGITER ) {
        if ( fge < f ) {
          *nu = ge;
          goto finish_sd;
        }
        else {
          result = PKN_LMT_NO_PROGRESS;
          goto way_out;
        }
      }
      if ( fgd >= fga )     { ga = gd;  fga = fgd; }
      else if ( fgd < fgc ) { gb = gc;  fgb = fgc;  gc = gd;  fgc = fgd; }
      else                  { ga = gd;  fga = fgd;  break; }
    }
  }
        /* find the second point in the bracket */
  if ( fga < fgb ) {
    gd = gc;  fgd = fgc;
    gc = GRDIV ( ga, gd );
    fgc = SDFUNC ( gc );
    RECORD_MIN ( gc, fgc );
  }
  else {
    gd = GRDIV ( gb, gc );
    fgd = SDFUNC ( gd );
    RECORD_MIN ( gd, fgd );
  }
        /* subdivision */
  for ( i = 0; i < MAXGITER;  i++ ) {
    if ( _pkn_DivideIntervalf ( &ga, &gc, &gd, &gb,
                                &fga, &fgc, &fgd, &fgb ) ) {
      fgc = SDFUNC ( gc );
      RECORD_MIN ( gc, fgc );
    }
    else {
      fgd = SDFUNC ( gd );
      RECORD_MIN ( gd, fgd );
    }
    fm = min ( fga-fge, fgb-fge );
    if ( fm < GRTHR*(f-fge) || fm < EPS2*fge )
      break;
  }
  if ( ge < MYINFINITY )
    *nu = ge;
  if ( fge >= f ) {
    result = PKN_SD_NO_PROGRESS;
    goto way_out;
  }
  else if ( !progress && incne < delta )
    result = PKN_SD_FOUND_ZEROGRAD;
finish_sd:
  memcpy ( x, minx, n*sizeof(float) );

way_out:
  pkv_SetScratchMemTop ( sp );
  return result;
} /*pkn_SDIterf*/

