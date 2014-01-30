
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

#define DEBUG

#define MYINFINITY 1.0e38
#define MAXNITER   10
#define MAXGITER   16
#define THR1        0.9
#define GRTHR       0.01
#define EPS2        1.0e-5

#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

/* ///////////////////////////////////////////////////////////////////////// */
static float _pkn_NLM_NuFuncf ( int n, void *usrdata,
                                pkn_NLMTevalfuncf funcf, pkn_NLMTtunnelfuncf tunnel,
                                float nu, float *x, float *grad,
                                float *incr, float *auxx,
                                float *hess, float *lhess,
                                boolean *went_out, float *incn )
{
  int     i, k;
  float   fnu;

  *incn = MYINFINITY;
  memcpy ( lhess, hess, (n*(n+1)/2)*sizeof(float) );
  for ( i = k = 0;  i < n;  i++, k += i+1 )
    lhess[k] += nu;
  if ( pkn_CholeskyDecompf ( n, lhess ) ) {
    pkn_LowerTrMatrixSolvef ( n, lhess, 1, 1, grad, 1, incr );
    pkn_UpperTrMatrixSolvef ( n, lhess, 1, 1, incr, 1, incr );
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
  }
  else {
    *went_out = false;
    return MYINFINITY;
  }
} /*_pkn_NLM_NuFuncf*/

static boolean _pkn_ComputeDeltaQf ( int n, const float *hess,
                                     const float *grad, const float *dcoeff,
                                     float *dq )
{
  void  *sp;
  float *aux;

  sp = pkv_GetScratchMemTop ();
  aux = pkv_GetScratchMemf ( n );
  if ( !aux )
    goto failure;
  pkn_SymMatrixMultf ( n, hess, 1, 1, dcoeff, 1, aux );
  pkn_MatrixLinCombf ( 1, n, 0, aux, 0.5, 0, grad, -1.0, 0, aux );
  *dq = pkn_ScalarProductf ( n, aux, dcoeff );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_pkn_ComputeDeltaQf*/

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_NLMIterf ( int n, void *usrdata, float *x,
                   pkn_NLMTevalfuncf funcf, pkn_NLMTevalfuncgf funcfg,
                   pkn_NLMTevalfuncghf funcfgh, pkn_NLMTtransfuncf trans,
                   pkn_NLMTtunnelfuncf tunnel,
                   float lowerbound, float eps, float delta,
                   float *nu )
{
  void    *sp;
  float   f, *grad, *incr, *minx, *auxx, *auxgr, *hess, *lhess;
  float   gn, gna, incn, incne, df, dq, rk, lmin, lmax;
  float   ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  int     i, j;
  boolean positive, progress, went_out;
  int     result;

#define LMFUNC(nu) \
  _pkn_NLM_NuFuncf ( n, usrdata, funcf, tunnel, \
                     nu, x, grad, incr, auxx, hess, lhess, &went_out, &incn )

/* this macrodefinition must always be preceded by the previous one */
#define RECORD_MIN(nu,fnu) \
  { if ( went_out && incn < delta ) { \
      result = PKN_LMT_FOUND_BARRIER; \
      goto finish_lmt; \
    } \
    if ( fnu < lowerbound ) { \
      result = PKN_LMT_CROSSED_LIMIT; \
      goto way_out; \
    } \
    if ( fnu < fge ) { \
      memcpy ( minx, auxx, n*sizeof(float) ); \
      ge = nu;  fge = fnu;  incne = incn; \
      progress = true; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
  result = PKN_LMT_CONTINUE_N;
  grad = pkv_GetScratchMemf ( 5*n );
  hess = pkv_GetScratchMemf ( n*(n+1) );
  if ( !grad || !hess ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  incr = &grad[n];
  minx = &incr[n];
  auxx = &minx[n];
  auxgr = &auxx[n];
  lhess = &hess[(n*(n+1))/2];

        /* pre-transformation */
  if ( trans )
    if ( !trans ( n, usrdata, x ) ) {
      result = PKN_LMT_ERROR;
      goto way_out;
    }
        /* compute the function, gradient and Hessian */
  progress = false;
  if ( !funcfgh ( n, usrdata, x, &f, grad, hess ) ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  gn = sqrt ( pkn_ScalarProductf ( n, grad, grad ) );
        /* try to decompose the Hessian */
  memcpy ( lhess, hess, ((n*(n+1))/2)*sizeof(float) );
  positive = pkn_CholeskyDecompf ( n, lhess );
  if ( !positive ) {
    if ( gn < eps ) {
      result = PKN_LMT_FOUND_ZEROGRAD;
      goto way_out;
    }
    else
      goto lm_trajectory;
  }
  if ( gn < eps ) {
    result = PKN_LMT_FOUND_MINIMUM;
    goto way_out;
  }
        /* try the Newton method step */
  pkn_LowerTrMatrixSolvef ( n, lhess, 1, 1, grad, 1, incr );
  pkn_UpperTrMatrixSolvef ( n, lhess, 1, 1, incr, 1, incr );
  pkn_SubtractMatrixf ( 1, n, 0, x, 0, incr, 0, auxx );
  incn = sqrt ( pkn_ScalarProductf ( n, incr, incr ) );
  if ( tunnel ) {
    if ( !tunnel ( n, usrdata, x, auxx, &went_out ) ) {
      result = PKN_LMT_ERROR;
      goto way_out;
    }
    if ( went_out ) {
      fga = MYINFINITY;
      goto lm_trajectory;
    }
  }
  if ( incn < delta ) {
    result = PKN_LMT_FOUND_MINIMUM;
    goto way_out;
  }
  if ( !funcfg ( n, usrdata, auxx, &fga, auxgr ) ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  if ( fga >= f )
    goto lm_trajectory;
  memcpy ( x, auxx, n*sizeof(float) );
  df = fga - f;
  if ( !_pkn_ComputeDeltaQf ( n, hess, grad, incr, &dq ) ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  rk = df/dq;
  f = fga;
  ge = -1.0;
  gna = sqrt ( pkn_ScalarProductf ( n, auxgr, auxgr ) );
  if ( gna < eps ) {
    result = PKN_LMT_FOUND_MINIMUM;
    goto way_out;
  }
  else if ( gna > gn ) {
    progress = true;
    goto lm_trajectory;
  }
  else if ( rk < THR1 /*|| gna > 0.5*gn*/ )
    goto way_out;
        /* additional Newton method steps */
  for ( i = 0; i < MAXNITER; i++ ) {
    gn = gna;
    pkn_LowerTrMatrixSolvef ( n, lhess, 1, 1, auxgr, 1, incr );
    pkn_UpperTrMatrixSolvef ( n, lhess, 1, 1, incr, 1, incr );
    pkn_SubtractMatrixf ( 1, n, 0, x, 0, incr, 0, auxx );
    incn = sqrt ( pkn_ScalarProductf ( n, incr, incr ) );
    if ( tunnel ) {
      if ( !tunnel ( n, usrdata, x, auxx, &went_out ) ) {
        result = PKN_LMT_ERROR;
        goto way_out;
      }
      if ( went_out )
        goto way_out;
    }
    if ( incn < delta ) {
      result = PKN_LMT_FOUND_MINIMUM;
      goto way_out;
    }
    if ( !funcfg ( n, usrdata, auxx, &fga, auxgr ) ) {
      result = PKN_LMT_ERROR;
      goto way_out;
    }
    if ( fga >= f )
      goto way_out;
    memcpy ( x, auxx, n*sizeof(float) );
    f = fga;
    gna = sqrt ( pkn_ScalarProductf ( n, auxgr, auxgr ) );
    if ( gna < eps ) {
      result = PKN_LMT_FOUND_MINIMUM;
      goto way_out;
    }
    else if ( gna > 0.5*gn )
      goto way_out;
  }
  goto way_out;

lm_trajectory:
        /* minimization along the Levenberg-Marquardt trajectory */
  result = positive ? PKN_LMT_CONTINUE_LM_P : PKN_LMT_CONTINUE_LM;
  ga = 0.0;              fga = MYINFINITY;
  gb = ge = MYINFINITY;  fgb = fge = f;
  memcpy ( minx, x, n*sizeof(float) );
  incne = MYINFINITY;
        /* get the initial nu for bracketing */
  if ( *nu <= 0.0 ) {
    pkn_SymMatFindEigenvalueIntervalf ( n, hess, &lmin, &lmax );
    gc = 0.001*(lmax-lmin);
  }
  else
    gc = *nu;
  fgc = LMFUNC ( gc );
  RECORD_MIN ( gc, fgc );
        /* bracketing */
  for ( i = 0; fgc >= fga; i++ ) {
    ga = gc;  fga = fgc;
    gc *= (float)(i+2);  /* trying the Fibonacci sequence */
    fgc = LMFUNC ( gc );
    RECORD_MIN ( gc, fgc );
    if ( i >= MAXGITER ) {
      if ( progress ) {
        *nu = ge;
        goto finish_lmt;
      }
      else {
        result = PKN_LMT_NO_PROGRESS;
        goto way_out;
      }
    }
  }
  gd = gc*(float)(i+2);
  fgd = LMFUNC ( gd );
  RECORD_MIN ( gd, fgd );
  if ( fgd < fgc ) {
    for ( j = 0; j < MAXGITER; j++ ) {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = gc*(float)(i+2);  i++;
      fgd = LMFUNC ( gd );
      RECORD_MIN ( gd, fgd );
      if ( fgd < fgc )
        goto cont1;
    }
    if ( progress ) {
      *nu = ge;
      goto finish_lmt;
    }
    else {
      result = PKN_LMT_NO_PROGRESS;
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
      fgd = LMFUNC ( gd );
      RECORD_MIN ( gd, fgd );
      if ( i >= MAXGITER ) {
        if ( fge < f ) {
          *nu = ge;
          goto finish_lmt;
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
    fgc = LMFUNC ( gc );
    RECORD_MIN ( gc, fgc );
  }
  else {
    gd = GRDIV ( gb, gc );
    fgd = LMFUNC ( gd );
    RECORD_MIN ( gd, fgd );
  }
        /* subdivision */
  for ( i = 0; i < MAXGITER;  i++ ) {
    if ( _pkn_DivideIntervalf ( &ga, &gc, &gd, &gb,
                                &fga, &fgc, &fgd, &fgb ) ) {
      fgc = LMFUNC ( gc );
      RECORD_MIN ( gc, fgc );
    }
    else {
      fgd = LMFUNC ( gd );
      RECORD_MIN ( gd, fgd );
    }
    fm = min ( fga-fge, fgb-fge );
    if ( fm < GRTHR*(f-fge) || fm < EPS2*fge )
      break;
  }
  if ( ge < MYINFINITY )
    *nu = ge;
  if ( fge >= f ) {
    result = PKN_LMT_NO_PROGRESS;
    goto way_out;
  }
  else if ( !progress && incne < delta )
    result = positive ? PKN_LMT_FOUND_ZEROGRAD_P : PKN_LMT_FOUND_ZEROGRAD;
finish_lmt:
  memcpy ( x, minx, n*sizeof(float) );

way_out:
  pkv_SetScratchMemTop ( sp );
  return result;
} /*pkn_NLMIterf*/

