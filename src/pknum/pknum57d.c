
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
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

#define MYINFINITY 1.0e308
#define MAXNITER   10
#define MAXGITER   16
#define THR1        0.9
#define GRTHR       0.01
#define EPS2        1.0e-7

#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

/* ///////////////////////////////////////////////////////////////////////// */
static double _pkn_NLM_NuFuncd ( int n, void *usrdata,
                                 pkn_NLMTevalfuncd funcf, pkn_NLMTtunnelfuncd tunnel,
                                 double nu, double *x, double *grad,
                                 double *incr, double *auxx,
                                 double *hess, double *lhess,
                                 boolean *went_out, double *incn )
{
  int     i, k;
  double  fnu;

  *incn = MYINFINITY;
  memcpy ( lhess, hess, (n*(n+1)/2)*sizeof(double) );
  for ( i = k = 0;  i < n;  i++, k += i+1 )
    lhess[k] += nu;
  if ( pkn_CholeskyDecompd ( n, lhess ) ) {
    pkn_LowerTrMatrixSolved ( n, lhess, 1, 1, grad, 1, incr );
    pkn_UpperTrMatrixSolved ( n, lhess, 1, 1, incr, 1, incr );
    *incn = sqrt ( pkn_ScalarProductd ( n, incr, incr ) );
    pkn_SubtractMatrixd ( 1, n, 0, x, 0, incr, 0, auxx );
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
} /*_pkn_NLM_NuFuncd*/

static boolean _pkn_DivideIntervald ( double *ga, double *gc, double *gd, double *gb,
                                      double *fa, double *fc, double *fd, double *fb )
{
  double  x[4], y[4], a, b, c, x1, x2;
  int     i, j;
  boolean ok;
  
        /* construct the cubic polynomial of interpolation */
  x[0] = log(*ga);  x[1] = log(*gc);  x[2] = log(*gd);  x[3] = log(*gb);
  y[0] = *fa;       y[1] = *fc;       y[2] = *fd;       y[3] = *fb;
          /* compute the divided differences */
  for ( i = 1; i <= 3; i++ )
    for ( j = 3; j >= i; j-- )
      y[j] = (y[j]-y[j-1])/(x[j]-x[j-i]);
          /* find the derivative coefficients */
  a = 6.0*y[3];
  b = 2.0*y[2] - 4.0*y[3]*(x[0]+x[1]+x[2]);
  c = y[1] - y[2]*(x[0]+x[1]) + 2.0*y[3]*(x[0]*(x[1]+x[2])+x[1]*x[2]);
  ok = false;
  if ( a ) {
    if ( pkn_SolveSqEqd ( b/(2.0*a), c/a, &x1, &x2 ) ) {
      if ( a > 0.0 )
        x1 = x2;
      ok = true;
    }
  }
  else {
    x1 = -c/b;
    ok = true;
  }
  if ( ok ) {
    x2 = y[3];
    for ( i = 2; i >= 0; i-- )
      x2 = x2*(x1-x[i]) + y[i];
    if ( x2 <= 0.25*(*fa+*fb) )
      ok = false;
  }
  if ( *fc < *fd ) {
    *gb = *gd;  *fb = *fd;
    if ( ok && x1 > x[0] && x1 < x[2] ) {
      if ( x1 > x[1] ) {
        *gd = exp ( x1 );
        return false;
      }
      *gd = *gc;  *fd = *fc;
      if ( x1 < x[1] )
        *gc = exp ( x1 );
      else
        *gc = GRDIV ( *ga, *gd );
    }
    else {
      *gd = *gc;  *fd = *fc;
      *gc = GRDIV ( *ga, *gd );
    }
    return true;
  }
  else {
    *ga = *gc;  *fa = *fc;
    if ( ok && x1 > x[1] && x1 < x[3] ) {
      if ( x1 < x[2] ) {
        *gc = exp ( x1 );
        return true;
      }
      *gc = *gd;  *fc = *fd;
      if ( x1 > x[2] )
        *gd = exp ( x1 );
      else
        *gd = GRDIV ( *gb, *gc );
    }
    else {
      *gc = *gd;  *fc = *fd;
      *gd = GRDIV ( *gb, *gc );
    }
    return false;
  }
} /*_pkn_DivideIntervald*/

static boolean _pkn_ComputeDeltaQd ( int n, const double *hess,
                                     const double *grad, const double *dcoeff,
                                     double *dq )
{
  void   *sp;
  double *aux;

  sp = pkv_GetScratchMemTop ();
  aux = pkv_GetScratchMemd ( n );
  if ( !aux )
    goto failure;
  pkn_SymMatrixMultd ( n, hess, 1, 1, dcoeff, 1, aux );
  pkn_MatrixLinCombd ( 1, n, 0, aux, 0.5, 0, grad, -1.0, 0, aux );
  *dq = pkn_ScalarProductd ( n, aux, dcoeff );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_pkn_ComputeDeltaQd*/

/* ///////////////////////////////////////////////////////////////////////// */
int pkn_NLMIterd ( int n, void *usrdata, double *x,
                   pkn_NLMTevalfuncd funcf, pkn_NLMTevalfuncgd funcfg,
                   pkn_NLMTevalfuncghd funcfgh, pkn_NLMTtransfuncd trans,
                   pkn_NLMTtunnelfuncd tunnel,
                   double lowerbound, double eps, double delta,
                   double *nu )
{
  void    *sp;
  double  f, *grad, *incr, *minx, *auxx, *auxgr, *hess, *lhess;
  double  gn, gna, incn, incne, df, dq, rk, lmin, lmax;
  double  ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  int     i;
  boolean went_out;
  int     result;

#define LMFUNC(nu) \
  _pkn_NLM_NuFuncd ( n, usrdata, funcf, tunnel, \
                     nu, x, grad, incr, auxx, hess, lhess, &went_out, &incn )

/* this macrodefinition must always be preceded by the previous one */
#define RECORD_MIN(nu,fnu) \
  { if ( went_out && incn < delta ) { \
      result = PKN_LMT_FOUND_BARRIER; \
      goto finish_lmt; \
    } \
    if ( fnu < lowerbound ) { \
      result = PKN_LMT_ERROR; \
      goto way_out; \
    } \
    if ( fnu < fge ) { \
      memcpy ( minx, auxx, n*sizeof(double) ); \
      ge = nu;  fge = fnu;  incne = incn; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
  result = PKN_LMT_CONTINUE_N;
  grad = pkv_GetScratchMemd ( 5*n );
  hess = pkv_GetScratchMemd ( n*(n+1) );
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
        /* compute the functional, gradient and Hessian */
  if ( !funcfgh ( n, usrdata, x, &f, grad, hess ) ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  gn = sqrt ( pkn_ScalarProductd ( n, grad, grad ) );
        /* try to decompose the Hessian */
  memcpy ( lhess, hess, ((n*(n+1))/2)*sizeof(double) );
  if ( !pkn_CholeskyDecompd ( n, lhess ) ) {
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
  pkn_LowerTrMatrixSolved ( n, lhess, 1, 1, grad, 1, incr );
  pkn_UpperTrMatrixSolved ( n, lhess, 1, 1, incr, 1, incr );
  pkn_SubtractMatrixd ( 1, n, 0, x, 0, incr, 0, auxx );
  incn = sqrt ( pkn_ScalarProductd ( n, incr, incr ) );
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
  memcpy ( x, auxx, n*sizeof(double) );
  df = fga - f;
  if ( !_pkn_ComputeDeltaQd ( n, hess, grad, incr, &dq ) ) {
    result = PKN_LMT_ERROR;
    goto way_out;
  }
  rk = df/dq;
  f = fga;
  ge = -1.0;
  gna = sqrt ( pkn_ScalarProductd ( n, auxgr, auxgr ) );
  if ( gna < eps ) {
    result = PKN_LMT_FOUND_MINIMUM;
    goto way_out;
  }
  else if ( gna > gn )
    goto lm_trajectory;
  else if ( rk < THR1 /*|| gna > 0.5*gn*/ )
    goto way_out;
        /* additional Newton method steps */
  for ( i = 0; i < MAXNITER; i++ ) {
    gn = gna;
    pkn_LowerTrMatrixSolved ( n, lhess, 1, 1, auxgr, 1, incr );
    pkn_UpperTrMatrixSolved ( n, lhess, 1, 1, incr, 1, incr );
    pkn_SubtractMatrixd ( 1, n, 0, x, 0, incr, 0, auxx );
    incn = sqrt ( pkn_ScalarProductd ( n, incr, incr ) );
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
    memcpy ( x, auxx, n*sizeof(double) );
    f = fga;
    gna = sqrt ( pkn_ScalarProductd ( n, auxgr, auxgr ) );
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
  result = PKN_LMT_CONTINUE_LM;
  ga = 0.0;              fga = MYINFINITY;
  gb = ge = MYINFINITY;  fgb = fge = f;
  memcpy ( minx, x, n*sizeof(double) );
  incne = MYINFINITY;
        /* get the initial nu for bracketing */
  if ( *nu <= 0.0 ) {
    pkn_SymMatFindEigenvalueIntervald ( n, hess, &lmin, &lmax );
    gc = 0.001*(lmax-lmin);
  }
  else
    gc = *nu;
  fgc = LMFUNC ( gc );
  RECORD_MIN ( gc, fgc );
        /* bracketing */
  for ( i = 0; fgc >= fga; i++ ) {
    ga = gc;  fga = fgc;
    gc *= (double)(i+2);  /* trying the Fibonacci sequence */
    fgc = LMFUNC ( gc );
    RECORD_MIN ( gc, fgc );
  }
  gd = gc*(double)(i+2);
  fgd = LMFUNC ( gd );
  RECORD_MIN ( gd, fgd );
  if ( fgd < fgc ) {
    do {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = gc*(double)(i+2);  i++;
      fgd = LMFUNC ( gd );
      RECORD_MIN ( gd, fgd );
    } while ( fgd < fgc );
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
      if ( i >= MAXGITER && fge < f ) {
        *nu = ge;
        goto finish_lmt;
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
    if ( _pkn_DivideIntervald ( &ga, &gc, &gd, &gb,
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
  else if ( incne < delta )
    result = PKN_LMT_FOUND_ZEROGRAD;
finish_lmt:
  memcpy ( x, minx, n*sizeof(double) );

way_out:
  pkv_SetScratchMemTop ( sp );
  return result;
} /*pkn_NLMIterd*/

