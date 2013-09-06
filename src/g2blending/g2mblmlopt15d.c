
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
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
#include "g2mblmlprivated.h"

#define _DEBUG
/*#define __DEBUG*/

/* ///////////////////////////////////////////////////////////////////////// */
static double _g2mbl_MLAuxNuFuncd ( int nkn, double *qcoeff, double **Nitabs,
               double **Jac, int nv, point3d *mvcp, int nvcp, int *vncpi,
               int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
               double *ftab, int hsize, int *hprof, double **hrows, double **lhrows,
               double *grad, double *dcoeff, point3d *auxmvcp, double nu )
{
  int    nvars;
  double f;

  nvars = 3*nvcp;
#ifdef _DEBUG
printf ( "F" );
#endif
  if ( _g2bl_ShiftDecompHessiand ( nvars, hsize, hprof,
                                   hrows[0], lhrows[0], lhrows, nu ) ) {
    pkn_NRBLowerTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, dcoeff, 1, dcoeff );
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
    f = g2mbl_MLFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                        ndomel, domelind, domelem, domelcpind, true, ftab );
#ifdef _DEBUG
printf ( " nu = %10g, f = %10g\n", nu, f );
#endif
    return f;
  }
  else {
#ifdef _DEBUG
printf ( " nu = %10g, f = %10g\n", nu, MYINFINITY );
#endif
    return MYINFINITY;
  }
} /*_g2mbl_MLAuxNuFuncd*/

boolean g2mbl_MLOptBlockBd ( void *data, int bl )
{
  void            *sp;
  mesh_ml_optdata *d;
  mlblock_desc    *bd, *bd1;
  point3d         *mvcp, *auxmvcp;
  int             nv, nvcp, nvars, ndomel, *vncpi, *domelind, *domelcpind;
  meshdom_elem    *domelem;
  double          *ftab1, *ftab2, *gtab1, *gtab2, *htab, *coeff, *dcoeff, *grad;
  double          **Nitabs, **Jac, *qcoeff;
  int             *hprof, hsize;
  double          **hrows, **lhrows;
  int             nkn;
  int             i, j, k, kk;
  double          lco, ldco, func, gn, gna, eps, f0, gn0, thr;
  double          lmin, lmax, ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  boolean         recalc, fg_ok, advance;

#define MLFUNC(nu) \
_g2mbl_MLAuxNuFuncd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, nvcp, vncpi, \
        ndomel, domelind, domelem, domelcpind, ftab2, \
        hsize, hprof, hrows, lhrows, grad, dcoeff, auxmvcp, nu )

#define RECORD_MIN(g,fnu) \
  { if ( fnu < fge ) { \
      memcpy ( coeff, dcoeff, nvars*sizeof(double) ); \
      ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) ); \
      fge = fnu; \
      d->nu[1] = ge = g; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
        /* copy data to local variables */
  d = (mesh_ml_optdata*)data;
  nv         = d->nv;
  mvcp       = d->mvcp;
  domelem    = d->domelem;
  domelcpind = d->domelcpind;
  ftab1      = d->ftab1;
  gtab1      = d->gtab1;
  htab       = d->htab1;
  if ( d->currentblock ) {
    nkn    = d->nkn1;
    Nitabs = d->aNitabs;
    Jac    = d->aJac;
    qcoeff = d->aqcoeff;
    ftab2  = d->ftab1;
    gtab2  = d->gtab1;
  }
  else {
    nkn    = d->nkn2;
    Nitabs = d->bNitabs;
    Jac    = d->bJac;
    qcoeff = d->bqcoeff;
    ftab2  = d->ftab2;
    gtab2  = d->gtab2;
  }
  eps = bl == 0 ? EPS2 : EPS1;

  bd   = &d->bd[bl];
  nvcp     = bd->nvcp;
  nvars    = bd->nvars;
  ndomel   = bd->ndomel;
  vncpi    = bd->vncpi;
  domelind = bd->domelind;
  hprof    = bd->hprof;
  hsize    = bd->hsize;
  hrows    = bd->hrows;
  lhrows   = bd->lhrows;
  bd->fghflag &= ~FLAG_ADVANCE;
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvars );
  if ( !coeff ) {
printf ( "%s\n", ERRMSG_2 );
    goto failure;
  }
  dcoeff = &coeff[nvars];
  grad = &dcoeff[nvars];
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp ) {
printf ( "%s\n", ERRMSG_2 );
    goto failure;
  }

  ga = fgc = 0.0;
  if ( d->log_level > 1 ) {
    f0 = gn0 = gna = MYINFINITY;
  }
/*  d->nu[1] = -1.0; */
  fg_ok = advance = false;
  for ( k = kk = 0; ; k++ ) {
    memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
    lco = 0.0;
    for ( i = 0; i < nvcp; i++ )
      lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
    lco = sqrt ( lco );

        /* compute the functional value, gradient and Hessian */
    if ( !(bd->fghflag & FLAG_H) ) {
      if ( bl > 0 && !(bd->fghflag & FLAG_ADVANCE) ) {
        i = bl;
        do {
          i = (i-1)/2;
          bd1 = &d->bd[i];
        } while ( i > 0 && !bd1->iHbl );
        if ( !bd1->iHbl || !(bd1->fghflag & FLAG_H) )
          goto comp_Hessian;
        if ( !g2mbl_MLGetHessianRowsd ( nv, bd1->nvcp, bd1->vncpi,
                        bd1->nHbl, bd1->iHbl, bd1->cHbl, bd1->tHbl, d->Hbl,
                        nvcp, vncpi, hsize, hprof, hrows ) ) {
printf ( "%s\n", ERRMSG_25 );
          goto failure;
        }
        bd->fghflag &= ~FLAG_LH;
      }
      else {
comp_Hessian:
#ifdef _DEBUG
printf ( " H" );
#endif
#ifdef G2MBL_TIME_IT
        pkv_Tic ( NULL );
#endif
        if ( !g2mbl_MLFuncGradHessianBd ( d->nkn1, d->aqcoeff, d->aNitabs,
                     d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp,
                     ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                     true, ftab1, gtab1, htab,
                     &fga, dcoeff, hsize, hprof, hrows ) ) {
printf ( "%s\n", ERRMSG_27 );
          goto failure;
        }
#ifdef G2MBL_TIME_IT
        d->time_h += pkv_Toc ( NULL );
#endif
        bd->fghflag |= FLAG_H | FLAG_CH;
        bd->fghflag &= ~FLAG_LH;
        if ( nkn == d->nkn1 ) {
          func = fga;
          memcpy ( grad, dcoeff, nvars*sizeof(double) );
          bd->fghflag |= FLAG_F | FLAG_G;
        }
        else
          bd->fghflag &= ~(FLAG_F | FLAG_G);
      }
    }
    if ( !(bd->fghflag & FLAG_G) ) {
#ifdef _DEBUG
printf ( "G" );
#endif
      if ( bl != d->currentblock && k == 0 )
        recalc = !(d->bd[(bl-1)/2].fghflag & FLAG_G);  
      else
        recalc = true;
      if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                 ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                 recalc, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_28 );
        goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
    }
    else if ( k == 0 ) {
      if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                   ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                   false, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_28 );
          goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
    }
    else {  /* nothing to do, it is already computed */
    }
    fg_ok = true;

    gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
    if ( k == 0 ) {
      f0 = func;
      gn0 = gn;
    }
        /* try to decompose the Hessian */
    if ( !(bd->fghflag & FLAG_LH) ) {
      memcpy ( lhrows[0], hrows[0], hsize*sizeof(double) );
      if ( pkn_NRBSymCholeskyDecompd ( nvars, hprof, lhrows[0], lhrows, NULL ) )
        bd->fghflag |= FLAG_LH;
      else {
#ifdef _DEBUG
printf ( "!" );
#endif
        fga = MYINFINITY;
        goto lm_trajectory;
      }
    }
        /* try the Newton method */
#ifdef _DEBUG
printf ( "+" );
#endif
    pkn_NRBLowerTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, dcoeff, 1, dcoeff );
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    if ( ldco > BETA*lco ) {
      fga = MYINFINITY;
      goto lm_trajectory;
    }
    else if ( ldco <= EPS2*lco ) {
      advance = true;
      fg_ok = false;
      goto way_out;
    }

    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
    ga = 0.0;
    if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
               ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
               true, ftab2, gtab2, &fga, dcoeff ) ) {
printf ( "%s\n", ERRMSG_28 );
      goto failure;
    }
    bd->fghflag &= ~(FLAG_F | FLAG_G);  /* until acceptance */
    gna = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    if ( ldco >= eps*lco && (fga >= func || gna > 0.75*gn) )
      goto lm_trajectory;
    memcpy ( grad, dcoeff, nvars*sizeof(double) );
    fge = fga;
        /* more iterations of the Newton method */
    for ( j = 0; j < nvcp; j++ )
      mvcp[vncpi[j]] = auxmvcp[vncpi[j]];
    bd->fghflag |= FLAG_F | FLAG_G | FLAG_ADVANCE;  /* accepted */
    fg_ok = true;
    if ( ldco < eps*lco ) {
      advance = true;
      bd->fghflag |= FLAG_ADVANCE;
      goto way_out;
    }
    for ( i = 0, thr = 0.5; i < MAXNTN; i++, thr*= 0.5 ) {
      if ( fge >= func || gna > thr*gn ) {
        fge = func;
        gna = gn;
        fg_ok = true;
        goto next_iter;
      }
      for ( j = 0; j < nvcp; j++ )
        mvcp[vncpi[j]] = auxmvcp[vncpi[j]];
      fg_ok = true;
      bd->fghflag |= FLAG_F | FLAG_G;  /* accepted */
      if ( d->currentblock != 0 && gna <= 0.25*gn ) {
        advance = true;
        bd->fghflag |= FLAG_ADVANCE;
        goto way_out;
      }
      func = fge;
      gn = gna;
#ifdef _DEBUG
printf ( "," );
#endif
      pkn_NRBLowerTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, grad, 1, dcoeff );
      pkn_NRBUpperTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, dcoeff, 1, dcoeff );
      ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
      if ( ldco < eps*lco ) {
        advance = true;
        fg_ok = false;
        goto way_out;
      }
      _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
      if ( i < MAXNTN-1 ) {
#ifdef _DEBUG
printf ( "G" );
#endif
        if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                   ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                   true, ftab2, gtab2, &fge, grad ) ) {
printf ( "%s\n", ERRMSG_28 );
          goto failure;
        }
        bd->fghflag &= ~(FLAG_F | FLAG_G);  /* until acceptance */
        gna = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        fg_ok = false;
      }
    }
    goto next_iter;

lm_trajectory:
        /* minimization along the Levenberg-Marquardt trajectory */
    fg_ok = false;
    ga = 0.0;  fga = MYINFINITY;
    gb = ge = MYINFINITY;
    fgb = fge = func;
    memset ( coeff, 0, nvars*sizeof(double) );
          /* get the initial nu for bracketing */
    if ( d->nu[0] <= 0.0 ) {
      if ( !pkn_NRBSymFindEigenvalueIntervald ( nvars, hprof, hrows[0], hrows,
                                                &lmin, &lmax ) ) {
printf ( "%s\n", ERRMSG_31 );
        goto failure;
      }
      d->nu[0] = d->nu[1] = gc = 0.0001*(lmax-lmin);
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
      if ( _g2mbl_DivideIntervald ( &ga, &gc, &gd, &gb,
                                    &fga, &fgc, &fgd, &fgb ) ) {
        fgc = MLFUNC ( gc );
        RECORD_MIN ( gc, fgc )
      }
      else {
        fgd = MLFUNC ( gd );
        RECORD_MIN ( gd, fgd )
      }
      fm = min ( fga-fge, fgb-fge );
      if ( fm < GRTHR*(func-fge) ||
           fm < EPS2*fge )
        break;
    }
finish_it:
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)coeff, mvcp );
    bd->fghflag |= FLAG_ADVANCE;
    kk ++;
next_iter:
    bd->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H | FLAG_LH);

    if ( bl == d->currentblock )
      goto way_out;
  }

way_out:
  if ( d->log_level > 1 ) {
    if ( bl == d->currentblock ) {
      if ( !fg_ok ) {
#ifdef _DEBUG
        printf ( "G" );
#endif
        g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                            ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                            !fg_ok, ftab2, gtab2, &fge, grad );
        bd->fghflag |= FLAG_F | FLAG_G;
        gna = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        fg_ok = true;
      }
      printf ( "\n func %12g -> %12g, gn %12g -> %12g", f0, fge, gn0, gna );
    }
  }
  if ( bl == d->currentblock )
    _g2mbl_MLNextBlockNumd ( d, advance );

  d->lastblock = bl;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  d->lastblock = bl;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLOptBlockBd*/

