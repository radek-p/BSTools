
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/times.h>
#include <unistd.h>

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
#include "msgpool.h"

#define _DEBUG
/*#define __DEBUG*/

#define USE_PCG

/* ///////////////////////////////////////////////////////////////////////// */
static double _g2mbl_MLCAuxNuFuncd ( mesh_ml_optdata *d, int bl,
               int nkn, double *qcoeff, double **Nitabs,
               double **Jac, int nv, point3d *mvcp, int nvcp, int *vncpi,
               int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
               double *ftab,
               int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl,
               double *grad, double *dcoeff, point3d *auxmvcp, double nu )
{
  int             nvars, itm;
  double          f, gn2;
  mesh_ml_cg_data cgdata;
  boolean         result;

  nvars = 3*nvcp;
  _g2mbl_MLInvalSmallBlocks ( d, 2*bl+1 );
  _g2mbl_MLInvalSmallBlocks ( d, 2*bl+2 );
  d->bd[bl].fghflag &= ~(FLAG_CMH | FLAG_CMLH);
  if ( !_g2mbl_MLDecomposeBlockPrecond ( d, bl, nu, &result ) )
    return -1.0;
  if ( !result ) {
#ifdef _DEBUG
printf ( "\n" );
#endif
    return MYINFINITY;
  }
  cgdata.d = d;
  cgdata.bl = cgdata.cbl = bl;
  cgdata.nu = nu;
  cgdata.failure = false;
#ifdef USE_PCG
  gn2 = pkn_ScalarProductd ( nvars, grad, grad );
#ifdef G2MBL_TIME_IT
  pkv_Tic ( NULL );
#endif
  memset ( dcoeff, 0, nvars*sizeof(double) );
  result = pkn_PCGd ( nvars, &cgdata, grad, dcoeff,
                      _g2mbl_MLmultAxd, _g2mbl_MLmultQIxd,
                      MAXCGN, CG_EPS*gn2, CG_DELTA*gn2, &itm );
#ifdef G2MBL_TIME_IT
  d->time_cg += pkv_Toc ( NULL );
#endif
#ifdef _DEBUG
printf ( " %d ", itm );
#endif

#else
  result = _g2mbl_MLmultQIxd ( nvars, &cgdata, grad, dcoeff );
#endif
  if ( cgdata.failure )
    return -1.0;
  if ( result ) {
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "F" );
#endif
    f = g2mbl_MLFuncd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
                        ndomel, domelind, domelem, domelcpind, true, ftab );
#ifdef _DEBUG
printf ( " nu = %10g, f = %10g\n", nu, f );
#endif
  }
  else {
#ifdef _DEBUG
printf ( " nu = %10g, f = inf\n", nu );
#endif
    f = MYINFINITY;
  }
  return f;
} /*_g2mbl_MLCAuxNuFuncd*/

static boolean _g2mbl_H3x3FindEigenvalueIntervald ( int nvcp, int nHbl,
                       nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl,
                       double *aux, double *a, double *b )
{
  int    i, j, j0, j1;
  double *bH, *x1, *x2, aa, bb, c;

        /* the array aux is used to compute the radii */
        /* of the Gershgorin circles */
  memset ( aux, 0, 3*nvcp*sizeof(double) );
  for ( i = 0; i < nvcp; i++ ) {
    j0 = iHbl[i].firsthbl;
    j1 = j0 + iHbl[i].nhbl - 1;
    x1 = &aux[3*i];
    for ( j = j0; j < j1; j++ ) {
      bH = &Hbl[tHbl[j]];
      x2 = &aux[3*cHbl[j]];
          /* block below the diagonal */
      x1[0] += fabs(bH[0]) + fabs(bH[1]) + fabs(bH[2]);
      x1[1] += fabs(bH[3]) + fabs(bH[4]) + fabs(bH[5]);
      x1[2] += fabs(bH[6]) + fabs(bH[7]) + fabs(bH[8]);
          /* block above the diagonal */
      x2[0] += fabs(bH[0]) + fabs(bH[3]) + fabs(bH[6]);
      x2[1] += fabs(bH[1]) + fabs(bH[4]) + fabs(bH[7]);
      x2[2] += fabs(bH[2]) + fabs(bH[5]) + fabs(bH[8]);
    }
          /* the diagonal block */
    bH = &Hbl[tHbl[j1]];
    x1[0] += fabs(bH[3]) + fabs(bH[6]);
    x1[1] += fabs(bH[3]) + fabs(bH[7]);
    x1[2] += fabs(bH[6]) + fabs(bH[7]);
  }
        /* now find the interval with all eigenvalues */
  aa = bb = Hbl[tHbl[iHbl[0].firsthbl]];
  for ( i = 0; i < nvcp; i++ ) {
    j = iHbl[i].firsthbl + iHbl[i].nhbl - 1;
    bH = &Hbl[tHbl[j]];
    x1 = &aux[3*i];
    c = bH[0] - x1[0];  aa = min ( aa, c );
    c = bH[0] + x1[0];  bb = max ( bb, c );
    c = bH[4] - x1[1];  aa = min ( aa, c );
    c = bH[4] + x1[1];  bb = max ( bb, c );
    c = bH[8] - x1[2];  aa = min ( aa, c );
    c = bH[8] + x1[2];  bb = max ( bb, c );
  }
  *a = aa;
  *b = bb;
  return true;
} /*_H3x3FindEigenvalueIntervald*/

boolean g2mbl_MLOptBlockCd ( void *data, int bl )
{
  void            *sp;
  mesh_ml_optdata *d;
  mlblock_desc    *bd;
  point3d         *mvcp, *auxmvcp;
  int             nv, nvcp, nvars, ndomel, *vncpi, *domelind, *domelcpind;
  meshdom_elem    *domelem;
  double          *ftab1, *ftab2, *gtab1, *gtab2, *htab, *coeff, *dcoeff, *grad;
  double          **Nitabs, **Jac, *qcoeff, *Hbl;
  int             nkn, nHbl, *cHbl, *tHbl;
  nzHbl_rowdesc   *iHbl;
  int             i, j, k, kk, itm;
  double          lco, ldco, func, gn, gna, eps, f0, gn0;
  double          lmin, lmax, ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  mesh_ml_cg_data cgdata;
  boolean         recalc, fg_ok, advance, in_this, cg_ok, positive;

#define MLFUNC(nu) \
_g2mbl_MLCAuxNuFuncd ( d, bl, nkn, qcoeff, Nitabs, Jac, nv, mvcp, nvcp, vncpi, \
        ndomel, domelind, domelem, domelcpind, ftab2, \
        nHbl, iHbl, cHbl, tHbl, Hbl, grad, dcoeff, auxmvcp, nu )

#define RECORD_MIN(g,fnu) \
  { if ( fnu < 0.0 ) \
      goto failure; \
    if ( fnu < fge ) { \
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
  Hbl        = d->Hbl;
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
  nHbl     = bd->nHbl;
  iHbl     = bd->iHbl;
  cHbl     = bd->cHbl;
  tHbl     = bd->tHbl;
  bd->fghflag &= ~FLAG_ADVANCE;
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvars );
  if ( !coeff ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  dcoeff = &coeff[nvars];
  grad = &dcoeff[nvars];
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

  ga = fgc = 0.0;
  if ( d->log_level > 1 ) {
    f0 = gn0 = gna = MYINFINITY;
  }
/*  d->nu[0] = -1.0; */
  fg_ok = advance = false;
  for ( k = kk = 0; ; k++ ) {
    memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
    lco = 0.0;
    for ( i = 0; i < nvcp; i++ )
      lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
    lco = sqrt ( lco );

        /* compute the functional value, gradient and Hessian */
    if ( !(bd->fghflag & FLAG_H) && k == 0 &&
         !(bd->fghflag & FLAG_CH) && bl != d->currentblock )
      bd->fghflag = (bd->fghflag & ~FLAG_H) | (d->bd[(bl-1)/2].fghflag & FLAG_H);
    if ( !(bd->fghflag & FLAG_H) ) {
#ifdef _DEBUG
printf ( " H" );
#endif
#ifdef G2MBL_TIME_IT
      pkv_Tic ( NULL );
#endif
      if ( !g2mbl_MLFuncGradHessianAd ( d->nkn1, d->aqcoeff, d->aNitabs,
                     d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp,
                     ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                     true, ftab1, gtab1, htab,
                     &fga, dcoeff, nHbl, iHbl, cHbl, tHbl, Hbl ) ) {
        PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_19, ERRMSG_19 );
        goto failure;
      }
#ifdef G2MBL_TIME_IT
      d->time_h += pkv_Toc ( NULL );
#endif
      bd->fghflag |= FLAG_H | FLAG_CH;
      bd->fghflag &= ~(FLAG_CMH | FLAG_CMLH);
      in_this = true;
      if ( d->nkn1 == nkn ) {
        func = fga;
        memcpy ( grad, dcoeff, nvars*sizeof(double) );
        bd->fghflag |= FLAG_F | FLAG_G;
        fg_ok = true;
      }
      _g2mbl_MLInvalSmallBlocks ( d, 2*bl+1 );
      _g2mbl_MLInvalSmallBlocks ( d, 2*bl+2 );
    }
    else
      in_this = false;
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
          PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_20, ERRMSG_20 );
          goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
      fg_ok = true;
    }
    else if ( k == 0 ) {
      if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp,
                 ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                 false, ftab2, gtab2, &func, grad ) ) {
          PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_20, ERRMSG_20 );
          goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
      fg_ok = true;
    }
    gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
    if ( k == 0 ) {
      f0 = func;
      gn0 = gn;
    }

    cgdata.d = d;
    cgdata.bl = cgdata.cbl = bl;
    cgdata.nu = 0.0;
    cgdata.failure = false;
    if ( !_g2mbl_MLDecomposeBlockPrecond ( d, bl, 0.0, &positive ) )
      goto failure;
    if ( !positive ) {
#ifdef _DEBUG
printf ( "\n" );
#endif
      goto lm_trajectory;
    }
    memset ( dcoeff, 0, nvars*sizeof(double) );
#ifdef _DEBUG
printf ( "," );
#endif
#ifdef G2MBL_TIME_IT
    pkv_Tic ( NULL );
#endif
    cg_ok = pkn_PCGd ( nvars, &cgdata, grad, dcoeff,
                     _g2mbl_MLmultAxd, _g2mbl_MLmultQIxd,
                     MAXCGN, CG_EPS*gn*gn, CG_DELTA*gn*gn, &itm );
#ifdef G2MBL_TIME_IT
    d->time_cg += pkv_Toc ( NULL );
#endif
    if ( cgdata.failure )
      goto failure;
    if ( !cg_ok ) {
#ifdef _DEBUG
printf ( "\n" );
#endif
      goto lm_trajectory;
    }
#ifdef _DEBUG
printf ( " %d ", itm );
#endif
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    if ( ldco > BETA*lco ) { /* increment much too big */
      if ( k == 0 ) {
        if ( in_this )
          bd->fghflag |= FLAG_H;
        else
          bd->fghflag &= ~FLAG_H;
#ifdef _DEBUG
printf ( "\n" );
#endif
        break;
      }
      else {
        bd->fghflag &= ~FLAG_H;
        goto next_iter;
      }
    }
    else if ( ldco <= EPS3*lco ) {
      advance = true;
      goto next_iter;
    }
    _g2mbl_AddCPIncrement ( nvcp, vncpi, mvcp, (vector3d*)dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
    if ( !g2mbl_MLFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp,
               ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
               true, ftab2, gtab2, &fga, dcoeff ) ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_20, ERRMSG_20 );
      goto failure;
    }
    bd->fghflag &= ~(FLAG_F | FLAG_G);  /* until acceptance */
    fg_ok = false;

    if ( fga >= func ) {
        /* it may turn out that the maximal accuracy has been reached */
        /* and then the iterations must be terminated */
      if ( ldco <= EPS4*lco && (fga-func) <= DELTA4*func ) {
        advance = true;
        goto next_iter;
      }
      else if ( k > 0 ) {
        bd->fghflag &= ~FLAG_H;
        goto next_iter;
      }
      else {
        if ( in_this )
          bd->fghflag |= FLAG_H;
        else
          bd->fghflag &= ~FLAG_H;
#ifdef _DEBUG
printf ( "\n" );
#endif
        break;
      }
    }
    for ( j = 0; j < nvcp; j++ )
      mvcp[vncpi[j]] = auxmvcp[vncpi[j]];
    bd->fghflag |= FLAG_F | FLAG_G | FLAG_ADVANCE;  /* new point has been accepted */
    func = fga;
    memcpy ( grad, dcoeff, nvars*sizeof(double) );
    gna = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
    if ( ldco < eps*lco ) {
      fge = fga;
      fg_ok = true;
      advance = true;
      goto next_iter1;
    }
    if ( gna >= gn ) {
      if ( k > 0 ) {
        bd->fghflag &= ~FLAG_H;
        goto next_iter;
      }
      else {
        if ( in_this )
          bd->fghflag |= FLAG_H;
        else
          bd->fghflag &= ~FLAG_H;
#ifdef _DEBUG 
printf ( "\n" );
#endif 
        break;
      }
    }
    else if ( gna > 0.5*gn ) {
      bd->fghflag &= ~FLAG_H;
      if ( bl == d->currentblock )
        goto next_iter;
    }
    else if ( bl != 0 && gna < 0.25*gn ) {
      advance = true;
      goto next_iter;
    }
  }

lm_trajectory:
        /* minimization along the generalised Levenberg-Marquardt trajectory */
  fg_ok = false;
  ga = 0.0;  fga = MYINFINITY;
  gb = ge = MYINFINITY;
  fgb = fge = func;
  memset ( coeff, 0, nvars*sizeof(double) );
          /* get the initial nu for bracketing */
  if ( d->nu[0] <= 0.0 ) {
    if ( !_g2mbl_H3x3FindEigenvalueIntervald ( nvcp,
                     nHbl, iHbl, cHbl, tHbl, Hbl, dcoeff, &lmin, &lmax ) ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_23, ERRMSG_23 );
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
        /* find second point in the bracket to start the subdivision */
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
        /* subdivision */
#ifdef _DEBUG
printf ( "g\n" );
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
next_iter1:
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
      printf ( "\n func %10g -> %10g, gn %10g -> %10g", f0, fge, gn0, gna );
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
} /*g2mbl_MLOptBlockCd*/

/* ///////////////////////////////////////////////////////////////////////// */
#ifdef _DEBUG
static clock_t tic, tt;
#endif

boolean g2mbl_MLCOptIterd ( void *data, boolean *finished )
{
  void            *sp;
  mesh_ml_optdata *d;
  int             bl;

#ifdef _DEBUG
pkv_Tic ( &tic );
#endif
  sp = pkv_GetScratchMemTop ();
  d = (mesh_ml_optdata*)data;
  bl = d->currentblock;
  if ( bl < 0 || bl >= d->nblocks )
    goto failure;
#ifdef _DEBUG
printf ( "%3d:", bl );
#endif
  d->bd[bl].fghflag &= ~FLAG_H;
  if ( d->bd[bl].iHbl ) {  /* a "big" block */
    if ( !g2mbl_MLOptBlockCd ( data, bl ) )
      goto failure;
  }
  else {                   /* a "small", i.e. undivided block */
    if ( !g2mbl_MLOptBlockBd ( data, bl ) )
      goto failure;
  }
  *finished = d->currentblock == -1;

  pkv_SetScratchMemTop ( sp );
#ifdef _DEBUG
tt = pkv_Toc ( &tic );
printf ( " time = %5.2f\n", (double)pkv_Seconds ( tt ) );
#endif
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLCOptIterd*/

