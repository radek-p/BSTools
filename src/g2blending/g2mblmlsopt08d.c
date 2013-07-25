
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
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
static double _g2mbl_MLSCAuxNuFuncd ( mesh_ml_optdata *d, int bl,
               int nkn, double *qcoeff, double **Nitabs,
               double **Jac, int nv, point3d *mvcp, vector3d *mvcpn,
               int nvcp, int *vncpi,
               int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
               double *ftab,
               int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl,
               double *grad, double *dcoeff, point3d *auxmvcp, double nu )
{
  int             itm;
  double          f, gn2;
  mesh_ml_cg_data cgdata;
  boolean         result;

  _g2mbl_MLInvalSmallBlocks ( d, 2*bl+1 );
  _g2mbl_MLInvalSmallBlocks ( d, 2*bl+2 );
  d->bd[bl].fghflag &= ~(FLAG_CMH | FLAG_CMLH);
  if ( !_g2mbl_MLSDecomposeBlockPrecond ( d, bl, nu, &result ) )
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
  gn2 = pkn_ScalarProductd ( nvcp, grad, grad );
#ifdef G2MBL_TIME_IT
  pkv_Tic ( NULL );
#endif
  memset ( dcoeff, 0, nvcp*sizeof(double) );
  result = pkn_PCGd ( nvcp, &cgdata, grad, dcoeff,
                      (cg_mult)_g2mbl_MLSmultAxd, (cg_mult)_g2mbl_MLSmultQIxd,
                      MAXCGN, CG_EPS*gn2, CG_DELTA*gn2, &itm );
#ifdef G2MBL_TIME_IT
  d->time_cg += pkv_Toc ( NULL );
#endif
#ifdef _DEBUG
printf ( " %d ", itm );
#endif

#else
  result = _g2mbl_MLSmultQIxd ( nvcp, &cgdata, grad, dcoeff );
#endif
  if ( cgdata.failure )
    return -1.0;
  if ( result ) {
    _g2mbl_MLSAddCPIncrement ( nvcp, vncpi, mvcp, mvcpn, dcoeff, auxmvcp );
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
} /*_g2mbl_MLSCAuxNuFuncd*/

static boolean _g2mbl_H1x1FindEigenvalueIntervald ( int nvcp, int nHbl,
                       nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl,
                       double *aux, double *a, double *b )
{
  int    i, j, j0, j1;
  double aa, bb, c;

        /* the array aux is used to compute the radii */
        /* of the Gershgorin circles */
  memset ( aux, 0, nvcp*sizeof(double) );
  for ( i = 0; i < nvcp; i++ ) {
    j0 = iHbl[i].firsthbl;
    j1 = j0 + iHbl[i].nhbl - 1;
    for ( j = j0; j < j1; j++ ) {
      aa = fabs ( Hbl[tHbl[j]] );
          /* coefficient below the diagonal */
      aux[i] += aa;
          /* coefficient above the diagonal */
      aux[cHbl[j]] += aa;
    }
  }
        /* now find the interval with all eigenvalues */
  aa = bb = Hbl[tHbl[iHbl[0].firsthbl]];
  for ( i = 0; i < nvcp; i++ ) {
    j = tHbl[iHbl[i].firsthbl + iHbl[i].nhbl - 1];
    c = Hbl[j] - aux[i];  aa = min ( aa, c );
    c = Hbl[j] + aux[i];  bb = max ( bb, c );
  }
  *a = aa;
  *b = bb;
  return true;
} /*_g2mbl_H1x1FindEigenvalueIntervald*/

boolean g2mbl_MLSOptBlockCd ( void *data, int bl )
{
  void            *sp;
  mesh_ml_optdata *d;
  mlblock_desc    *bd;
  point3d         *mvcp, *auxmvcp;
  vector3d        *mvcpn;
  int             nv, nvcp, ndomel, *vncpi, *domelind, *domelcpind;
  meshdom_elem    *domelem;
  double          *ftab1, *ftab2, *gtab1, *gtab2, *htab, *coeff, *dcoeff, *grad;
  double          **Nitabs, **Jac, *qcoeff, *Hbl;
  int             nkn, nHbl, *cHbl, *tHbl;
  nzHbl_rowdesc   *iHbl;
  int             i, j, k, kk, itm;
  double          lco, ldco, func, gn, gna, eps, f0, gn0;
  double          lmin, lmax,ga, gb, gc, gd, ge, fga, fgb, fgc, fgd, fge, fm;
  mesh_ml_cg_data cgdata;
  boolean         recalc, fg_ok, advance, in_this, cg_ok, positive;

#define MLFUNC(nu) \
_g2mbl_MLSCAuxNuFuncd ( d, bl, nkn, qcoeff, Nitabs, Jac, nv, mvcp, mvcpn, \
        nvcp, vncpi, ndomel, domelind, domelem, domelcpind, ftab2, \
        nHbl, iHbl, cHbl, tHbl, Hbl, grad, dcoeff, auxmvcp, nu )

#define RECORD_MIN(g,fnu) \
  { if ( fnu < 0.0 ) \
      goto failure; \
    if ( fnu < fge ) { \
      memcpy ( coeff, dcoeff, nvcp*sizeof(double) ); \
      ldco = sqrt ( pkn_ScalarProductd ( nvcp, dcoeff, dcoeff ) ); \
      fge = fnu; \
      d->nu[1] = ge = g; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
        /* copy data to local variables */
  d = (mesh_ml_optdata*)data;
  nv         = d->nv;
  mvcp       = d->mvcp;
  mvcpn      = d->mvcpn;
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
  ndomel   = bd->ndomel;
  vncpi    = bd->vncpi;
  domelind = bd->domelind;
  nHbl     = bd->nHbl;
  iHbl     = bd->iHbl;
  cHbl     = bd->cHbl;
  tHbl     = bd->tHbl;
  bd->fghflag &= ~FLAG_ADVANCE;
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvcp );
  if ( !coeff ) {
printf ( "%s\n", ERRMSG_0 );
    goto failure;
  }
  dcoeff = &coeff[nvcp];
  grad = &dcoeff[nvcp];
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp ) {
printf ( "%s\n", ERRMSG_0 );
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
      if ( !g2mbl_MLSFuncGradHessianAd ( d->nkn1, d->aqcoeff, d->aNitabs,
                     d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp, mvcpn,
                     ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                     true, ftab1, gtab1, htab,
                     &fga, dcoeff, nHbl, iHbl, cHbl, tHbl, Hbl ) ) {
printf ( "%s\n", ERRMSG_19 );
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
        memcpy ( grad, dcoeff, nvcp*sizeof(double) );
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
      if ( !g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, mvcpn,
                 ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                 recalc, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_20 );
          goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
      fg_ok = true;
    }
    else if ( k == 0 ) {
      if ( !g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, mvcpn,
                 ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                 false, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_20 );
          goto failure;
      }
      bd->fghflag |= FLAG_F | FLAG_G;
      fg_ok = true;
    }
    gn = sqrt ( pkn_ScalarProductd ( nvcp, grad, grad ) );
    if ( k == 0 ) {
      f0 = func;
      gn0 = gn;
    }

    cgdata.d = d;
    cgdata.bl = cgdata.cbl = bl;
    cgdata.nu = 0.0;
    cgdata.failure = false;
    if ( !_g2mbl_MLSDecomposeBlockPrecond ( d, bl, 0.0, &positive ) )
      goto failure;
    if ( !positive ) {
#ifdef _DEBUG   
printf ( "\n" );
#endif
      goto lm_trajectory;
    }
    memset ( dcoeff, 0, nvcp*sizeof(double) );
#ifdef _DEBUG
printf ( "," );
#endif
#ifdef G2MBL_TIME_IT
    pkv_Tic ( NULL );
#endif
    cg_ok = pkn_PCGd ( nvcp, &cgdata, grad, dcoeff,
                     (cg_mult)_g2mbl_MLSmultAxd, (cg_mult)_g2mbl_MLSmultQIxd,
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
    ldco = sqrt ( pkn_ScalarProductd ( nvcp, dcoeff, dcoeff ) );
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
    _g2mbl_MLSAddCPIncrement ( nvcp, vncpi, mvcp, mvcpn, dcoeff, auxmvcp );
#ifdef _DEBUG
printf ( "G" );
#endif
    if ( !g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, auxmvcp, mvcpn,
               ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
               true, ftab2, gtab2, &fga, dcoeff ) ) {
printf ( "%s\n", ERRMSG_20 );
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
    memcpy ( grad, dcoeff, nvcp*sizeof(double) );
    gna = sqrt ( pkn_ScalarProductd ( nvcp, grad, grad ) );
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
  memset ( coeff, 0, nvcp*sizeof(double) );
          /* get the initial nu for bracketing */
  if ( d->nu[0] <= 0.0 ) {
    if ( !_g2mbl_H1x1FindEigenvalueIntervald ( nvcp,
                     nHbl, iHbl, cHbl, tHbl, Hbl, dcoeff, &lmin, &lmax ) ) {
printf ( "%s\n", ERRMSG_23 );
      goto failure;
    }
    d->nu[0] = d->nu[1] = gc = 0.001*(lmax-lmin);
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
  _g2mbl_MLSAddCPIncrement ( nvcp, vncpi, mvcp, mvcpn, coeff, mvcp );
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
        g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, mvcpn,
                            ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                            !fg_ok, ftab2, gtab2, &fge, grad );
        bd->fghflag |= FLAG_F | FLAG_G;
        gna = sqrt ( pkn_ScalarProductd ( nvcp, grad, grad ) );
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
} /*g2mbl_MLSOptBlockCd*/

/* ///////////////////////////////////////////////////////////////////////// */
#ifdef _DEBUG
static clock_t tic, tt;
#endif

boolean g2mbl_MLSCOptIterd ( void *data, boolean *finished )
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
    if ( !g2mbl_MLSOptBlockCd ( data, bl ) )
      goto failure;
  }
  else {                   /* a "small", i.e. undivided block */
    if ( !g2mbl_MLSOptBlockBd ( data, bl ) )
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
} /*g2mbl_MLSCOptIterd*/

