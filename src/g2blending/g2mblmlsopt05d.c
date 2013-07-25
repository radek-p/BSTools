
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
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

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_MLSAddCPIncrement ( int nvcp, int *vncpi,
                                point3d *mvcp, vector3d *mvcpn,
                                double *incr, point3d *omvcp )
{
  int i, j;

  for ( i = 0; i < nvcp; i++ ) {
    j = vncpi[i];
    AddVector3Md ( &mvcp[j], &mvcpn[j], -incr[i], &omvcp[j] );
  }
} /*_g2mbl_MLSAddCPIncrement*/

boolean g2mbl_MLSOptBlockAd ( void *data, int bl )
{
  void            *sp;
  mesh_ml_optdata *d;
  mlblock_desc    *bd, *bd1, *bd2;
  point3d         *mvcp, *auxmvcp;
  vector3d        *mvcpn;
  int             nv, nvcp, ndomel, *vncpi, *domelind, *domelcpind;
  meshdom_elem    *domelem;
  double          *ftab1, *ftab2, *gtab1, *gtab2, *htab, *coeff, *dcoeff, *grad, *Hbl;
  double          **Nitabs, **Jac, *qcoeff;
  int             nkn, nHbl, *cHbl, *tHbl;
  nzHbl_rowdesc   *iHbl;
  int             i, j, k, l, b1, b2, itm;
  double          lco, ldco, func, gn, fga, gna, eps, gn0, f0;
  mesh_ml_cg_data cgdata;
  boolean         in_this, recalc, fg_ok, advance, positive;

  sp = pkv_GetScratchMemTop ();
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
        /* numbers of subblocks */
  b1 = 2*bl+1;
  b2 = b1+1;
  bd1 = &d->bd[b1];
  bd2 = &d->bd[b2];
        /* allocate arrays */
  coeff = pkv_GetScratchMemd ( 3*nvcp );
  if ( !coeff ) {
printf ( "%s\n", ERRMSG_0 );
    goto failure;
  }
  dcoeff = &coeff[nvcp];
  grad   = &dcoeff[nvcp];
  auxmvcp = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !auxmvcp ) {
printf ( "%s\n", ERRMSG_0 );
    goto failure;
  }

  if ( d->log_level > 1 ) {
    f0 = gn0 = gna = MYINFINITY;
  }
  advance = fg_ok = false;
        /* this loop is forced to execute once, */
        /* which seems to be the right thing */
  for ( l = 0; l < 1; l++ ) {
        /* the first stage is using the Newton method for the block,    */
        /* the system of equations with the Hessian is solved using     */
        /* the CG method; first it is necessary to compute the gradient */
    for ( k = 0; ; k++ ) {
      memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
      lco = 0.0;
      for ( i = 0; i < nvcp; i++ )
        lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
      lco = sqrt ( lco );

        /* compute the Hessian */
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
        if ( !g2mbl_MLSFuncGradHessianAd ( d->nkn1, d->aqcoeff,
                     d->aNitabs, d->aNijtabs, d->aMijtabs, d->aJac, nv, mvcp, mvcpn,
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
        _g2mbl_MLInvalSmallBlocks ( d, b1 );
        _g2mbl_MLInvalSmallBlocks ( d, b2 );
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
        if ( !g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs,
                     Jac, nv, mvcp, mvcpn,
                     ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                     recalc, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_20 );
          goto failure;
        }
        bd->fghflag |= FLAG_F | FLAG_G;
        fg_ok = true;
      }
      else if ( k == 0 ) {
        if ( !g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs,
                     Jac, nv, mvcp, mvcpn,
                     ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                     false, ftab2, gtab2, &func, grad ) ) {
printf ( "%s\n", ERRMSG_20 );
          goto failure;
        }
        bd->fghflag |= FLAG_F | FLAG_G;
        fg_ok = true;
      }
      else {  /* these have been computed in the previous loop execution */
        func = fga;
        memcpy ( grad, dcoeff, nvcp*sizeof(double) );
        gn = gna;
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
      memset ( dcoeff, 0, nvcp*sizeof(double) );
      if ( !_g2mbl_MLSDecomposeBlockPrecond ( d, bl, 0.0, &positive ) )
        goto failure;
      if ( !positive )
        break;
#ifdef _DEBUG
printf ( "," );
#endif
#ifdef G2MBL_TIME_IT
pkv_Tic ( NULL );   
#endif
      if ( !pkn_PCGd ( nvcp, &cgdata, grad, dcoeff,
                       (cg_mult)_g2mbl_MLSmultAxd, (cg_mult)_g2mbl_MLSmultQIxd,
                       MAXCGN, CG_EPS*gn*gn, CG_DELTA*gn*gn, &itm ) ) {
        if ( cgdata.failure )
          goto failure;
        break;
      }
#ifdef G2MBL_TIME_IT
d->time_cg += pkv_Toc ( NULL );
#endif

#ifdef _DEBUG
printf ( " %d ", itm );
#endif
      ldco = sqrt ( pkn_ScalarProductd ( nvcp, dcoeff, dcoeff ) );
      if ( ldco > BETA*lco ) {  /* increment much too big */
        if ( k == 0 ) {
          if ( in_this )
            bd->fghflag |= FLAG_H;
          else
            bd->fghflag &= ~FLAG_H;
          break;
        }
        else {
          bd->fghflag &= ~FLAG_H;
          goto way_out;
        }
      }
      else if ( ldco <= EPS3*lco ) {
        advance = true;
        goto way_out;
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
          goto way_out;
        }
        else if ( k > 0 ) {
          bd->fghflag &= ~FLAG_H;
          goto way_out;
        }
        else {
          if ( in_this )
            bd->fghflag |= FLAG_H;
          else
            bd->fghflag &= ~FLAG_H;
          break;
        }
      }
      for ( j = 0; j < nvcp; j++ )
        mvcp[vncpi[j]] = auxmvcp[vncpi[j]];
      bd->fghflag |= FLAG_F | FLAG_G | FLAG_ADVANCE;  /* new point has been accepted */
      fg_ok = true;
      if ( ldco < eps*lco ) {
        advance = true;
        goto way_out;
      }
      func = fga;
      gna = sqrt ( pkn_ScalarProductd ( nvcp, dcoeff, dcoeff ) );
      if ( gna >= gn ) {
        if ( k > 0 ) {
          bd->fghflag &= ~FLAG_H;
          goto way_out;
        }
        else {
          if ( in_this )
            bd->fghflag |= FLAG_H;
          else
            bd->fghflag &= ~FLAG_H;
          break;
        }
      }
      else if ( gna > 0.5*gn ) {
        bd->fghflag &= ~FLAG_H;
        if ( bl == d->currentblock )
          goto way_out;
      }
      else if ( bl != 0 && gna < 0.25*gn ) {
        advance = true;
        goto way_out;
      }
    }

        /* the second stage is optimizing the two subblocks until */
        /* the Hessian seems good enough, i.e. it has not been    */
        /* recomputed by the subblock optimization algorithm      */
        /* actually - this is done also only once                 */
    for ( k = 0; k < 1; k++ ) {
      fg_ok = false;
      memcpy ( auxmvcp, mvcp, nv*sizeof(point3d) );
      lco = 0.0;
      for ( i = 0; i < nvcp; i++ )
        lco += DotProduct3d ( &mvcp[vncpi[i]], &mvcp[vncpi[i]] );
      lco = sqrt ( lco );

        /* compute the Hessian */
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
                     &func, grad, nHbl, iHbl, cHbl, tHbl, Hbl ) ) {
printf ( "%s\n", ERRMSG_19 );
          goto failure;
        }
#ifdef G2MBL_TIME_IT
d->time_h += pkv_Toc ( NULL );
#endif
        bd->fghflag |= FLAG_H | FLAG_CH;
        bd->fghflag &= ~FLAG_LH;
        if ( nkn == d->nkn1 )
          bd->fghflag |= FLAG_F | FLAG_G;
        else
          bd->fghflag &= ~(FLAG_F | FLAG_G);
        _g2mbl_MLInvalSmallBlocks ( d, b1 );
        _g2mbl_MLInvalSmallBlocks ( d, b2 );
      }
#ifdef _DEBUG
printf ( " (%d:", b2 );
#endif
      if ( bd2->iHbl ) {
        if ( !g2mbl_MLSOptBlockAd ( data, b2 ) )
          goto failure;
      }
      else {
        if ( !g2mbl_MLSOptBlockBd ( data, b2 ) )
          goto failure;
      }
      if ( bd2->fghflag & FLAG_CH ) {
        bd->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H);
        bd1->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H);
      }
      if ( bd2->fghflag & FLAG_ADVANCE )
        bd->fghflag &= ~(FLAG_F | FLAG_G);
#ifdef _DEBUG
printf ( ")(%d:", b1 );
#endif
      if ( bd1->iHbl ) {
        if ( !g2mbl_MLSOptBlockAd ( data, b1 ) )
          goto failure;
      }
      else {
        if ( !g2mbl_MLSOptBlockBd ( data, b1 ) )
          goto failure;
      }
      if ( bd1->fghflag & FLAG_CH ) {
        bd->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H);
        bd2->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H);
      }
      if ( bd1->fghflag & FLAG_ADVANCE )
        bd->fghflag &= ~(FLAG_F | FLAG_G);
#ifdef _DEBUG
printf ( ")" );
#endif

      bd->fghflag |= (bd1->fghflag | bd2->fghflag) & FLAG_ADVANCE;
      if ( bl == d->currentblock &&
          !((bd1->fghflag | bd2->fghflag) & FLAG_ADVANCE) ) {
printf ( "%s\n", ERRMSG_22 );
        goto failure;
      }
      else if ( !(bd1->fghflag & FLAG_CH) && !(bd2->fghflag & FLAG_CH) ) {
        bd->fghflag &= ~(FLAG_F | FLAG_G | FLAG_H);
        break;
      }
    }

    if ( bl == d->currentblock ) {
      bd->fghflag &= ~FLAG_H;
      goto way_out;
    }
  }

way_out:
  if ( d->log_level > 1 ) {
    if ( bl == d->currentblock ) {
#ifdef _DEBUG
if ( !fg_ok )
  printf ( "G" );
#endif
      g2mbl_MLSFuncGradd ( nkn, qcoeff, Nitabs, Jac, nv, mvcp, mvcpn,
                          ndomel, domelind, domelem, domelcpind, nvcp, vncpi,
                          true, ftab2, gtab2, &func, grad );
      bd->fghflag |= FLAG_F | FLAG_G;
      gna = sqrt ( pkn_ScalarProductd ( nvcp, grad, grad ) );
#ifdef _DEBUG
      printf ( "\n func %12g -> %12g, gn %12g -> %12g", f0, func, gn0, gna );
#endif
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
} /*g2mbl_MLSOptBlockAd*/

