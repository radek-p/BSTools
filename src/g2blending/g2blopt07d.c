
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"

#include "g2blprivated.h"

#define DEBUG
#define _DEBUG
/*#define __DEBUG */
#define COUNT

/* ////////////////////////////////////////////////////////////////////////// */
#define EPS0     5.0e-10
#define EPS1     1.0e-8
#define MAXNTN  16
#define MAXGTN  16
#define MAXSTN  30
#define MAXBTN  10
#define MAXCTN  20
#define THR1     0.9
#define THR2     0.75
#define THR3     0.25
#define DELTA   0.01
#define RHO     0.05
#define SIGMA   0.1
#define GRTHR   0.01

/* ////////////////////////////////////////////////////////////////////////// */
/* private data structure with variables preserved between iterations         */
typedef struct {
    int     lastknotu, lastknotv;
    int     pitch, pitch1, pitch2;
    point3d *cp;
    int     nkn1, nkn2;
    int     nsq;
    int     nip, nvars, niip;
    int     hsize;
    double  C;
    boolean newpoint, accurate;
    int     itn;
    double  nu;

    double  func, gn0;

    char    *dirtypt, *dirtysq;
    int     *prof;
    point3d *hcp;

    double  *ftab, *gtab, *htab,
            *hessian, **hrows, *Lhessian, **Lhrows;

    double  *aqcoeff, *aNitab, *aNijtab, *aMijtab, *bqcoeff, *bNitab;
#ifdef COUNT
    int     nLM, nN, nF, nG, nH;
#endif
  } lmt_cloptdata;

/* ////////////////////////////////////////////////////////////////////////// */
void g2bl_ClosedOptLMTDeallocated ( void **data )
{
  lmt_cloptdata *d;

  if ( *data ) {
    d = *data;
    if ( d->dirtypt ) PKV_FREE ( d->dirtypt );
    if ( d->ftab )    PKV_FREE ( d->ftab );
    if ( d->aqcoeff ) PKV_FREE ( d->aqcoeff );

    PKV_FREE ( *data )
  }
} /*g2bl_ClosedOptLMTDeallocated*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2bl_ClosedInitBlSurfaceOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data )
{
  void          *sp;
  lmt_cloptdata *d;
  int           size, s1, s2, s3, s4;
  int           nsq, nip, niip, nvars, hsize, pitch1;
  char          *dirtypt, *dirtysq;
  int           *prof;
  double        *ftab, *gtab, *htab,
                *hessian, **hrows, *Lhessian, **Lhrows;
  double        *aqknots, *aqcoeff, *abf, *adbf, *addbf, *adddbf,
                *bqknots, *bqcoeff, *bbf, *bdbf, *bddbf, *bdddbf;
  double        *aNitab, *aNijtab, *aMijtab, *bNitab;

  sp = pkv_GetScratchMemTop ();

  *data = d = NULL;
  if ( lastknotu < 10 || lastknotv < 10 ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_5, ERRMSG_5 );
    goto failure;
  }
  if ( nkn1 < 3 || nkn1 > 10 || nkn2 < nkn1 || nkn2 > 10 ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_5, ERRMSG_5 );
    goto failure;
  }

  PKV_MALLOC ( *data, sizeof(lmt_cloptdata) )
  d = *data;
  if ( !d ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  memset ( d, 0, sizeof(lmt_cloptdata) );  /* clear all pointer fields */

        /* numbers of squares, unknown control points and unknown variables */
  d->nsq = nsq = (lastknotu-6)*(lastknotv-6);
  d->nip = nip = (lastknotu-6)*(lastknotv-9);
  d->niip = niip = (lastknotu-3)*(lastknotv-9);
  d->nvars = nvars = 3*nip;

  d->lastknotu = lastknotu;
  d->lastknotv = lastknotv;
  d->pitch = pitch;
  d->pitch1 = pitch1 = 3*(lastknotv-3);
  d->pitch2 = 3*(lastknotv-9);
  d->cp = cp;
  d->nkn1 = nkn1;
  d->nkn2 = nkn2;

        /* allocate memory */
  size = (niip+nsq)*sizeof(char) + nvars*sizeof(int) +
         (lastknotu-3)*(lastknotv-3)*sizeof(point3d);
  PKV_MALLOC ( dirtypt, size )
  if ( !dirtypt ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  d->dirtypt = dirtypt;
  d->dirtysq = dirtysq = &dirtypt[niip];
  d->prof = prof = (int*)&dirtysq[nsq];
  d->hcp = (point3d*)&prof[nvars];

  d->hsize = hsize = _g2bl_SetupClosedHessian1Profile ( lastknotu, lastknotv, prof );

#ifdef _DEBUG
printf ( "nip = %d,  nsq = %d,  n = %d\n", nip, nsq, nvars );
printf ( "hsize = %d\n", hsize );
printf ( "nkn1 = %d, nkn2 = %d\n", nkn1, nkn2 );
#endif

  size = (nsq*(1 + SQUAREGDS + SQUAREHDS) + 2*hsize)*sizeof(double) +
         2*nvars*sizeof(double*);
  PKV_MALLOC ( ftab, size );
  if ( !ftab ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  d->ftab = ftab;
  d->gtab = gtab = &ftab[nsq];
  d->htab = htab = &gtab[nsq*SQUAREGDS];
  d->hessian = hessian = &htab[nsq*SQUAREHDS];
  d->Lhessian = Lhessian = &hessian[hsize];
  d->hrows = hrows = (double**)&Lhessian[hsize];
  d->Lhrows = Lhrows = &hrows[nvars];
  pkn_NRBFindRowsd ( nvars, prof, hessian, hrows );
  pkn_NRBFindRowsd ( nvars, prof, Lhessian, Lhrows );

        /* evaluate the basis functions */
  s1 = g2bl_NiSize ( nkn1 );
  s2 = g2bl_NijSize ( nkn1 );
  s3 = g2bl_MijSize ( nkn1 );
  s4 = g2bl_NiSize ( nkn2 );
  size = (nkn1 + s1 + s2 + s3 + nkn2 + s4)*sizeof(double);
  PKV_MALLOC ( aqcoeff, size );
  if ( !aqcoeff ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  d->aqcoeff = aqcoeff;
  d->aNitab = aNitab = &aqcoeff[nkn1];
  d->aNijtab = aNijtab = &aNitab[s1];
  d->aMijtab = aMijtab = &aNijtab[s2];
  d->bqcoeff = bqcoeff = &aMijtab[s3];
  d->bNitab = bNitab = &bqcoeff[nkn2];
          /* at the knots of the quadrature of lower order */
  if ( !_g2bl_TabBasisFuncd ( nkn1, &aqknots, &aqcoeff,
                              &abf, &adbf, &addbf, &adddbf ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_11, ERRMSG_11 );
    goto failure;
  }
  memcpy ( d->aqcoeff, aqcoeff, nkn1*sizeof(double) );
  g2bl_TabNid ( nkn1, abf, adbf, addbf, adddbf, aNitab );
  g2bl_TabNijd ( nkn1, abf, adbf, addbf, aNijtab );
  g2bl_TabMijd ( nkn1, abf, adbf, addbf, adddbf, aMijtab );
          /* at the knots of the quadrature of higher order */
  if ( !_g2bl_TabBasisFuncd ( nkn2, &bqknots, &bqcoeff,
                              &bbf, &bdbf, &bddbf, &bdddbf ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_11, ERRMSG_11 );
    goto failure;
  }
  memcpy ( d->bqcoeff, bqcoeff, nkn2*sizeof(double) );
  g2bl_TabNid ( nkn2, bbf, bdbf, bddbf, bdddbf, bNitab );
  
        /* compute the functional factor */
          /* square of the domain diameter */
  if ( dO <= 0.0 )
    dO = (lastknotu-6)*(lastknotu-6) + (lastknotv-6)*(lastknotv-6);
          /* square of the surface diameter */
  if ( dM <= 0.0 )
    dM = g2bl_ClosedSurfNetDiameterSqd ( lastknotu, lastknotv, pitch, cp );
  if ( dM <= 0.0 )
    goto failure;
          /* this is to make the construction invariant with homotetiae */
          /* of the domain and the surface control net */
  d->C = C*dO*dO/(dM*dM);
#ifdef _DEBUG
printf ( "C = %f\n", d->C );
#endif

        /* initialise data before the first iteration */
  d->newpoint = true;
  d->accurate = nkn1 == nkn2;
  memset ( dirtysq, DIRTY_HESS, nsq );
  pkv_Selectd ( lastknotu-3, pitch1, pitch, pitch1, &cp[0].x, &d->hcp[0].x );
  d->itn = 0;
  d->func = MYINFINITY;
  d->nu = -1.0;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2bl_OptLMTDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_ClosedInitBlSurfaceOptLMTd*/

boolean g2bl_ClosedIterBlSurfaceOptLMTd ( void *data, boolean *finished )
{
  void          *sp;
  lmt_cloptdata *d;
  point3d       *cp, *acp, *hcp;
  char          *dirtypt, *dirtysq;
  double        *coeff, *dcoeff, *grad;
  double        *aqcoeff, *aNitab, *aNijtab, *aMijtab, *qcoeff, *Nitab;
  double        *ftab, *gtab, *htab;
  int           *prof;
  double        *hessian, **hrows, *Lhessian, **Lhrows;
  int           lastknotu, lastknotv, pitch, pitch1, pitch2;
  int           nkn, nkn1, nvars, nsq, hsize;
  double        func, gn, gn1, lco, ldco, df, dq, rk;
  int           ktn, i;
  double        g0, g1, g2, ga, gb, gc, gd, fg1, fga, fgb, fgc, fgd;
  boolean       positive, all;

#define MLFUNC(nu) \
_g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian, \
                     Lhessian, Lhrows, nu, lastknotu, lastknotv, pitch, cp, acp, \
                     grad, dcoeff, d->C, ftab )

#ifdef _DEBUG
#define SUCCESS \
  { \
printf ( "\nfunc = %12.7e, gn = %12.7e", func, gn ); \
    *finished = true; \
    goto next_iter; \
  }
#else
#define SUCCESS \
  { \
    *finished = true; \
    goto next_iter; \
  }
#endif

#ifdef _DEBUG
#define SWITCH_ACCURACY \
  { \
printf ( "*" ); \
    d->accurate = d->newpoint = true; \
    d->func = MYINFINITY;  /* this is to prevent effects of changing */ \
                           /* quadrature approximation error */ \
    goto next_iter; \
  }
#else
#define SWITCH_ACCURACY \
  { \
    d->accurate = d->newpoint = true; \
    d->func = MYINFINITY;  /* this is to prevent effects of changing */ \
                           /* quadrature approximation error */ \
    goto next_iter; \
  }
#endif

#define RECORD_MIN(g,fnu) \
  { if ( fnu < d->func ) {\
      memcpy ( coeff, dcoeff, nvars*sizeof(double) ); \
      d->func = fnu; \
      d->nu = g; \
      d->newpoint = true; \
    } \
  }

  sp = pkv_GetScratchMemTop ();
        /* extract data to local variables */
  d = data;
#ifdef _DEBUG
printf ( "%2d: ", d->itn );
#endif
  fga = fgb = fgc = d->func;
  lastknotu = d->lastknotu;
  lastknotv = d->lastknotv;
  pitch = d->pitch;
  pitch1 = d->pitch1;
  pitch2 = d->pitch2;
  cp = d->cp;
  hcp = d->hcp;
  nvars = d->nvars;
  nsq = d->nsq;
  nkn1 = d->nkn1;
  aNitab = d->aNitab;
  aNijtab = d->aNijtab;
  aMijtab = d->aMijtab;
  aqcoeff = d->aqcoeff;
  if ( d->accurate ) {
    nkn = d->nkn2;
    Nitab = d->bNitab;
    qcoeff = d->bqcoeff;
  }
  else {
    nkn = d->nkn1;
    Nitab = d->aNitab;
    qcoeff = d->aqcoeff;
  }
  dirtypt = d->dirtypt;
  dirtysq = d->dirtysq;
  prof = d->prof;
  hessian = d->hessian;
  hrows = d->hrows;
  Lhessian = d->Lhessian;
  Lhrows = d->Lhrows;
  ftab = d->ftab;
  gtab = d->gtab;
  htab = d->htab;
  hsize = d->hsize;

        /* allocate necessary arrays */
  acp = (point3d*)pkv_GetScratchMemd ( pitch1*(lastknotu-3) );
  grad = pkv_GetScratchMemd ( nvars );
  coeff = pkv_GetScratchMemd ( nvars );
  dcoeff = pkv_GetScratchMemd ( nvars );
  if ( !acp || !grad || !coeff || !dcoeff ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

  pkv_Selectd ( lastknotu-3, pitch1, pitch, pitch1, &cp[0].x, &acp[0].x );
  lco = sqrt ( pkn_ScalarProductd ( pitch1*(lastknotu-3),
               &acp[0].x, &acp[0].x ) );
  if ( d->itn ) {
    if ( !_g2bl_ClosedLazyHessiand ( lastknotu, lastknotv, d->niip, nsq,
                                     acp, hcp, dirtypt, dirtysq, &all ) )
      goto failure;
  }
  else
    all = true;

#ifdef __DEBUG
g2bl_ClosedUFuncGradd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                        NULL, d->C, ftab, gtab, &func, grad );
#endif

  if ( d->newpoint ) {
        /* evaluate the functional, gradient and Hessian */
#ifdef _DEBUG
printf ( "H" );
#endif
#ifdef COUNT
    d->nH ++;
#endif
    g2bl_ClosedUFuncGradHessiand ( nkn1, aqcoeff, aNitab, aNijtab, aMijtab,
                                   lastknotu, lastknotv, pitch1, hcp, dirtysq,
                                   d->C, ftab, gtab, htab,
                                   &func, grad, hsize, prof, hrows );
    if ( !all || nkn != nkn1 )
      g2bl_ClosedUFuncGradd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                              NULL, d->C, ftab, gtab, &func, grad );
    d->newpoint = false;
  }
  else {
    if ( d->accurate ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, 14, ERRMSG_14 );
#ifdef _DEBUG
func = d->func;
gn = 1.0/0.0;
#endif
      SUCCESS          /* perhaps no success, but the algorithm */
                       /* failed to improve the point */
    }
    else
      SWITCH_ACCURACY
  }
  d->func = min ( func, d->func );
  gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
  if ( !d->itn )
    d->gn0 = gn;
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
        /* first, try the Newton method */
  memcpy ( Lhessian, hessian, hsize*sizeof(double) );
  positive = pkn_NRBSymCholeskyDecompd ( nvars, prof, Lhessian, Lhrows, NULL );
  if ( positive ) {
#ifdef _DEBUG
printf ( "+" );
#endif
    pkn_NRBLowerTrSolved ( nvars, prof, Lhessian, Lhrows, 1, 1, grad, 1, dcoeff );
    pkn_NRBUpperTrSolved ( nvars, prof, Lhessian, Lhrows, 1, 1, dcoeff, 1, dcoeff );
    pkn_SubtractMatrixd ( lastknotu-6, pitch2, pitch, &cp[3].x,
                          pitch2, dcoeff, pitch1, &acp[3].x );
    memcpy ( &acp[(lastknotu-6)*(lastknotv-3)], acp, 3*pitch1*sizeof(double) );
    ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
    if ( !d->accurate ) {
      if ( ldco <= 5.0*EPS1*lco )
        SWITCH_ACCURACY
    }
    else if ( ldco <= 5.0*EPS1*lco )
      SUCCESS

#ifdef _DEBUG
printf ( "F" );
#endif
#ifdef COUNT
    d->nF ++;
#endif
    func = g2bl_UFuncd ( nkn, qcoeff, Nitab,
                         lastknotu, lastknotv, pitch1,
                         acp, NULL, d->C, ftab );
    if ( func >= d->func )
      goto switch_to_lmt;

#ifdef COUNT
    d->nN ++;
#endif
    pkv_Selectd ( lastknotu-3, pitch2, pitch1, pitch, &acp[3], &cp[3].x );
    df = func-d->func;
    d->func = func;
    d->newpoint = true;
          /* test if the functional value decrement is sufficient */
    if ( !_g2bl_ComputeDeltaQd ( nvars, prof, hrows, grad, dcoeff, &dq ) )
      goto failure;
    rk = df/dq;
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
    d->nG ++;
#endif
    g2bl_ClosedUFuncGradd ( nkn, qcoeff, Nitab,
                            lastknotu, lastknotv, pitch1, acp, NULL,
                            d->C, ftab, gtab, &func, grad );
    gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
    if ( !d->accurate ) {
      if ( gn < 5.0*EPS0*d->gn0 )
        SWITCH_ACCURACY
    }
    else if ( gn <= EPS0*d->gn0 )
      SUCCESS

    if ( rk > THR1 ) {
        /* make additional iterations of the Newton method with the same */
        /* Hessian matrix */
      for ( ktn = 0; ktn < MAXNTN; ktn++ ) {
#ifdef _DEBUG
printf ( "," );
#endif
        gn1 = gn;
        pkn_NRBLowerTrSolved ( nvars, prof, Lhessian, Lhrows, 1, 1, grad, 1, dcoeff );
        pkn_NRBUpperTrSolved ( nvars, prof, Lhessian, Lhrows, 1, 1, dcoeff, 1, dcoeff );
        pkn_SubtractMatrixd ( lastknotu-6, pitch2, pitch, &cp[3].x,
                              pitch2, dcoeff, pitch1, &acp[3].x );
        memcpy ( &acp[(lastknotu-6)*(lastknotv-3)], acp, 3*pitch1*sizeof(double) );
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
        d->nG ++;
#endif
        g2bl_ClosedUFuncGradd ( nkn, qcoeff, Nitab, lastknotu, lastknotv,
                                pitch1, acp, NULL, d->C, ftab, gtab,
                                &func, grad );
        gn = sqrt ( pkn_ScalarProductd ( nvars, grad, grad ) );
        if ( func < d->func ) {
          d->func = func;
          d->newpoint = true;
          pkv_Selectd ( lastknotu-3, 3*(lastknotv-9), 3*(lastknotv-3),
                        pitch, &acp[3], &cp[3].x );
          if ( !d->accurate ) {
            if ( gn <= 5.0*EPS0*d->gn0 )
              SWITCH_ACCURACY
          }
          else if ( gn <= 5.0*EPS0*d->gn0 )
            SUCCESS
          ldco = sqrt ( pkn_ScalarProductd ( nvars, dcoeff, dcoeff ) );
          if ( !d->accurate ) {
            if ( ldco <= 5.0*EPS1*lco )
              SWITCH_ACCURACY
          }
          else if ( ldco <= EPS1*lco )
            SUCCESS
        }
        else {
            /* there may be func > d->func due to rounding errors, */
            /* hence a way out is necessary */
          if ( !d->accurate ) {
            if ( gn <= 5.0*EPS1*d->gn0 )
              SWITCH_ACCURACY
          }
          else if ( gn <= EPS0*d->gn0 )
            SUCCESS
        }
        if ( gn > THR3*gn1 )
          goto next_iter;
      }
    }
    goto next_iter;
  }
#ifdef _DEBUG
else
  printf ( "!" );
#endif
switch_to_lmt:
#ifdef COUNT
    d->nLM ++;
#endif
        /* minimization along the Levenberg-Marquardt trajectory */
          /* bracketing */
  if ( d->nu <= 0.0 ) {
    if ( !pkn_NRBSymFindEigenvalueIntervald ( nvars, prof, hessian, hrows,
                                              &gc, &gd ) )
      goto failure;
    g1 = 0.01*(gd-gc);
  }
  else
    g1 = d->nu;

  g0 = 0.0;
  for ( i = 0; ; i++ ) {
    fg1 = _g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian,
                   Lhessian, Lhrows, g1, lastknotu, lastknotv, pitch, cp, acp,
                   grad, dcoeff, d->C, ftab );
    RECORD_MIN ( g1, fg1 )
    if ( fg1 < 10.0*func )
      break;
    g1 *= (double)(i+2);
  }
  ga = gb = gc = gd = g1;
  fga = fgc = fgd = fg1;
  do {
    g1 = (1.0-0.01)*g0 + 0.01*ga;
    g2 = sqrt ( g0*ga );
    g1 = max ( g1, g2 );
    fg1 = _g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian,
                   Lhessian, Lhrows, g1, lastknotu, lastknotv, pitch, cp, acp,
                   grad, dcoeff, d->C, ftab );
    RECORD_MIN ( g1, fg1 )
    if ( fg1 < 10.0*func ) {
      gb = gd;  fgb = fgd;
      gd = gc;  fgd = fgc;
      gc = ga;  fgc = fga;
      ga = g1;  fga = fg1;
    }
    else
      g0 = g1;
  } while ( fga <= fgc );

        /* bracketing - second stage */
  for ( i = 0; gc == gd || fgb <= fgd; i++ ) {
    if ( fgc > fgd ) { ga = gc;  fga = fgc; }
    gc = gd;  fgc = fgd;
    gd = gb;  fgd = fgb;
    gb *= (double)(i+2);  /* trying the Fibonacci sequence */
    if ( fgd <= fgc ) {
      fgb = _g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian,
                     Lhessian, Lhrows, gb, lastknotu, lastknotv, pitch, cp, acp,
                     grad, dcoeff, d->C, ftab );
      RECORD_MIN ( gb, fgb )
    }
    else
      fgb *= 2.0;  /* will be rejected at once */
  }

          /* golden ratio minimization */
#ifdef _DEBUG
printf ( "g" );
#endif
  for ( ktn = 0; ktn < MAXGTN; ktn++ ) {
/*
    if ( fgc < fgd ) {
      gb = gd;  fgb = fgd;
      gd = gc;  fgd = fgc;
      gc = GRDIV ( ga, gd );
      fgc = _g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian,
                     Lhessian, Lhrows, gc, lastknotu, lastknotv, pitch, cp, acp,
                     grad, dcoeff, d->C, ftab );
      RECORD_MIN ( gc, fgc )
    }
    else {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = GRDIV ( gb, gc );
      fgd = _g2bl_ClosedAuxNuFuncd ( nkn, qcoeff, Nitab, nvars, hsize, prof, hessian,
                     Lhessian, Lhrows, gd, lastknotu, lastknotv, pitch, cp, acp,
                     grad, dcoeff, d->C, ftab );
      RECORD_MIN ( gd, fgd )
    }
*/
    if ( _g2mbl_DivideIntervald ( &ga, &gc, &gd, &gb,
                                  &fga, &fgc, &fgd, &fgb ) ) {
      fgc = MLFUNC ( gc );
      RECORD_MIN ( gc, fgc )
    }
    else {
      fgd = MLFUNC ( gd );
      RECORD_MIN ( gd, fgd )
    }
    if ( fga-d->func < GRTHR*(func-fgb) )
      break;
  }
  pkn_SubtractMatrixd ( lastknotu-6, pitch2, pitch, &cp[3].x,
                        pitch2, coeff, pitch, &cp[3].x );
  pkv_Selectd ( 3, pitch2, pitch, pitch,
                &cp[3].x, &cp[(lastknotu-6)*pitch/3+3].x );

next_iter:
  d->itn ++;
#ifdef _DEBUG
printf ( "\n" );
#endif
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_ClosedIterBlSurfaceOptLMTd*/

boolean g2bl_ClosedFindBlSurfaceLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 )
{
  void          *sp;
  lmt_cloptdata *data;
  int           itn;
  boolean       finished;

  sp = pkv_GetScratchMemTop ();

  data = NULL;
  if ( !g2bl_ClosedInitBlSurfaceOptLMTd ( lastknotu, lastknotv, pitch, cp,
                                          C, dO, dM, nkn1, nkn2, (void**)&data ) )
    goto failure;

  finished = false;
  for ( itn = 0; itn < maxit; itn++ ) {
    if ( !g2bl_ClosedIterBlSurfaceOptLMTd ( data, &finished ) )
      goto failure;
    if ( finished )
      break;
  }

#ifdef COUNT
printf ( "nF = %d, nG = %d, nH = %d, nN = %d, nLM = %d\n",
         data->nF, data->nG, data->nH, data->nN, data->nLM );
#endif
  g2bl_ClosedOptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2bl_OptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_ClosedFindBlSurfaceLMTd*/

