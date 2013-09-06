
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
#define COUNT

/* ////////////////////////////////////////////////////////////////////////// */
#define EPS0     5.0e-10
#define EPS1     1.0e-8
#define MAXNTN  16
#define MAXGTN  24
#define MAXSTN  30
#define MAXBTN  10
#define MAXCTN  20
#define THR1     0.9
#define THR2     0.75
#define THR3     0.25
#define DELTA   0.01
#define RHO     0.05
#define SIGMA   0.1
#define GRTHR   0.02

/* ////////////////////////////////////////////////////////////////////////// */
/* private data structure with variables preserved between iterations         */
typedef struct {
    int     lastknotu, lastknotv;
    int     pitch, pitch1, pitch2;
    point3d *cp;
    int     nkn1, nkn2;
    int     nsq;
    int     nip, nvars, nconstr, neqsw;
    int     hsize;
    double  C;
    boolean newpoint, accurate;
    int     itn;
    double  nu;

    double  func, gn0;

    char    *dirtypt, *dirtysq;
    int     *prof;
    point3d *hcp;
    double  *cmat, *acmat, *y, *y2, *D1;

    double  *ftab, *gtab, *htab,
            *hessian, **hrows;

    double  *aqcoeff, *aNitab, *aNijtab, *aMijtab,
            *bqcoeff, *bNitab;
#ifdef COUNT
    int     nLM, nN, nF, nG, nH;
#endif
  } lmt_constroptdata;

/* ////////////////////////////////////////////////////////////////////////// */
void g2bl_ConstrOptLMTDeallocated ( void **data )
{
  lmt_constroptdata *d;

  if ( *data ) {
    d = *data;
    if ( d->dirtypt ) PKV_FREE ( d->dirtypt );
    if ( d->ftab )    PKV_FREE ( d->ftab );
    if ( d->aqcoeff ) PKV_FREE ( d->aqcoeff );

    PKV_FREE ( *data )
  }
} /*g2bl_ConstrOptLMTDeallocated*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2bl_InitBlSurfaceConstrOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          int nconstr, double *constrmat,
                                          double *constrrhs,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data )
{
  void        *sp;
  lmt_constroptdata *d;
  int         size, s1, s2, s3, s4;
  int         i, j;
  int         ncp, nsq, nip, nvars, neqsw, hsize, pitch1, pitch2;
  char        *dirtypt, *dirtysq;
  int         *prof;
  double      *ftab, *gtab, *htab, *hessian, **hrows;
  double      *cmat, *acmat, *D1, *y, *coeff;
  double      *aqknots, *aqcoeff, *abf, *adbf, *addbf, *adddbf,
              *bqknots, *bqcoeff, *bbf, *bdbf, *bddbf, *bdddbf;
  double      *aNitab, *aNijtab, *aMijtab, *bNitab;

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
  if ( nconstr < 1 ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_5, ERRMSG_5 );
    goto failure;
  }

  PKV_MALLOC ( *data, sizeof(lmt_constroptdata) )
  d = *data;
  if ( !d ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  memset ( d, 0, sizeof(lmt_constroptdata) );  /* clear all pointer fields */

        /* numbers of squares, unknown control points and unknown variables */
  ncp = (lastknotu-3)*(lastknotv-3);
  d->nsq = nsq = (lastknotu-6)*(lastknotv-6);
  d->nip = nip = (lastknotu-9)*(lastknotv-9);
  d->nvars = nvars = 3*nip;
  d->nconstr = nconstr;
  d->neqsw = neqsw = nvars-nconstr;

  d->lastknotu = lastknotu;
  d->lastknotv = lastknotv;
  d->pitch = pitch;
  d->pitch1 = pitch1 = 3*(lastknotv-3);
  d->pitch2 = pitch2 = 3*(lastknotv-9);
  d->cp = cp;
  d->nkn1 = nkn1;
  d->nkn2 = nkn2;

        /* allocate memory */
  size = (nip+nsq)*sizeof(char) + nvars*sizeof(int) + ncp*sizeof(point3d) +
         (nconstr*(nvars+2) + nvars + (nconstr*(nconstr+1))/2)*sizeof(double);
  PKV_MALLOC ( dirtypt, size )
  if ( !dirtypt ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  d->dirtypt = dirtypt;
  d->dirtysq = dirtysq = &dirtypt[nip];
  d->prof = prof = (int*)&dirtysq[nsq];
  d->hcp = (point3d*)&prof[nvars];
  d->cmat = cmat = (double*)&d->hcp[ncp];  d->acmat = acmat = &cmat[nconstr*nvars];
  d->y = y = &cmat[nconstr*(nvars+2)];     d->y2 = &y[nconstr];
  d->D1 = D1 = &d->y[nvars];

  d->hsize = hsize = _g2bl_SetupHessian1Profile ( lastknotu, lastknotv, prof );

#ifdef _DEBUG
printf ( "nip = %d,  nsq = %d,  n = %d\n", nip, nsq, nvars );
printf ( "hsize = %d\n", hsize );
printf ( "nkn1 = %d, nkn2 = %d\n", nkn1, nkn2 );
#endif

  size = (nsq*(1 + SQUAREGDS + SQUAREHDS) + hsize)*sizeof(double) +
         nvars*sizeof(double*);
  PKV_MALLOC ( ftab, size );
  if ( !ftab ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  d->ftab = ftab;
  d->gtab = gtab = &ftab[nsq];
  d->htab = htab = &gtab[nsq*SQUAREGDS];
  d->hessian = hessian = &htab[nsq*SQUAREHDS];
  d->hrows = hrows = (double**)&hessian[hsize];
  pkn_NRBFindRowsd ( nvars, prof, hessian, hrows );

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
    dM = g2bl_SurfNetDiameterSqd ( lastknotu, lastknotv, pitch, cp );
  if ( dM <= 0.0 )
    goto failure;
          /* this is to make the construction invariant with homotetiae */
          /* of the domain and the surface control net */
  d->C = C*dO*dO/(dM*dM);
#ifdef _DEBUG
printf ( "C = %f\n", d->C );
#endif

        /* decompose the matrix of constraint equations */
  pkv_TransposeMatrixd ( nconstr, nvars, nvars, constrmat, nconstr, cmat );
  if ( !pkn_QRDecomposeMatrixd ( nvars, nconstr, cmat, acmat ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_15, ERRMSG_15 );
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cmat[i*nconstr+j];

        /* find y1 by solving D1*y1 = constrrhs */
  pkn_LowerTrMatrixSolved ( nconstr, D1, 1, 1, constrrhs, 1, y );

        /* extract the initial approximation of the unknown variables */
        /* from the control net */
  coeff = pkv_GetScratchMemd ( nvars );
  if ( !coeff ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pkv_Selectd ( lastknotu-9, pitch2, pitch, pitch2, &cp[pitch+3], coeff );
        /* change the coordinate system */
  pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
        /* enforce the constraints and go back to the original variables */
  memcpy ( coeff, y, nconstr*sizeof(double) );
  memcpy ( &y[nconstr], &coeff[nconstr], neqsw*sizeof(double) );
  pkn_multiInvReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
  pkv_Selectd ( lastknotu-9, pitch2, pitch2, pitch, coeff, &cp[pitch+3].x );

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
  g2bl_ConstrOptLMTDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_InitBlSurfaceConstrOptLMTd*/

/* ////////////////////////////////////////////////////////////////////////// */
static double _g2bl_AuxConstrNuFuncd ( int nknots,
                          const double *qcoeff, double *Nitab,
                          int neqs, int nconstr,
                          double *cmat, double *acmat,
                          int qsize, int *prof, double *QA22,
                          double *LM22, double **LM22rows, double nu,
                          int lastknotu, int lastknotv,
                          int pitch, point3d *cp, point3d *acp,
                          double *y, double *g2, double *coeff, double *dcoeff,
                          double tC, double *ftab )
{
  int neqsw, pitch1, pitch2;

  neqsw = neqs-nconstr;
  pitch1 = 3*(lastknotv-3);
  pitch2 = 3*(lastknotv-9);
  if ( _g2bl_ShiftDecompHessiand ( neqsw, qsize, prof, QA22,
                                   LM22, LM22rows, nu ) ) {
    pkv_Selectd ( lastknotu-9, pitch2, pitch, pitch2, &cp[pitch+3], coeff );
    pkn_multiReflectVectord ( neqs, nconstr, cmat, acmat, 1, 1, coeff );
    pkn_NRBLowerTrSolved ( neqsw, prof, LM22, LM22rows, 1, 1, g2,
                           1, dcoeff );
    pkn_NRBUpperTrSolved ( neqsw, prof, LM22, LM22rows, 1, 1, dcoeff,
                           1, dcoeff );
    pkn_SubtractMatrixd ( 1, neqsw, 0, &coeff[nconstr], 0, dcoeff,
                          0, &coeff[nconstr] );
    memcpy ( coeff, y, nconstr*sizeof(double) );
    pkn_multiInvReflectVectord ( neqs, nconstr, cmat, acmat, 1, 1, coeff );
    pkv_Selectd ( lastknotu-9, pitch2, pitch2, pitch1, coeff, &acp[pitch1+3].x );
#ifdef _DEBUG
printf ( "F" );
#endif
    return g2bl_UFuncd ( nknots, qcoeff, Nitab, lastknotu, lastknotv,
                         pitch1, acp, NULL, tC, ftab );
  }
  else
    return MYINFINITY;
} /*_g2bl_AuxConstrNuFuncd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2bl_IterBlSurfaceConstrOptLMTd ( void *data, boolean *finished )
{
  void              *sp;
  lmt_constroptdata *d;
  point3d           *cp, *acp, *hcp;
  char              *dirtypt, *dirtysq;
  double            C;
  double            *coeff, *coeff2, *dcoeff, *mcoeff, *grad, *g, *gr2, *y, *y2;
  double            *aqcoeff, *aNitab, *aNijtab, *aMijtab, *qcoeff, *Nitab;
  double            *ftab, *gtab, *htab;
  int               *prof, *qa11prof, *qa22prof, qsize;
  double            **hrows, **QA11rows, **QA22rows, *QA21,
                    *LM22, **LM22rows;
  double            *cmat, *acmat;
  int               lastknotu, lastknotv, pitch, pitch1, pitch2;
  int               nkn, nkn1, nvars, nconstr, neqsw, nip, nsq, hsize;
  double            func, gn, gn1, lco, ldco, df, dq, rk;
  int               ktn, i;
  double            g0, g1, g2, ga, gb, gc, gd, fg1, fga, fgb, fgc, fgd;
  boolean           positive, all;

#define MLFUNC(nu) \
_g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr, \
                     cmat, acmat, qsize, qa22prof, QA22rows[0], \
                     LM22, LM22rows, nu, lastknotu, lastknotv, pitch, cp, acp, \
                     y, gr2, coeff, dcoeff, C, ftab )

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
  { if ( fnu < d->func ) { \
      memcpy ( mcoeff, coeff, nvars*sizeof(double) ); \
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
  C = d->C;
  lastknotu = d->lastknotu;
  lastknotv = d->lastknotv;
  pitch = d->pitch;
  pitch1 = d->pitch1;
  pitch2 = d->pitch2;
  cp = d->cp;
  hcp = d->hcp;
  nsq = d->nsq;
  nip = d->nip;
  nvars = d->nvars;
  nconstr = d->nconstr;
  neqsw = d->neqsw;
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
  hrows = d->hrows;
  ftab = d->ftab;
  gtab = d->gtab;
  htab = d->htab;
  hsize = d->hsize;
  cmat = d->cmat;
  acmat = d->acmat;
  y = d->y;
  y2 = d->y2;

        /* allocate necessary arrays */
  acp = (point3d*)pkv_GetScratchMemd ( pitch1*(lastknotu-3) );
  grad = pkv_GetScratchMemd ( nvars );
  g = pkv_GetScratchMemd ( nvars );
  coeff = pkv_GetScratchMemd ( nvars );
  dcoeff = pkv_GetScratchMemd ( neqsw );
  mcoeff = pkv_GetScratchMemd ( nvars );
  QA11rows = pkv_GetScratchMem ( nvars*sizeof(double*) );
  qa11prof = pkv_GetScratchMemi ( nvars );
  LM22rows = pkv_GetScratchMem ( neqsw*sizeof(double*) );
  if ( !acp || !grad || !g || !coeff || !dcoeff || !mcoeff ||
       !QA11rows || !qa11prof || !LM22rows ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  coeff2 = &coeff[nconstr];
  QA11rows[0] = NULL;
  QA22rows = &QA11rows[nconstr];
  qa22prof = &qa11prof[nconstr];
  gr2 = &g[nconstr];

  pkv_Selectd ( lastknotu-3, pitch1, pitch, pitch1, &cp[0].x, &acp[0].x );
  lco = sqrt ( pkn_ScalarProductd ( pitch1*(lastknotu-3),
               &acp[0].x, &acp[0].x ) );
  pkv_Selectd ( lastknotu-9, pitch2, pitch, pitch2, &cp[pitch+3].x, coeff );
  pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
  memcpy ( y2, coeff2, neqsw*sizeof(double) );

  if ( d->itn ) {
    if ( !_g2bl_LazyHessiand ( lastknotu, lastknotv, nip, nsq,
                               acp, hcp, dirtypt, dirtysq, &all ) )
      goto failure;
  }
  else
    all = true;

  if ( d->newpoint ) {
        /* evaluate the functional, gradient and Hessian */
#ifdef _DEBUG
printf ( "H" );
#endif
#ifdef COUNT
    d->nH ++;
#endif
    g2bl_UFuncGradHessiand ( nkn1, aqcoeff, aNitab, aNijtab, aMijtab,
                             lastknotu, lastknotv, pitch1, hcp, dirtysq,
                             C, ftab, gtab, htab,
                             &func, grad, hsize, prof, hrows );
    if ( !all || nkn != nkn1 )
      g2bl_UFuncGradd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                        NULL, C, ftab, gtab, &func, grad );

    d->newpoint = false;
  }
  else {
    if ( d->accurate ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_14, ERRMSG_14 );
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
       /* variable elimination */
#ifdef _DEBUG
printf ( "e" );
#endif
  if ( !pkn_NRBComputeQTSQbld ( nvars, prof, hrows[0], hrows,
                                nconstr, cmat, acmat,
                                qa11prof, QA11rows, qa22prof, QA22rows,
                                &QA21  ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_16, ERRMSG_16 );
    goto failure;
  }
  qsize = pkn_NRBArraySize ( neqsw, qa22prof );
#ifdef _DEBUG
if ( !d->itn )
  printf ( " qsize = %d", qsize );
#endif

        /* find g = Q^T*grad */
  memcpy ( g, grad, nvars*sizeof(double) );
  pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, g );

  d->func = min ( func, d->func );
  gn = sqrt ( pkn_ScalarProductd ( neqsw, gr2, gr2 ) );
  if ( !d->itn )
    d->gn0 = gn;
#ifdef _DEBUG
printf ( " func = %12.7e, gn = %12.7e\n", func, gn );
#endif
        /* first, try the Newton method */
  LM22 = pkv_GetScratchMemd ( qsize );
  if ( !LM22 ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  if ( !pkn_NRBFindRowsd ( neqsw, qa22prof, LM22, LM22rows ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_17, ERRMSG_17 );
    goto failure;
  }
  memcpy ( LM22, QA22rows[0], qsize*sizeof(double) );
  positive = pkn_NRBSymCholeskyDecompd ( neqsw, qa22prof, LM22, LM22rows, NULL );
  if ( positive ) {
#ifdef _DEBUG
printf ( "+" );
#endif
    pkn_NRBLowerTrSolved ( neqsw, qa22prof, LM22, LM22rows,
                           1, 1, gr2, 1, dcoeff );
    pkn_NRBUpperTrSolved ( neqsw, qa22prof, LM22, LM22rows,
                           1, 1, dcoeff, 1, dcoeff );
    ldco = sqrt ( pkn_ScalarProductd ( neqsw, dcoeff, dcoeff ) );
    pkv_Selectd ( lastknotu-9, pitch2, pitch, pitch2, &cp[pitch+3].x, coeff );
    pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
    memcpy ( coeff, y, nconstr*sizeof(double) );
    pkn_SubtractMatrixd ( 1, neqsw, 0, coeff2, 0, dcoeff, 0, coeff2 );
    pkn_multiInvReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
    pkv_Selectd ( lastknotu-9, pitch2, pitch2, pitch1, coeff, &acp[pitch1+3].x );

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
    func = g2bl_UFuncd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch1,
                         acp, NULL, C, ftab );
    if ( func >= d->func )
      goto switch_to_lmt;

#ifdef COUNT
    d->nN ++;
#endif
    pkv_Selectd ( lastknotu-9, pitch2, pitch1,
                  pitch, &acp[pitch1+3], &cp[pitch+3].x );
    df = func-d->func;
    d->func = func;
    d->newpoint = true;
          /* test if the functional value decrement is sufficient */
    if ( !_g2bl_ComputeDeltaQd ( neqsw, qa22prof, QA22rows, gr2, dcoeff, &dq ) )
      goto failure;
    rk = df/dq;
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
    d->nG ++;
#endif
    g2bl_UFuncGradd ( nkn, qcoeff, Nitab,
                      lastknotu, lastknotv, pitch1, acp, NULL,
                      C, ftab, gtab, &func, grad );
    memcpy ( g, grad, nvars*sizeof(double) );
    pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, g );
    gn = sqrt ( pkn_ScalarProductd ( neqsw, gr2, gr2 ) );
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
        pkn_NRBLowerTrSolved ( neqsw, qa22prof, LM22, LM22rows,
                               1, 1, gr2, 1, dcoeff );
        pkn_NRBUpperTrSolved ( neqsw, qa22prof, LM22, LM22rows,
                               1, 1, dcoeff, 1, dcoeff );
        ldco = sqrt ( pkn_ScalarProductd ( neqsw, dcoeff, dcoeff ) );
        pkv_Selectd ( lastknotu-9, pitch2, pitch, pitch2, &cp[pitch+3].x, coeff );
        pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
        memcpy ( coeff, y, nconstr*sizeof(double) );
        pkn_SubtractMatrixd ( 1, neqsw, 0, coeff2, 0, dcoeff, 0, coeff2 );
        pkn_multiInvReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, coeff );
        pkv_Selectd ( lastknotu-9, pitch2, pitch2, pitch1, coeff, &acp[pitch1+3].x );
#ifdef _DEBUG
printf ( "G" );
#endif
#ifdef COUNT
        d->nG ++;
#endif
        g2bl_UFuncGradd ( nkn, qcoeff, Nitab, lastknotu, lastknotv,
                          pitch1, acp, NULL, C, ftab, gtab,
                          &func, grad );
        memcpy ( g, grad, nvars*sizeof(double) );
        pkn_multiReflectVectord ( nvars, nconstr, cmat, acmat, 1, 1, g );
        gn = sqrt ( pkn_ScalarProductd ( neqsw, gr2, gr2 ) );
        if ( func < d->func ) {
          d->func = func;
          d->newpoint = true;
          pkv_Selectd ( lastknotu-9, 3*(lastknotv-9), 3*(lastknotv-3),
                        pitch, &acp[3*(lastknotv-3)+3], &cp[pitch+3].x );
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
    if ( !pkn_NRBSymFindEigenvalueIntervald ( neqsw, qa22prof,
                             QA22rows[0], QA22rows, &gc, &gd ) )
      goto failure;
    g1 = 0.01*(gd-gc);
  }
  else
    g1 = d->nu;

  g0 = 0.0;
  for ( i = 0; ; i++ ) {
    fg1 = _g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr,
                   cmat, acmat, qsize, qa22prof, QA22rows[0],
                   LM22, LM22rows, g1, lastknotu, lastknotv, pitch, cp, acp,
                   y, gr2, coeff, dcoeff, C, ftab );
    RECORD_MIN ( g1, fg1 );
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
    fg1 = _g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr,
                   cmat, acmat, qsize, qa22prof, QA22rows[0],
                   LM22, LM22rows, g1, lastknotu, lastknotv, pitch, cp, acp,
                   y, gr2, coeff, dcoeff, C, ftab );
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
      fgb = _g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr,
                     cmat, acmat, qsize, qa22prof, QA22rows[0],
                     LM22, LM22rows, gb, lastknotu, lastknotv, pitch, cp, acp,
                     y, gr2, coeff, dcoeff, C, ftab );
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
      fgc = _g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr,
                     cmat, acmat, qsize, qa22prof, QA22rows[0],
                     LM22, LM22rows, gc, lastknotu, lastknotv, pitch, cp, acp,
                     y, gr2, coeff, dcoeff, C, ftab );
      RECORD_MIN ( gc, fgc );
    }
    else {
      ga = gc;  fga = fgc;
      gc = gd;  fgc = fgd;
      gd = GRDIV ( gb, gc );
      fgd = _g2bl_AuxConstrNuFuncd ( nkn, qcoeff, Nitab, nvars, nconstr,
                     cmat, acmat, qsize, qa22prof, QA22rows[0],
                     LM22, LM22rows, gd, lastknotu, lastknotv, pitch, cp, acp,
                     y, gr2, coeff, dcoeff, C, ftab );
      RECORD_MIN ( gd, fgd );
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
  if ( d->newpoint )
    pkv_Selectd ( lastknotu-9, pitch2, pitch2, pitch, mcoeff, &cp[pitch+3] );

next_iter:
  d->itn ++;
#ifdef _DEBUG
printf ( "\n" );
#endif
  if ( QA11rows[0] ) PKV_FREE ( QA11rows[0] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( QA11rows[0] ) PKV_FREE ( QA11rows[0] );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_IterBlSurfaceConstrOptLMTd*/

boolean g2bl_FindBlSurfaceConstrLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       int nconstr, double *constrmat,
                                       double *constrrhs,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 )
{
  void              *sp;
  lmt_constroptdata *data;
  int               itn;
  boolean           finished;

  sp = pkv_GetScratchMemTop ();

  data = NULL;
  if ( !g2bl_InitBlSurfaceConstrOptLMTd ( lastknotu, lastknotv, pitch, cp,
                                    nconstr, constrmat, constrrhs,
                                    C, dO, dM, nkn1, nkn2, (void**)&data ) )
    goto failure;

  finished = false;
  for ( itn = 0; itn < maxit; itn++ ) {
    if ( !g2bl_IterBlSurfaceConstrOptLMTd ( data, &finished ) )
      goto failure;
    if ( finished )
      break;
  }

#ifdef COUNT
printf ( "nF = %d, nG = %d, nH = %d, nN = %d, nLM = %d\n",
         data->nF, data->nG, data->nH, data->nN, data->nLM );
#endif
  g2bl_ConstrOptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2bl_ConstrOptLMTDeallocated ( (void**)&data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_FindBlSurfaceConstrLMTd*/

