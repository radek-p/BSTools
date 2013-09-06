
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* multiplication of curves by scalar spline functions */

/* The procedure below computes nproducts products of univariate scalar  */
/* polynomials and vector curves. The data and result are represented in */
/* scaled Bernstein bases. This procedure is static, as the application  */
/* interface should be of higher level.                                  */
static void _mbs_multiMultScaledBezf ( int nproducts, int spdimen,
                                       int degpoly, int polypitch,
                                       const float *polycoef,
                                       int degcurve, int curvepitch,
                                       const float *curvecoef,
                                       int outpitch, float *outcoef )
{
  int   i, j, k, l, outdeg;
  float p;

  outdeg = degpoly+degcurve;
  for ( k = 0;
        k < nproducts;
        k++,
        polycoef += polypitch, curvecoef += curvepitch, outcoef += outpitch ) {
    memset ( outcoef, 0, (outdeg+1)*spdimen*sizeof(float) );
    for ( i = 0; i <= degpoly; i++ ) {
      p = polycoef[i];
      for ( j = 0; j <= degcurve; j++ )
        for ( l = 0; l < spdimen; l++ )
          outcoef[spdimen*(i+j)+l] += p*curvecoef[spdimen*j+l];
    }
  }
} /*_mbs_multiMultScaledBezf*/

/* The two procedures below may be called either with nscf equal to nvecf */
/* or with one of the two parameters equal to one. The number of products    */
/* is equal to the greater parameter. The degree of the product is the sum   */
/* of the degrees of the arguments, so outpitch should not be less than      */
/* spdimen*(degscf+degvecf). The arguments are represented in Bernstein   */
/* polynomial bases of degrees specified by the parameters discussed above.  */
void mbs_multiMultBezCf ( int nscf, int degscf, int scfpitch,
                          const float *scfcoef,
                          int spdimen,
                          int nvecf, int degvecf, int vecfpitch,
                          const float *vecfcp,
                          int *degprod, int prodpitch, float *prodcp )
{
  void  *stp;
  float *auxscf, *auxvecf;
  int   nprod;

  if ( nscf < nvecf ) { scfpitch = 0;   nscf = 1; }
  if ( nvecf > nscf ) { vecfpitch = 0;  nvecf = 1; }
  nprod = max ( nscf, nvecf );

/* Allocate working memory area and copy data. */
  stp = pkv_GetScratchMemTop();
  auxscf  = pkv_GetScratchMemf ( nscf*(degscf+1) );
  auxvecf = pkv_GetScratchMemf ( nvecf*(degvecf+1)*spdimen );
  if ( !auxscf || !auxvecf ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }
  pkv_Selectf ( nscf, degscf+1, scfpitch, degscf+1, scfcoef, auxscf );
  pkv_Selectf ( nvecf, (degvecf+1)*spdimen, vecfpitch,
                (degvecf+1)*spdimen, vecfcp, auxvecf );

/* Convert data to scaled Bernstein bases */
  mbs_multiBezScalef ( degscf, nscf, 1, 1, degscf+1, auxscf );
  mbs_multiBezScalef ( degvecf, nvecf, 1, spdimen, (degvecf+1)*spdimen,
                       auxvecf );

/* Do the actual multiplication and convert the result to Bernstein basis. */
  _mbs_multiMultScaledBezf ( nprod, spdimen, degscf, degscf+1, auxscf,
                             degvecf, (degvecf+1)*spdimen, auxvecf,
                             prodpitch, prodcp );
  *degprod = degscf+degvecf;
  mbs_multiBezUnscalef ( *degprod, 1, nprod, spdimen, prodpitch, prodcp );

  pkv_SetScratchMemTop ( stp );
} /*mbs_multiMultBezCf*/

/* The value of the prodlastknot parameter must be computed before calling  */
/* this procedure, by the instruction like                                  */
/*    prodlastknot = mbs_BSProdRepSizef ( degscf, scflastknot, scfknots,    */
/*                             degvecf, vecflastknot, vecfknots );          */
/* The caller must ensure that the length of the array pointed by prodknots */
/* is at least prodlastknot+1. The procedure below does not repeat this     */
/* computation.                                                             */
void mbs_multiMultBSCf ( int nscf, int degscf,
                         int scflastknot, const float *scfknots,
                         int scfpitch, const float *scfcoef,
                         int spdimen,
                         int nvecf, int degvecf,
                         int vecflastknot, const float *vecfknots,
                         int vecfpitch, const float *vecfcp,
                         int *degprod, int *prodlastknot, float *prodknots,
                         int prodpitch, float *prodcp )
{
  void  *stp;
  int   nprod, lastpknot, auxlastkn, npp, degpr, ascpitch, avecpitch;
  int   i, j;
  float *auxkn, *sauxc, *vauxc, *akn, *acp, *auxprod, *auxpr, *aakn;

  if ( nscf < nvecf ) { scfpitch = 0;   nscf = 1; }
  if ( nvecf > nscf ) { vecfpitch = 0;  nvecf = 1; }
  nprod = max ( nscf, nvecf );

  stp = pkv_GetScratchMemTop ();

/* The first step is to create the final knot sequence of the product */
/* representation. It is necessary to do it now, as it determines the */
/* amount of scratch memory needed for intermediate results.          */
  mbs_SetBSProdKnotsf ( degscf, scflastknot, scfknots,
                        degvecf, vecflastknot, vecfknots,
                        degprod, &lastpknot, prodknots );
  if ( lastpknot > *prodlastknot ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    exit ( 1 );  /* apparently the array for the result knots is too short */
  }
  *prodlastknot = lastpknot;

  *degprod = degpr = degscf+degvecf;
  npp = mbs_NumKnotIntervalsf ( degpr, lastpknot, prodknots );
  auxkn = pkv_GetScratchMemf ( (npp+1)*(degpr+1) );
  aakn = pkv_GetScratchMemf ( degpr+1 );
  if ( !auxkn || !aakn ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }

/* Find piecewise Bezier representations of the arguments, with use of the  */
/* Oslo coefficients matrices. Unfortunately, we have to copy data, as they */
/* may be represented by knots of multiplicity greater than degree+1.       */
  akn = pkv_GetScratchMemf ( max(scflastknot,vecflastknot)+1 );
  i = (scflastknot-degscf)*nscf;
  j = (vecflastknot-degvecf)*nvecf*spdimen;
  acp = pkv_GetScratchMemf ( max(i,j) );
  if ( !akn || !acp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }

  memcpy ( akn, scfknots, (scflastknot+1)*sizeof(float) );
  pkv_Selectf ( nscf, scflastknot-degscf, scfpitch, scflastknot-degscf,
                scfcoef, acp );
  if ( nscf > 1 )
    scfpitch = scflastknot-degscf;
  else
    scfpitch = 0;
  mbs_multiRemoveSuperfluousKnotsf ( nscf, 1, degscf, &scflastknot,
                                     akn, scfpitch, scfpitch, acp );
  if ( akn[1] != akn[degscf] ) {
    for ( i = 0; i <= degscf; i++ )
      aakn[i] = akn[degscf];
    mbs_multiBSChangeLeftKnotsf ( nscf, 1, degscf, akn, scfpitch, acp, aakn );
  }
  if ( akn[scflastknot-degscf] != akn[scflastknot-1] ) {
    for ( i = 0; i <= degscf; i++ )
      aakn[i] = akn[scflastknot-degscf];
    mbs_multiBSChangeRightKnotsf ( nscf, 1, degscf, scflastknot, akn,
                                   scfpitch, acp, aakn );
  }
  mbs_SetKnotPatternf ( lastpknot-2, &prodknots[1], degscf+1,
                        &auxlastkn, auxkn );
  auxkn[0] = akn[0];
  auxkn[auxlastkn] = akn[scflastknot];
  ascpitch = auxlastkn-degscf;
  sauxc = pkv_GetScratchMemf ( ascpitch*nscf );
  if ( !sauxc ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }
  mbs_multiOsloInsertKnotsf ( nscf, 1, degscf, scflastknot, akn,
                              scfpitch, acp, auxlastkn, auxkn,
                              ascpitch, sauxc );

/* After that we do not need knots of this representation any more.  */
/* Now do the same thing with the vector argument of multiplication. */
  memcpy ( akn, vecfknots, (vecflastknot+1)*sizeof(float) );
  pkv_Selectf ( nvecf, (vecflastknot-degvecf)*spdimen, vecfpitch,
                (vecflastknot-degvecf)*spdimen, vecfcp, acp );
  if ( nvecf > 1 )
    vecfpitch = (vecflastknot-degvecf)*spdimen;
  else
    vecfpitch = 0;
  mbs_multiRemoveSuperfluousKnotsf ( nvecf, spdimen, degvecf, &vecflastknot,
                                     akn, vecfpitch, vecfpitch, acp );
  if ( akn[1] != akn[degvecf] ) {
    for ( i = 0; i <= degvecf; i++ )
      aakn[i] = akn[degvecf];
    mbs_multiBSChangeLeftKnotsf ( nvecf, spdimen, degvecf, akn,
                                  vecfpitch, acp, aakn );
  }
  if ( akn[vecflastknot-degvecf] != akn[vecflastknot-1] ) {
    for ( i = 0; i <= degvecf; i++ )
      aakn[i] = akn[vecflastknot-degvecf];
    mbs_multiBSChangeRightKnotsf ( nvecf, spdimen, degvecf, vecflastknot, akn,
                                   vecfpitch, acp, aakn );
  }
  mbs_SetKnotPatternf ( lastpknot-2, &prodknots[1], degvecf+1,
                        &auxlastkn, auxkn );
  auxkn[0] = akn[0];
  auxkn[auxlastkn] = akn[vecflastknot];
  avecpitch = (auxlastkn-degvecf)*spdimen;
  vauxc = pkv_GetScratchMemf ( avecpitch*nvecf );
  if ( !vauxc ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }
  mbs_multiOsloInsertKnotsf ( nvecf, spdimen, degvecf, vecflastknot, akn,
                              vecfpitch, acp, auxlastkn, auxkn,
                              avecpitch, vauxc );

/* Now the actual multiplication, preceded by conversion to scaled Bernstein */
/* bases and followed by conversion of the product to Bernstein basis.       */
  mbs_multiBezScalef ( degscf, npp, nscf, 1, scfpitch, sauxc );
  mbs_multiBezScalef ( degvecf, npp, nvecf, spdimen, vecfpitch, vauxc );
  auxprod = pkv_GetScratchMemf ( nprod*npp*(degpr+1)*spdimen );
  if ( !auxprod ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }
  for ( i = 0, auxpr = auxprod;
        i < nprod;
        i++, sauxc += scfpitch, vauxc += vecfpitch,
        auxpr += npp*(degpr+1)*spdimen )
    _mbs_multiMultScaledBezf ( npp, spdimen, degscf, degscf+1, sauxc,
                               degvecf, (degvecf+1)*spdimen, vauxc,
                               (degpr+1)*spdimen, auxpr );
  mbs_multiBezUnscalef ( degpr, npp, nprod, spdimen,
                         (degpr+1)*spdimen, auxprod );

/* The last thing to do is to generate the knot sequence of this piecewise */
/* Bezier representation of the product and then knot removal, with use    */
/* of the Oslo algorithm. */
  mbs_SetKnotPatternf ( lastpknot-2, &prodknots[1], degpr+1,
                        &auxlastkn, auxkn );
  auxkn[0] = prodknots[0];
  auxkn[auxlastkn] = prodknots[lastpknot];
  mbs_multiOsloRemoveKnotsLSQf ( nprod, spdimen, degpr, auxlastkn, auxkn,
                                 npp*(degpr+1)*spdimen, auxprod,
                                 lastpknot, prodknots,
                                 prodpitch, prodcp );

  pkv_SetScratchMemTop ( stp );
} /*mbs_multiMultBSCf*/

