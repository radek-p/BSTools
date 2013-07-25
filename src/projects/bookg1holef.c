
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

#include "bookg1holef.h"


void (*G1OutCentralPointf) ( point3f *p ) = NULL;
void (*G1OutAuxCurvesf) ( int ncurves, int degree, const point3f *accp,
                          float t ) = NULL;
void (*G1OutStarCurvesf) ( int ncurves, int degree,
                           const point3f *sccp ) = NULL;
void (*G1OutAuxPatchesf) ( int npatches, int degu, int degv,
                           const point3f *apcp ) = NULL;


/* ////////////////////////////////////////////////////////////////////////// */
static float delta, csd, snd, u1, u2;

static void GetSurroundingPatches ( int hole_k,
                                    point3f*(*GetBezp)(int i, int j),
                                    point3f **npatch )
{
  int     i;
  point3f *spp, *spq;

  /* There are 2k patches surrounding the hole, each of degree (3,3). */
  /* We need only two rows, adjacent to the hole edge. */
  *npatch = (point3f*)pkv_GetScratchMem ( 2*hole_k*8*sizeof(point3f) );
  if ( !(*npatch) )
    exit ( 1 );
  for ( i = 0, spp = *npatch;  i < hole_k;  i++ ) {
    spq = GetBezp ( i, 1 );
    memcpy ( spp, spq, 8*sizeof(point3f) );
    spp += 8;
    spq = GetBezp ( i, 2 );
    memcpy ( spp, spq, 8*sizeof(point3f) );
    spp += 8;
  }
  if ( G1OutAuxPatchesf )
    G1OutAuxPatchesf ( 2*hole_k, 1, 3, *npatch );
} /*GetSurroundingPatches*/

static void Der1BC3 ( const point3f *bc, vector3f *der1 )
{
  SubtractPoints3f ( &bc[1], &bc[0], der1 );
  MultVector3f ( 3.0, der1, der1 );
} /*Der1BC3*/

static void Der2BC3 ( const point3f *bc, vector3f *der2 )
{
  AddVector3f ( &bc[0], &bc[2], der2 );
  AddVector3Mf ( der2, &bc[1], -2.0, der2 );
  MultVector3f ( 6.0, der2, der2 );
} /*Der2BC3*/

static float ConstructAuxCurves ( int hole_k, float beta1,
                                  const point3f *npatch,
                                  point3f **auxcurve )
{
  int      i, j, hk, nac;
  vector3f v0, v1, v2;
  float    lgt, auxct;
  point3f  *auxc;

  if ( hole_k & 0x1 ) {   /* k is odd */
    hk = (hole_k-1)/2;
    *auxcurve = (point3f*)pkv_GetScratchMem ( hole_k*4*sizeof(point3f) );
    if ( !(*auxcurve) )
      exit ( 1 );
    u1 = (float)((1.0+csd/3.0)*beta1);
    u2 = (float)((2.0*cos(0.5*delta))*beta1);
    lgt = u1+u2;
    for ( i = 0, auxc = *auxcurve;  i < hole_k;  i++, auxc += 4 ) {
      SubtractPoints3f ( &npatch[16*i], &npatch[16*i+4], &v0 );
      auxc[0] = npatch[16*i];
      AddVector3Mf ( &auxc[0], &v0, lgt, &auxc[1] );
      j = (i+hk+1) % hole_k;
      auxc[3] = npatch[16*j+3];
      SubtractPoints3f ( &auxc[3], &npatch[16*j+7], &v1 );
      j = (i+hk) % hole_k;
      SubtractPoints3f ( &auxc[3], &npatch[16*j+12], &v2 );
      AddVector3f ( &v1, &v2, &v0 );
      AddVector3Mf ( &auxc[3], &v0, lgt/u2, &auxc[2] );
    }
    nac = hole_k;
    auxct = u1/lgt;
  }
  else {                  /* k is even */
    hk = hole_k/2;
    *auxcurve = (point3f*)pkv_GetScratchMem ( hk*4*sizeof(point3f) );
    if ( !(*auxcurve) )
      exit ( 1 );
    u1 = u2 = (float)((1.0+csd/3.0)*beta1);
    lgt = u1+u2;
    for ( i = 0, auxc = *auxcurve; i < hk; i++, auxc += 4 ) {
      SubtractPoints3f ( &npatch[16*i], &npatch[16*i+4], &v0 );
      auxc[0] = npatch[16*i];
      AddVector3Mf ( &auxc[0], &v0, lgt, &auxc[1] );
      j = (i+hk) % hole_k;
      SubtractPoints3f ( &npatch[16*j], &npatch[16*j+4], &v0 );
      auxc[3] = npatch[16*j];
      AddVector3Mf ( &auxc[3], &v0, lgt, &auxc[2] );
    }
    nac = hk;
    auxct = 0.5;
  }
  if ( G1OutAuxCurvesf )
    G1OutAuxCurvesf ( nac, 3, *auxcurve, auxct );

  return auxct;
} /*ConstructAuxCurves*/

static void EvenCurvatureCorrection ( int hole_k, vector3f *normal_vector,
                                      point3f *starcurve )
{
#define EPS 1.0e-5
  void     *sp;
  int      i;
  vector3f w, *v;
  float    t, nw, snw, snv, *nv;
  point3f  *starc;

  sp = pkv_GetScratchMemTop ();
  v  = (vector3f*)pkv_GetScratchMem ( hole_k*sizeof(vector3f) );
  nv = pkv_GetScratchMemf ( hole_k );
  if ( !v || !nv )
    exit ( 1 );

  snw = snv = 0.0;
  for ( i = 0, starc = starcurve;  i < hole_k;  i++, starc += 4 ) {
    SubtractPoints3f ( &starc[1], &starc[2], &w );
    nw = (float)DotProduct3f ( normal_vector, &w );
    if ( i & 0x1 ) snw -= nw;
      else         snw += nw;
    SubtractPoints3f ( &starc[2], &starc[3], &v[i] );
    t = nv[i] = (float)DotProduct3f ( normal_vector, &v[i] );
    snv += t*t;
  }
  if ( fabs(snv) < EPS )
    return;
  snv = snw/snv;
  for ( i = 0, starc = starcurve;  i < hole_k;  i++, starc += 4 ) {
    if ( i & 0x1 ) t = -nv[i]*snv;
      else         t = +nv[i]*snv;
    AddVector3Mf ( &starc[2], &v[i], t, &starc[2] );
  }
  pkv_SetScratchMemTop ( sp );
#undef EPS
} /*EvenCurvatureCorrection*/

static void ConstructStarCurves ( int hole_k, point3f *auxcurve,
                                  float auxct, float beta2,
                                  point3f *central_point,
                                  vector3f *normal_vector,
                                  point3f **starcurve, vector3f **starpcd )
{
  void     *sp;
  int      i, hk, nac;
  trans2f  t;
  vector2f v;
  vector3f x, y, w;
  point3f  *auxc, *starc;
  vector3f *wvect, *starp;

        /* memory allocation */
  if ( hole_k & 0x1 )
    nac = hole_k;
  else
    nac = hole_k/2;
  *starcurve = (point3f*)pkv_GetScratchMem ( 8*hole_k*sizeof(point3f) );
  sp = pkv_GetScratchMemTop ();
  wvect = pkv_GetScratchMem ( nac*sizeof(vector3f) );
  if ( !wvect || !(*starcurve) )
    exit ( 1 );
  *starpcd = &(*starcurve)[4*hole_k];

        /* construction of the central point */
  SetVector3f ( &x, 0.0, 0.0, 0.0 );
  for ( i = 0, auxc = auxcurve;  i < nac;  i++, auxc += 4 ) {
    mbs_BCHornerDerC3f ( 3, auxc, auxct, &y, &wvect[i] );
    AddVector3f ( &x, &y, &x );
  }
  MultVector3f ( 1.0/(float)nac, &x, central_point );
  if ( G1OutCentralPointf )
    G1OutCentralPointf ( central_point );

        /* construction of star curves */
  SetVector3f ( &x, 0.0, 0.0, 0.0 );
  y = x;
  IdentTrans2f ( &t );
  RotTrans2f ( &t, (float)(2.0*PI/(float)hole_k) );
  if ( hole_k & 0x1 ) {  /* k is odd */
    SetVector2f ( &v, (float)((2.0/hole_k)*auxct/3.0*beta2), 0.0 );
    for ( i = 0; i < hole_k; i++ ) {
      AddVector3Mf ( &x, &wvect[i], v.x, &x );
      AddVector3Mf ( &y, &wvect[i], v.y, &y );
      TransVector2f ( &t, &v, &v );
    }
    CrossProduct3f ( &x, &y, normal_vector );
    SetVector2f ( &v, 1.0, 0.0 );
    for ( i = 0, auxc = auxcurve, starc = *starcurve, starp = *starpcd;
          i < hole_k;
          i++, auxc += 4, starc += 4, starp += 4 ) {
      starc[3] = auxc[0];
      InterPoint3f ( &auxc[0], &auxc[1], auxct, &starc[2] );
      starc[0] = *central_point;
      SetVector3f ( &w, -(v.x*x.x+v.y*y.x),
                        -(v.x*x.y+v.y*y.y),
                        -(v.x*x.z+v.y*y.z) );
      AddVector3f ( central_point, &w, &starc[1] );
      SetVector3f ( &starp[0], snd*(-v.y*x.x+v.x*y.x)/u1,
                               snd*(-v.y*x.y+v.x*y.y)/u1,
                               snd*(-v.y*x.z+v.x*y.z)/u1 );
      TransVector2f ( &t, &v, &v );
    }
  }
  else {            /* k is even */
    hk = hole_k/2;
    SetVector2f ( &v, (float)((4.0/hole_k)*auxct/3.0*beta2), 0.0 );
    for ( i = 0; i < hk; i++ ) {
      AddVector3Mf ( &x, &wvect[i], v.x, &x );
      AddVector3Mf ( &y, &wvect[i], v.y, &y );
      TransVector2f ( &t, &v, &v );
    }
    CrossProduct3f ( &x, &y, normal_vector );
    SetVector2f ( &v, 1.0, 0.0 );
    for ( i = 0, auxc = auxcurve, starc = *starcurve, starp = *starpcd;
          i < hk;
          i++, auxc += 4, starc += 4, starp += 4 ) {
      starc[3] = auxc[0];
      starc[4*hk+3] = auxc[3];
      InterPoint3f ( &auxc[0], &auxc[1], auxct, &starc[2] );
      InterPoint3f ( &auxc[3], &auxc[2], auxct, &starc[4*hk+2] );
      starc[0] = starc[4*hk] = *central_point;
      SetVector3f ( &w, -(v.x*x.x+v.y*y.x),
                        -(v.x*x.y+v.y*y.y),
                        -(v.x*x.z+v.y*y.z) );
      AddVector3f ( central_point, &w, &starc[1] );
      SubtractPoints3f ( central_point, &w, &starc[4*hk+1] );
      SetVector3f ( &starp[0], snd*(-v.y*x.x+v.x*y.x)/u1,
                               snd*(-v.y*x.y+v.x*y.y)/u1,
                               snd*(-v.y*x.z+v.x*y.z)/u1 );
      SetVector3f ( &starp[4*hk], -starp[0].x, -starp[0].y, -starp[0].z );
      TransVector2f ( &t, &v, &v );
    }
    if ( hole_k != 4 )
      EvenCurvatureCorrection ( hole_k, normal_vector, *starcurve );
  }

        /* dispose the wvect array */
  pkv_SetScratchMemTop ( sp );
  if ( G1OutStarCurvesf )
    G1OutStarCurvesf ( hole_k, 3, *starcurve );
} /*ConstructStarCurves*/

static void ConstructAuxPatches ( int hole_k, point3f *npatch,
                                  point3f *starcurve, vector3f *starpcd,
                                  float **cval1 )
{
#define EPS 1.0e-5
  void     *sp;
  int      i, j;
  point3f  *auxpatch, *auxp, *np;
  vector3f *starp, *spcd, *cc, *xx, *yy;
  vector3f v1, v2, tw;
  float    m;

      /* memory allocation */
  *cval1 = pkv_GetScratchMemf ( hole_k );
  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (hole_k+1)*sizeof(vector3f) );
  if ( !cc || !cval1 )
    exit ( 1 );

      /* setting appropriate twists at the hole edge centres */
  for ( i = 0, starp = starcurve, spcd = starpcd, np = npatch;
        i < hole_k;
        i++, starp += 4, spcd += 4, np += 16 ) {
    SubtractPoints3f ( &np[1], &np[0], &spcd[3] );
    SubtractPoints3f ( &np[4], &np[0], &v1 );
    SubtractPoints3f ( &np[5], &np[1], &v2 );
    SubtractPoints3f ( &v1, &v2, &tw );
    SubtractPoints3f ( &starp[2], &starp[3], &v2 );
    m = (float)(fabs(v1.x)+fabs(v1.y)+fabs(v1.z));
    if ( m >= EPS ) m = (float)((fabs(v2.x)+fabs(v2.y)+fabs(v2.z))/m);
      else          m = 1.0;
    AddVector3Mf ( &spcd[3], &tw, -m, &spcd[2] );
    (*cval1)[i] = m;
  }
      /* solving compatibility conditions equations at the centre */
        /* computing second order derivatives */
  m = (float)(-csd/(3.0*u1));
  for ( i = 0, starp = starcurve;  i < hole_k;  i++, starp += 4 ) {
    Der2BC3 ( &starp[0], &cc[i] );
    MultVector3f ( m, &cc[i], &cc[i] );
  }
        /* computing the right hand side */
  cc[hole_k] = cc[0];
  for ( i = 0; i < hole_k; i++ )
    SubtractPoints3f ( &cc[i+1], &cc[i], &cc[i] );

  if ( hole_k & 0x1 ) {  /* k odd, a definite system to solve */
    for ( i = 0; i < hole_k-1; i++ ) {
      if ( i & 0x1 ) AddVector3f ( &cc[hole_k-1], &cc[i], &cc[hole_k-1] );
        else    SubtractPoints3f ( &cc[hole_k-1], &cc[i], &cc[hole_k-1] );
    }
    MultVector3f ( 0.5, &cc[hole_k-1], &cc[hole_k-1] );
    for ( i = hole_k-2; i >= 0; i-- )
      SubtractPoints3f ( &cc[i], &cc[i+1], &cc[i] );
  }
  else {                 /* k even, a dual least squares problem */
                           /* solved k times for symmetry */
    xx = pkv_GetScratchMem ( 2*hole_k*sizeof(vector3f) );
    if ( !xx )
      exit ( 1 );
    yy = &xx[hole_k];

    memset ( xx, 0, (hole_k+1)*sizeof(vector3f) );
    for ( j = 0; j < hole_k; j++ ) {

          /* renumbering the equations and unknowns */
      for ( i = 0; i < hole_k; i++ )
        yy[i] = cc[(i+j) % hole_k];

          /* now solve it */
      for ( i = 0; i < hole_k-2; i++ )
        AddVector3Mf ( &yy[i+1], &yy[i], -(float)(i+1)/(float)(i+2),
                       &yy[i+1] );
      MultVector3f ( (float)(hole_k-1)/(float)hole_k, &yy[hole_k-2],
                     &yy[hole_k-2] );
      for ( i = hole_k-3; i >= 0; i-- ) {
        SubtractPoints3f ( &yy[i], &yy[i+1], &yy[i] );
        MultVector3f ( (float)(i+1)/(float)(i+2), &yy[i], &yy[i] );
      }
      yy[hole_k-1] = yy[hole_k-2];
      for ( i = hole_k-2; i > 0; i-- )
        AddVector3f ( &yy[i-1], &yy[i], &yy[i] );

          /* sum to obtain the average solution */
      for ( i = 0; i < hole_k; i++ )
        AddVector3f ( &xx[i], &yy[(i+hole_k-j) % hole_k], &xx[i] );
    }
    m = (float)(1.0/(float)hole_k);
    for ( i = 0; i < hole_k; i++ )
      MultVector3f ( m, &xx[i], &cc[i] );
  }

      /* compute mixed partial derivatives at the centre */
  for ( i = 0, spcd = starpcd;
        i < hole_k;
        i++, spcd += 4 )
    AddVector3Mf ( &spcd[0], &cc[i], 1.0/3.0, &spcd[1] );

  if ( G1OutAuxPatchesf ) {
    pkv_SetScratchMemTop ( sp );
    auxpatch = (point3f*)pkv_GetScratchMem ( hole_k*8*sizeof(point3f) );
    if ( !auxpatch )
      exit ( 1 );
    for ( i = 0, auxp = auxpatch, spcd = starpcd;
          i < hole_k;
          i++, auxp += 8, spcd += 4 ) {
      memcpy ( &auxp[0], &starcurve[4*i], 4*sizeof(point3f) );
      for ( j = 0; j <= 3; j++ )
        AddVector3f ( &auxp[j], &spcd[j], &auxp[j+4] );
    }
    G1OutAuxPatchesf ( hole_k, 1, 3, auxpatch );
  }
  pkv_SetScratchMemTop ( sp );
#undef EPS
} /*ConstructAuxPatches*/

static void Solve3x2 ( float *v1, float *v2, float *v3 )
{
/* solve a consistent system of 3 linear equations with 2 unknowns */
  float a, b, c;

  a = (float)fabs(v1[0]);  b = (float)fabs(v1[1]);  c = (float)fabs(v1[2]);
  if ( b > a && b > c ) {
    a = v1[0];  v1[0] = v1[1];  v1[1] = a;
    a = v2[0];  v2[0] = v2[1];  v2[1] = a;
    a = v3[0];  v3[0] = v3[1];  v3[1] = a;
  }
  else if ( c > a ) {
    a = v1[0];  v1[0] = v1[2];  v1[2] = a;
    a = v2[0];  v2[0] = v2[2];  v2[2] = a;
    a = v3[0];  v3[0] = v3[2];  v3[2] = a;
  }

  a = v1[1]/v1[0];  v2[1] -= a*v2[0];  v3[1] -= a*v3[0];
  a = v1[2]/v1[0];  v2[2] -= a*v2[0];  v3[2] -= a*v3[0];

  a = (float)fabs(v2[1]);  b = (float)fabs(v2[2]);
  if ( b > a ) {
    a = v2[1];  v2[1] = v2[2];  v2[2] = a;
    a = v3[1];  v3[1] = v3[2];  v3[2] = a;
  }
  v3[1] /= v2[1];
  v3[0] = (v3[0]-v2[0]*v3[1])/v1[0];
} /*Solve3x2*/

static void CompFinalPoints ( const float *b, const float *c,
                              const vector3f *cc, const vector3f *cd,
                              const vector3f *qq,
                              point3f *pk0, point3f *pk1, int pitch )
{
  float   b2i[3] = {1.0, 2.0, 1.0};
  float   b5i[6] = {1.0, 5.0, 10.0, 10.0, 5.0, 1.0};
  int     i, j, k;
  point3f *pka, *pkb, *pkc;

        /* degree elevation of the boundary curve */
        /* all input data are represented in scaled bases */
  for ( k = 0, pka = pk0;  k <= 5;  k++, pka += pitch )
    memset ( pka, 0, sizeof(point3f) );
  for ( i = 0; i <= 2; i++ )
    for ( j = 0; j <= 3; j++ ) {
      pkc = &pk0[(i+j)*pitch];
      AddVector3Mf ( pkc, &cc[j], b2i[i], pkc );
    }

        /* multiplication of polynomials and curves */
  for ( k = 0, pka = pk0, pkb = pk1;
        k <= 5;
        k++, pka += pitch, pkb += pitch )
    *pkb = *pka;

  for ( i = 0; i <= 3; i++ )
    for ( j = 0; j <= 2; j++ ) {
      pkc = &pk1[(i+j)*pitch];
      AddVector3Mf ( pkc, &cd[j], 0.2*b[i], pkc );
      AddVector3Mf ( pkc, &qq[i], 0.2*c[j], pkc );
    }

        /* conversion from scaled to Bernstein basis */
  for ( k = 1, pka = pk0+pitch, pkb = pk1+pitch;
        k < 5;
        k++, pka += pitch, pkb += pitch ) {
    MultVector3f ( 1.0/b5i[k], pka, pka );
    MultVector3f ( 1.0/b5i[k], pkb, pkb );
  }
} /*CompFinalPoints*/

static void ConstructFinalPatches ( int hole_k, point3f *npatch,
                                    point3f *starcurve, vector3f *starpcd,
                                    float *cval1,
                                    point3f *finalpatch )
{
  void     *sp;
  int      i, im1, j, k;
  point3f  *stci, *stcim1, *spci, *spcim1, *npi, *npim1, *fp;
  vector3f *q0b, *q1b, *r0b, *r1b, *q0c, *q1c, *r0c, *r1c,
           *qc0, *qc1, *rc0, *rc1;
  float    *b0, *b1, *c0, *c1, *f0, *f1, *g0, *g1;
  float    dcc, dgg, df, db;
  vector3f v1, v2, v3, v4, v5;

    /* memory allocation */
  sp = pkv_GetScratchMemTop ();
  q0b = (point3f*)pkv_GetScratchMem ( 44*sizeof(point3f)+28*sizeof(float) );
  if ( !q0b )
    exit ( 1 );
  q1b = &q0b[3];  r0b = &q1b[3];  r1b = &r0b[3];
  q0c = &r1b[3];  q1c = &q0c[4];  r0c = &q1c[4];  r1c = &r0c[4];
  qc0 = &r1c[4];  qc1 = &qc0[4];  rc0 = &qc1[4];  rc1 = &rc0[4];
  b0 = (float*)(&rc1[4]);  f0 = &b0[4];  b1 = &f0[4];  f1 = &b1[4];
  c0 = &f1[4];  g0 = &c0[3];  c1 = &g0[3];  g1 = &c1[3];

    /* the construction */

      /* for k odd all final patches are constructed with the same */
      /* polynomials, which are constructed outsize the main loop. */
      /* for k even only some of the coefficients vary. */

  b0[0] =              csd;  b0[3] =   0.0;
  c0[0] =  (float)(3.0*u1);  c0[2] =   3.0;
  b1[0] =              0.0;  b1[3] =   0.0;
                             c1[2] =  -3.0;
  f0[0] =              csd;  f0[3] =   0.0;
  g0[0] = (float)(-3.0*u1);  g0[2] =  -3.0;
  f1[0] =              0.0;  f1[3] =   0.0;
                             g1[2] =  -3.0;

  c0[1] = c0[0]+c0[2];
  g0[1] = g0[0]+g0[2];

  if ( hole_k & 0x1 ) {

          /* c1 and g1 are the same only if k is odd */
    c1[0] = (float)(-3.0*u1);  g1[0] = (float)(-3.0*u1);
    c1[1] = c1[0]+c1[2];
    g1[1] = g1[0]+g1[2];

        /* the polynomials f and g are obtained from compatibility */
        /* conditions */
    dcc = (float)(c0[2]/c0[0]-1.0);   dgg = (float)(g0[2]/g0[0]-1.0);
    df = dcc+dgg*f0[0];                db = dcc*b0[0]+dgg;
    b0[1] = (float)(b0[0]+db/3.0);  f0[1] = (float)(f0[0]+df/3.0);

    dcc = (float)(1.0-c0[0]/c0[2]);   dgg = (float)(g1[2]/g1[0]-1.0);
    df = dcc+dgg*f1[0];                db = dcc*b0[3]+dgg;
    b0[2] = (float)(b0[3]-db/3.0);  f1[1] = (float)(f1[0]+df/3.0);

    dcc = (float)(c1[2]/c1[0]-1.0);   dgg = (float)(1.0-g0[0]/g0[2]);
    df = dcc+dgg*f0[3];                db = dcc*b1[0]+dgg;
    b1[1] = (float)(b1[0]+db/3.0);  f0[2] = (float)(f0[3]-df/3.0);

    dcc = (float)(1.0-c1[0]/c1[2]);   dgg = (float)(1.0-g1[0]/g1[2]);
    df = dcc+dgg*f1[3];                db = dcc*b1[3]+dgg;
    b1[2] = (float)(b1[3]-db/3.0);  f1[2] = (float)(f1[3]-df/3.0);

    mbs_multiBezScalef ( 3, 4, 1, 1, 0, b0 );
  }

      /* now construct the subsequent patches */
  for ( i = 0, im1 = hole_k-1,
        stci = starcurve, stcim1 = &starcurve[(hole_k-1)*4],
        spci = starpcd,   spcim1 = &starpcd[(hole_k-1)*4],
        npi = npatch, npim1 = &npatch[(hole_k-1)*16], fp = finalpatch;
        i < hole_k;
        im1 = i, i++,
        stcim1 = stci, stci += 4,
        spcim1 = spci, spci += 4,
        npim1 = npi, npi += 16, fp += 36 ) {

        /* get boundary curves and cross derivatives */
    for ( j = 0; j <= 3; j++ ) {
      qc0[j] = stci[j];
      rc0[j] = stcim1[j];
      qc1[j] = npim1[11-j];
      rc1[j] = npi[j];
      q0c[j] = spci[j];
      r0c[j] = spcim1[j];
      SubtractPoints3f ( &npi[4+j], &npi[j], &r1c[j] );
      SubtractPoints3f ( &npim1[15-j], &npim1[11-j], &q1c[j] );
    }
        /* compute the derivatives of boundary curves */
    for ( j = 0; j < 3; j++ ) {
      SubtractPoints3f ( &qc0[j+1], &qc0[j], &q0b[j] );
      MultVector3f ( 3.0, &q0b[j], &q0b[j] );
      SubtractPoints3f ( &rc0[j+1], &rc0[j], &r0b[j] );
      MultVector3f ( 3.0, &r0b[j], &r0b[j] );
      SubtractPoints3f ( &qc1[j+1], &qc1[j], &q1b[j] );
      MultVector3f ( 3.0, &q1b[j], &q1b[j] );
      SubtractPoints3f ( &rc1[j+1], &rc1[j], &r1b[j] );
      MultVector3f ( 3.0, &r1b[j], &r1b[j] );
    }

    if ( !(hole_k & 0x1) ) {
        /* compute the first and last polynomial coefficients */
      c1[0] = (float)(-3.0*cval1[im1]);
      g1[0] = (float)(-3.0*cval1[i]);

        /* now other coefficients; the polynomials c and g are of degree 1 */
      c1[1] = c1[0]+c1[2];
      g1[1] = g1[0]+g1[2];

        /* compute the mixed partial and second order derivatives */
        /* at the centre */
      Der1BC3 ( r0c, &v1 );
      Der1BC3 ( q0c, &v2 );
      Der2BC3 ( rc0, &v4 );
      Der2BC3 ( qc0, &v5 );
      MultVector3f ( -b0[0], &v5, &v3 );
      AddVector3Mf ( &v3, &v2, -c0[0], &v3 );
      AddVector3Mf ( &v3, &v4, f0[0], &v3 );
      AddVector3Mf ( &v3, &v1, g0[0], &v3 );  /* now v3 is the v00 vector */
      AddVector3Mf ( &v3, &q0c[0], -(c0[2]-c0[0]), &v3 );
      AddVector3Mf ( &v3, &r0c[0], +(g0[2]-g0[0]), &v3 );
      Der1BC3 ( qc0, &v1 );
      Der1BC3 ( rc0, &v2 );
      Solve3x2 ( &v1.x, &v2.x, &v3.x );
      db = v3.x;  df = -v3.y;
      b0[1] = (float)(b0[0]+db/3.0);   f0[1] = (float)(f0[0]+df/3.0);

      dcc = (float)(1.0-c0[0]/c0[2]);  dgg = (float)(g1[2]/g1[0]-1.0);
      df = dcc+dgg*f1[0];              db = dcc*b0[3]+dgg;
      b0[2] = (float)(b0[3]-db/3.0);   f1[1] = (float)(f1[0]+df/3.0);

      dcc = (float)(c1[2]/c1[0]-1.0);  dgg = (float)(1.0-g0[0]/g0[2]);
      df = dcc+dgg*f0[3];              db = dcc*b1[0]+dgg;
      b1[1] = (float)(b1[0]+db/3.0);   f0[2] = (float)(f0[3]-df/3.0);

      dcc = (float)(1.0-c1[0]/c1[2]);  dgg = (float)(1.0-g1[0]/g1[2]);
      df = dcc+dgg*f1[3];              db = dcc*b1[3]+dgg;
      b1[2] = (float)(b1[3]-db/3.0);   f1[2] = (float)(f1[3]-df/3.0);

      mbs_multiBezScalef ( 3, 4, 1, 1, 0, b0 );
    }
        /* convert the curves b to scaled bases */
    mbs_multiBezScalef ( 2, 4, 1, 3, 0, (float*)q0b );
    mbs_multiBezScalef ( 3, 8, 1, 3, 0, (float*)q0c );

    memset ( fp, 0, 36*sizeof(point3f) );
    CompFinalPoints ( b0, c0, qc0, q0b, q0c, &fp[0], &fp[1], 6 );
    CompFinalPoints ( b1, c1, qc1, q1b, q1c, &fp[5], &fp[4], 6 );
    CompFinalPoints ( f0, g0, rc0, r0b, r0c, &fp[0], &fp[6], 1 );
    CompFinalPoints ( f1, g1, rc1, r1b, r1c, &fp[30], &fp[24], 1 );

        /* construct the innermost points */
    for ( j = 2; j <= 3; j++ ) {
      q0b[0] = fp[6*j];
      InterPoint3f ( &q0b[0], &fp[6*j+1], 5.0/3.0, &q0b[1] );
      q0b[3] = fp[6*j+5];
      InterPoint3f ( &q0b[3], &fp[6*j+4], 5.0/3.0, &q0b[2] );
      mbs_BCDegElevC3f ( 3, q0b, 2, &k, r0b );
      fp[6*j+2] = r0b[2];  fp[6*j+3] = r0b[3];
    }
    for ( j = 2; j <= 3; j++ ) {
      q0b[0] = fp[j];
      InterPoint3f ( &q0b[0], &fp[j+6], 5.0/3.0, &q0b[1] );
      q0b[3] = fp[j+30];
      InterPoint3f ( &q0b[3], &fp[j+24], 5.0/3.0, &q0b[2] );
      mbs_BCDegElevC3f ( 3, q0b, 2, &k, r0b );
      MidPoint3f ( &fp[j+12], &r0b[2], &fp[j+12] );
      MidPoint3f ( &fp[j+18], &r0b[3], &fp[j+18] );
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*ConstructFinalPatches*/

boolean FillG1Holef ( int hole_k, point3f*(*GetBezp)(int i, int j),
                      float beta1, float beta2,
                      point3f *hpcp )
{
  void     *sp;
  point3f  *npatch;
  point3f  *auxcurve, *starcurve;
  vector3f *starpcd;
  float    auxct;
  float    *cval1;
  point3f  central_point;
  vector3f normal_vector;

  sp = pkv_GetScratchMemTop ();

  delta = (float)(2.0*PI/(float)hole_k);
  csd = (float)cos(delta);
  snd = (float)sin(delta);

  GetSurroundingPatches ( hole_k, GetBezp, &npatch );
  auxct = ConstructAuxCurves ( hole_k, beta1, npatch, &auxcurve );
  ConstructStarCurves ( hole_k, auxcurve, auxct, beta2,
                        &central_point, &normal_vector, &starcurve, &starpcd );
  ConstructAuxPatches ( hole_k, npatch, starcurve, starpcd, &cval1 );
  ConstructFinalPatches ( hole_k, npatch, starcurve, starpcd, cval1, hpcp );

  pkv_SetScratchMemTop ( sp );
  return true;
} /*FillG1Holef*/

