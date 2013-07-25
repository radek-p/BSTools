
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* computing normal vector Bezier patches      */

char mbs_BezP3Normalf ( int degreeu, int degreev, const point3f *ctlpoints,
                        int *ndegu, int *ndegv, vector3f *ncp )
{
  vector3f *du, *dv;
  int      nn, mm, nm, n1, m1;
  int      i, j, k, l, in, iu, iv;
  void     *ScratchPtr;

  if ( degreeu < 1 || degreev < 1 ) {
    *ndegu = *ndegv = 0;
    SetVector3f ( ncp, 0.0, 0.0, 0.0 );
    return 0;
  }

  nm = degreeu*degreev;
  mbs_BezP3NormalDeg ( degreeu, degreev, ndegu, ndegv );
  nn = *ndegu;     mm = *ndegv;
  n1 = degreeu+1;  m1 = degreev+1;
  ScratchPtr = pkv_GetScratchMemTop ();
  du = pkv_GetScratchMem ( degreeu*(degreev+1)*sizeof(vector3f) );
  dv = pkv_GetScratchMem ( (degreeu+1)*degreev*sizeof(vector3f) );
  if ( !du || !dv ) {
    pkv_SetScratchMemTop ( ScratchPtr );
    return 0;
  }

  pkn_SubtractMatrixf ( 1, 3*degreeu*m1,
                        0, (float*)&ctlpoints[m1], 0, (float*)ctlpoints,
                        0, (float*)du );
  for ( i = 0; i < 3*degreeu*m1; i++ )
    ((float*)du)[i] *= nm;
  mbs_multiBezScalef ( degreeu-1, 1, 1, 3*m1, 0, (float*)du );
  mbs_multiBezScalef ( degreev, degreeu, 1, 3, 0, (float*)du );
  pkn_SubtractMatrixf ( n1, 3*degreev,
                        3*m1, (float*)&ctlpoints[1], 3*m1, (float*)ctlpoints,
                        3*degreev, (float*)dv );
  mbs_multiBezScalef ( degreeu, 1, 1, 3*degreev, 0, (float*)dv );
  mbs_multiBezScalef ( degreev-1, n1, 1, 3, 0, (float*)dv );

  memset ( ncp, 0, (nn+1)*(mm+1)*sizeof(vector3f) );
  for ( i = 0; i < degreeu; i++ )
    for ( j = 0; j <= degreev; j++ ) {
      iu = m1*i+j;
      for ( k = 0; k <= degreeu; k++ )
        for ( l = 0; l < degreev; l++ ) {
          iv = degreev*k+l;
          in = (mm+1)*(i+k)+j+l;
          ncp[in].x += du[iu].y*dv[iv].z - du[iu].z*dv[iv].y;
          ncp[in].y += du[iu].z*dv[iv].x - du[iu].x*dv[iv].z;
          ncp[in].z += du[iu].x*dv[iv].y - du[iu].y*dv[iv].x;
        }
    }

  mbs_multiBezUnscalef ( nn, 1, 1, 3*(mm+1), 0, (float*)ncp );
  mbs_multiBezUnscalef ( mm, nn+1, 1, 3, 0, (float*)ncp );

  pkv_SetScratchMemTop ( ScratchPtr );
  return 1;
} /*mbs_BezP3Normalf*/


/* /////////////////////////////////////////// */

typedef struct {
    double d12, d13, d14, d23, d24, d34;
  } outerp4x4d;

static void _OuterProdfd ( const point4f *u, const point4f *v, outerp4x4d *d )
{
  d->d12 += u->x*v->y - u->y*v->x;
  d->d13 += u->x*v->z - u->z*v->x;
  d->d14 += u->x*v->w - u->w*v->x;
  d->d23 += u->y*v->z - u->z*v->y;
  d->d24 += u->y*v->w - u->w*v->y;
  d->d34 += u->z*v->w - u->w*v->z;
} /*_OuterProdfd*/

static void _OuterProdIfd ( int i, const point4f *u, const point4f *v,
                            outerp4x4d *d )
{
  point4f r;

  switch ( i ) {
case  0:
    break;

case  1:
    _OuterProdfd ( u, v, d );
    break;

case -1:
    _OuterProdfd ( v, u, d );
    break;

default: 
    MultVector4f ( (double)i, u, &r );
    _OuterProdfd ( &r, v, d );
    break;
  }
} /*_OuterProdIfd*/

char mbs_BezP3RNormalf ( int degreeu, int degreev, const point4f *ctlpoints,
                         int *ndegu, int *ndegv, vector3f *ncp )
{
  int        m1, n1, m2, n2, mm, nn;
  int        i, j, k, l, im, ni, mj, km, lk, ll, niml, nkmj, ikm, mi, mk, ik;
  void       *ScratchPtr;
  vector4f   *r, *c, *v;
  outerp4x4d *d, *ad;
  vector3f   *nip;

  if ( degreeu <= 0 || degreev <= 0 ) {  /* degree overflow ignored for now */
    *ndegu = *ndegv = 0;
    SetPoint3f ( ncp, 0.0, 0.0, 0.0 );
    return 0;
  }

  mbs_BezP3RNormalDeg ( degreeu, degreev, ndegu, ndegv );
  nn = *ndegu;     mm = *ndegv;
  n1 = degreeu+1;  m1 = degreev+1;
  n2 = 2*degreeu;  m2 = 2*degreev;
  ScratchPtr = pkv_GetScratchMemTop ();
  r  = pkv_GetScratchMem ( n1*m1*sizeof(point4f) );
  d  = pkv_GetScratchMem ( n2*m2*sizeof(outerp4x4d) );
  c  = pkv_GetScratchMem ( n1*degreev*sizeof(point4f) );
  if ( !r || !d || !c ) {
    pkv_SetScratchMemTop ( ScratchPtr );
    return 0;
  }

  memcpy ( r, ctlpoints, n1*m1*sizeof(point4f) );
  mbs_multiBezScalef ( degreeu, 1, 1, 4*m1, 0, (float*)r );
  mbs_multiBezScalef ( degreev, n1, 1, 4, 0, (float*)r );
  pkn_SubtractMatrixf ( n1, 4*degreev,
                        4*m1, (float*)&ctlpoints[1], 4*m1, (float*)ctlpoints,
                        4*degreev, (float*)c );
  mbs_multiBezScalef ( degreeu, 1, 1, 4*degreev, 0, (float*)c );
  mbs_multiBezScalef ( degreev-1, n1, 1, 4, 0, (float*)c );

  memset ( d, 0, n2*m2*sizeof(outerp4x4d) );
  for ( i = im = 0;  i <= degreeu;  i++, im += m1 ) {
    for ( j = 0, ni = degreeu-i;  j <= degreev;  j++ ) {
      for ( k = km = 0, lk = min(degreeu,i+j), mj = degreev-j;
            k <= lk;
            k++, km += m1 ) {
        nkmj = (degreeu-k)*mj;
        ikm = m2*(i+k);
        for ( l = 0, ll = min(degreev,i+j-k-1); l <= ll; l++ ) {
          niml = ni*(degreev-l);
          _OuterProdIfd ( niml-nkmj, &r[im+j], &r[km+l], &d[ikm+j+l] );
        }
        l = i+j-k;
        if ( k < i && l <= degreev ) {
          niml = ni*(degreev-l);
          _OuterProdIfd ( niml-nkmj, &r[im+j], &r[km+l], &d[ikm+j+l] );
        }
      }
    }
  }

  memset ( ncp, 0, (nn+1)*(mm+1)*sizeof(vector3f) );
  for ( i = mi = 0;  i < n2;  i++, mi += m2 ) {
    for ( k = mk = 0;  k <= degreeu;  k++, mk += degreev ) {
      if ( i+k > 0 ) {
        for ( j = 0, ik = (mm+1)*(i+k-1);  j < m2;  j++ ) {
          for ( l = 0, ad = &d[mi+j];  l < degreev;  l++ ) {
            v   = &c[mk+l];
            nip = &ncp[ik+j+l];
            nip->x += (float)(- ad->d34*v->y + ad->d24*v->z - ad->d23*v->w);
            nip->y += (float)(+ ad->d34*v->x - ad->d14*v->z + ad->d13*v->w);
            nip->z += (float)(- ad->d24*v->x + ad->d14*v->y - ad->d12*v->w);
          }
        }
      }
    }
  }

  mbs_multiBezUnscalef ( nn, 1, 1, 3*(mm+1), 0, (float*)ncp );
  mbs_multiBezUnscalef ( mm, nn+1, 1, 3, 0, (float*)ncp );

  pkv_SetScratchMemTop ( ScratchPtr );
  return 1;
} /*mbs_BezP3RNormalf*/

