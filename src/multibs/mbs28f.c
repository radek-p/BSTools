
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
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

/* ////////////////////////////////////////////////////////////// */
/* Computing intersections of lines with a trimmed patch boundary */

/* The boundary processed by the procedure mbs_FindBoundLineIntersectionsf */
/* must be described by a code produced by mbs_CompileTrimPatchBoundf.     */
 
static float _get_a3 ( const point2f *p0, const vector2f *normal,
                       const point3f *cp )
{
  vector2f v;
  float    w;

  w = cp->z;
  SetVector2f ( &v, cp->x - w*p0->x, cp->y - w*p0->y );
  return (float)DotProduct2f ( normal, &v );
} /*_get_a3*/

static float _get_a2 ( const point2f *p0, const vector2f *normal,
                       const point2f *cp )
{
  point3f q;

  SetPoint3f ( &q, cp->x, cp->y, 1.0 );
  return _get_a3 ( p0, normal, &q );
} /*_get_a2*/

static float _get_lparamf ( const point2f *p0, float t0,
                            const point2f *p1, float t1,
                            const point2f *q )
{
  vector2f a, b;
  float    v;

  SubtractPoints2f ( p1, p0, &a );
  SubtractPoints2f (  q, p0, &b );
  v = (float)(DotProduct2f ( &a, &b )/DotProduct2f ( &a, &a ));
  return t0 + v*(t1-t0);
} /*_get_lparamf*/

static void _mbs_AddLineLineIntersf ( const point2f *p0, float t0,
                                      const point2f *p1, float t1,
                                      mbs_signpoint1f *inters,
                                      char dim, const float *vert,
                                      float u, char sa )
{
  point2f p;
  point3f q;

  if ( dim == 2 )
    InterPoint2f ( (point2f*)vert, (point2f*)&vert[2], u, &p );
  else {
    InterPoint3f ( (point3f*)vert, (point3f*)&vert[3], u, &q );
    Point3to2f ( &q, &p );
  }
  inters->t = _get_lparamf ( p0, t0, p1, t1, &p );
  inters->sign1 = sa;
} /*_mbs_AddLineLineIntersf*/

static float _mbs_FindLinePolylineIntersf ( const point2f *p0, float t0,
                                            const point2f *p1, float t1,
                                            const vector2f *normal,
                                            char dim, short nlines,
                                            const float *vert,
                                            float a,
                                            mbs_signpoint1f *inters, int *ninters )
{
  int   i, ni, maxinters;
  float b;
  char  sa, sb;

  maxinters = *ninters;
  *ninters = ni = 0;
  sa = pkv_signf ( a );
  for ( i = 1; i <= nlines; i++ ) {
    if ( dim == 2 )
      b = _get_a2 ( p0, normal, (point2f*)(&vert[2*i]) );
    else
      b = _get_a3 ( p0, normal, (point3f*)(&vert[3*i]) );
    sb = pkv_signf ( b );
    if ( sb != sa ) {   /* there is perhaps an intersection */
      if ( ni >= maxinters )
        goto out;
      _mbs_AddLineLineIntersf ( p0, t0, p1, t1, &inters[ni],
                                dim, &vert[dim*(i-1)], -a/(b-a), sa );
      ni++;
    }
    a = b;
    sa = sb;
  }
  *ninters = ni;
  return a;

out:
  *ninters = -1;
  return 0.0;
} /*_mbs_FindLinePolylineIntersf*/

typedef struct {
    int   deg;  /* polynomial degree, for function _mbs_apoly et al. */
    float *ac;  /* polynomial coefficients */
    int   sp;   /* stack pointer, 0 if stack empty */
    float *sd;  /* stack data pointer */
  } _mbs_apolystrf;

static float _mbs_apoly ( void *usrptr, float t )
{
  float value;

  mbs_BCHornerC1f ( ((_mbs_apolystrf*)usrptr)->deg,
                    ((_mbs_apolystrf*)usrptr)->ac, t, &value );
  return value;
} /*_mbs_apoly*/

static boolean _mbs_varsign ( int deg, float *ac )
{
  int  i;
  char s, z;

  s = pkv_signf ( ac[0] );
  for ( i = 1; i <= deg; i++ ) {
    z = pkv_signf ( ac[i] );
    if ( z != s )
      return true;
  }
  return false;
} /*_mbs_varsign*/

static boolean _mbs_onesignch ( int deg, float *ac, boolean *nonzero )
{
  int  i, k;
  char sa, sb;

  k = 0;
  sa = pkv_signf ( ac[0] );
  *nonzero = (boolean)(sa != 0);
  for ( i = 1; i <= deg; i++ ) {
    sb = pkv_signf ( ac[i] );
    *nonzero = (boolean)(*nonzero || sb != 0);
    if ( sb != sa ) {
      sa = sb;
      k++;
      if ( k > 1 )
        return false;
    }
  }
  return true;
} /*_mbs_onesignch*/

static boolean _mbs_pushp ( _mbs_apolystrf *ap, float t0, float t1, float *ac )
{
  int   deg;
  float *sd;

  deg = ap->deg;
  sd = ap->sd = pkv_GetScratchMemf ( deg+3 );
  if ( !sd ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
    return false;
  }
  sd[0] = t0;
  sd[1] = t1;
  memcpy ( &sd[2], ac, (deg+1)*sizeof(float) );
  ap->sp ++;
  return true;
} /*_mbs_pushp*/

static boolean _mbs_popp ( _mbs_apolystrf *ap, float *t0, float *t1, float *ac )
{
  int deg;

  if ( !ap->sp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
    return false;
  }
  deg = ap->deg;
  ap->sp--;
  *t0 = ap->sd[0];
  *t1 = ap->sd[1];
  memcpy ( ac, &ap->sd[2], (deg+1)*sizeof(float) );
  ap->sd -= deg+3;
  pkv_FreeScratchMemf ( deg+3 );
  return true;
} /*_mbs_popp*/

static void _mbs_AddLineBezcIntersf ( const point2f *p0, float t0,
                                      const point2f *p1, float t1,
                                      mbs_signpoint1f *inters,
                                      char dim, int deg, const float *cpts,
                                      float u, char as )
{
  point2f q;

  if ( dim == 2 )
    mbs_BCHornerC2f ( deg, (point2f*)cpts, u, &q );
  else
    mbs_BCHornerC2Rf ( deg, (point3f*)cpts, u, &q );
  inters->t = _get_lparamf ( p0, t0, p1, t1, &q );
  inters->sign1 = as;
} /*_mbs_AddLineBezcIntersf*/

static float _mbs_FindLineBezcIntersf ( const point2f *p0, float t0,
                                        const point2f *p1, float t1,
                                        const vector2f *normal,
                                        char dim, short deg,
                                        const float *cpts,
                                        float a,
                                        mbs_signpoint1f *inters, int *ninters )
{
#define eps 1.0e-6
  void           *ssp;
  int            i, ni, maxinters;
  float          *ac, *bc;
  float          u0, u1, u;
  _mbs_apolystrf ap;
  boolean        nonzero, error;

  ssp = pkv_GetScratchMemTop ();

                        /* compute the polynomial coefficients */
  ac = pkv_GetScratchMemf ( 3*(deg+1) );
  if ( !ac )
    goto out;
  ap.ac = &ac[deg+1];
  bc  = &ap.ac[deg+1];

  ac[0] = a;
  if ( dim == 2 ) {
    for ( i = 1; i <= deg; i++ )
      ac[i] = _get_a2 ( p0, normal, (point2f*)&cpts[2*i] );
  }
  else { /* dim == 3 */
    for ( i = 1; i <= deg; i++ )
      ac[i] = _get_a3 ( p0, normal, (point3f*)&cpts[3*i]);
  }
  a = ac[deg];

                        /* find the polynomial zeros */
  maxinters = *ninters;
  *ninters = ni = 0;
  ap.deg = deg;
  ap.sp = 0;
  if ( _mbs_varsign ( deg, ac ) ) {
    if ( !_mbs_pushp ( &ap, 0.0, 1.0, ac ) )
      goto out;
    do {
      if ( !_mbs_popp ( &ap, &u0, &u1, ap.ac ) )
        goto out;
      if ( _mbs_onesignch ( deg, ap.ac, &nonzero ) ) {
        if ( nonzero ) {
          u = pkn_Illinoisf ( _mbs_apoly, (void*)&ap, 0.0, 1.0, eps, &error );
          u = u0 + u*(u1-u0);

          if ( ni >= maxinters )
            goto out;
          _mbs_AddLineBezcIntersf ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, u, pkv_signf(ap.ac[0]) );
          ni++;
        }
        else {
          if ( ni >= maxinters-1 )
            goto out;
          _mbs_AddLineBezcIntersf ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, u0, 0 );
          _mbs_AddLineBezcIntersf ( p0, t0, p1, t1, &inters[ni+1],
                                    dim, deg, cpts, u1, 0 );
          ni += 2;
        }
      }
      else {
        if ( u1-u0 > eps ) {
          mbs_BisectBC1f ( deg, ap.ac, bc );
          u = (float)(0.5*(u0+u1));
          if ( _mbs_varsign ( deg, ap.ac ) ) {
            if ( !_mbs_pushp ( &ap, u, u1, ap.ac ) )
              goto out;
          }
          if ( _mbs_varsign ( deg, bc ) ) {
            if ( !_mbs_pushp ( &ap, u0, u, bc ) )
              goto out;
          }
        }
        else {
          if ( ni >= maxinters )
            goto out;
          _mbs_AddLineBezcIntersf ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, (float)(0.5*(u0+u1)),
                                    pkv_signf(ap.ac[0]) );
          ni++;
        }
      }
    } while ( ap.sp );
  }

  pkv_SetScratchMemTop ( ssp );
  *ninters = ni;
  return a;

out:
  pkv_SetScratchMemTop ( ssp );
  *ninters = -1;
  return 0.0;
#undef eps
} /*_mbs_FindLineBezcIntersf*/

static boolean _mbs_GetBoundPointf ( const char *bufp, short k,
                                     const point2f *p0, const vector2f *normal,
                                     point3f *point, float *a )
{
  point2f *q;

/* it computes the number a, in order to avoid */
/* the nasty consequences of rounding errors.  */

  if ( bufp[1] == 2 ) {
    q = ((point2f*)&bufp[2+sizeof(short)+2*k*sizeof(float)]);
    SetPoint3f ( point, q->x, q->y, 1.0 );
  }
  else if ( bufp[1] == 3 )
    *point = *(point3f*)&bufp[2+sizeof(short)+3*k*sizeof(float)];
  else
    return false;
  *a = _get_a3 ( p0, normal, point );
  return true;
} /*_mbs_GetBoundPointf*/

static void _mbs_CleanLoopf ( mbs_signpoint1f *inters, int *ninters )
{
  int i, j, ni;

  if ( !(ni = *ninters) )
    return;
  for ( i = 0; i < ni-1; i++ )
    inters[i].sign2 = inters[i+1].sign1;
  inters[ni-1].sign2 = inters[0].sign1;

  for ( i = 1, j = 0; i < ni; ) {
    if ( inters[i].t == inters[j].t && !inters[j].sign2 && !inters[i].sign1 ) {
      inters[j].sign2 = inters[i++].sign2;
    }
    else
      inters[++j] = inters[i++];
  }
  if ( inters[0].t == inters[j].t && !inters[j].sign2 && !inters[0].sign1 ) {
    inters[0].sign1 = inters[j].sign1;
    ni = j;
  }
  else
    ni = j+1;
  for ( i = j = 0; i < ni; i++ )
    if ( inters[i].sign1 != inters[i].sign2 )
      inters[j++] = inters[i];

  *ninters = j;
} /*_mbs_CleanLoopf*/

void mbs_FindBoundLineIntersectionsf ( const void *bound,
                                       const point2f *p0, float t0,
                                       const point2f *p1, float t1,
                                       mbs_signpoint1f *inters, int *ninters )
{
  int      maxinters, ni, lni;
  short    N;
  vector2f normal;
  point3f  slp[3];
#define startp slp[0]
#define lastp  slp[1]
#define nextp  slp[2]
  char     *bufp, *bp, dim;
  float    *points, a, b = 0.0;

  SubtractPoints2f ( p1, p0, &normal );
  SetVector2f ( &normal, -normal.y, normal.x );

  maxinters = *ninters;
  *ninters = lni = 0;

  bufp = (char*)bound;
  if ( !_mbs_GetBoundPointf ( bufp, 0, p0, &normal, &startp, &a ) )
    goto out;

  while ( bufp[0] != 4 ) {
    dim    = bufp[1];
    N      = *((short*)(&bufp[2]));
    if ( dim < 2 || dim > 3 || N < 1 )
      goto out;
    points = ((float*)(&bufp[2+sizeof(short)]));

    ni = maxinters-*ninters;
    switch ( bufp[0] ) {
 case 0:
 case 1:
      a = _mbs_FindLinePolylineIntersf ( p0, t0, p1, t1, &normal, dim, N,
                                         points, a, &inters[*ninters], &ni );
      break;

 case 2:
 case 3:
      a = _mbs_FindLineBezcIntersf ( p0, t0, p1, t1, &normal, dim, N,
                                     points, a, &inters[*ninters], &ni );
      break;

 default:
      goto out;          /* signal invalid data */
    }
    if ( ni < 0 )
      goto out;
    *ninters += ni;

                         /* now deal with the joining line */
    if ( !_mbs_GetBoundPointf ( bufp, N, p0, &normal, &lastp, &a ) )
      goto out;
    bp = &bufp[2+sizeof(short)+(N+1)*dim*sizeof(float)];
    if ( bufp[0] == 1 || bufp[0] == 3 || bp[0] == 4 ) {
      nextp = startp;
      if ( bp[0] != 4 ) {
        if ( !_mbs_GetBoundPointf ( bp, 0, p0, &normal, &startp, &b ) )
          goto out;
      }
    }
    else {
      if ( !_mbs_GetBoundPointf ( bp, 0, p0, &normal, &nextp, &b ) )
        goto out;
    }

    ni = maxinters-*ninters;
    a = _mbs_FindLinePolylineIntersf ( p0, t0, p1, t1, &normal, 3, 1,
             (float*)&lastp, a, &inters[*ninters], &ni );
    if ( ni < 0 )
      goto out;
    *ninters += ni;
                         /* cleanup the closed outline */
                         /* it is necessary to review the signs */
    if ( bufp[0] == 1 || bufp[0] == 3 || bp[0] == 4 ) {
      ni = *ninters - lni;
      _mbs_CleanLoopf ( &inters[lni], &ni );
      *ninters = (lni += ni);
      a = b;
    }
                         /* and move to the next item */
    bufp = bp;
  }
                         /* finally, sort the intersection points */
  if ( *ninters > 1 ) {
    if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(mbs_signpoint1f),
                        0, *ninters, inters ) != SORT_OK ) {
out:
      *ninters = -1;
    }
  }
#undef startp
#undef lastp
#undef nextp
} /*mbs_FindBoundLineIntersectionsf*/

