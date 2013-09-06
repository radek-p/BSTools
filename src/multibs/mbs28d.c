
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

/* ////////////////////////////////////////////////////////////// */
/* Computing intersections of lines with a trimmed patch boundary */

/* The boundary processed by the procedure mbs_FindBoundLineIntersectionsf */
/* must be described by a code produced by mbs_CompileTrimPatchBoundf.     */
 
static double _get_a3 ( const point2d *p0, const vector2d *normal,
                        const point3d *cp )
{
  vector2d v;
  double    w;

  w = cp->z;
  SetVector2d ( &v, cp->x - w*p0->x, cp->y - w*p0->y );
  return DotProduct2d ( normal, &v );
} /*_get_a3*/

static double _get_a2 ( const point2d *p0, const vector2d *normal,
                        const point2d *cp )
{
  point3d q;

  SetPoint3d ( &q, cp->x, cp->y, 1.0 );
  return _get_a3 ( p0, normal, &q );
} /*_get_a2*/

static double _get_lparamd ( const point2d *p0, double t0,
                             const point2d *p1, double t1,
                             const point2d *q )
{
  vector2d a, b;
  double    v;

  SubtractPoints2d ( p1, p0, &a );
  SubtractPoints2d (  q, p0, &b );
  v = DotProduct2d ( &a, &b )/DotProduct2d ( &a, &a );
  return t0 + v*(t1-t0);
} /*_get_lparamd*/

static void _mbs_AddLineLineIntersd ( const point2d *p0, double t0,
                                      const point2d *p1, double t1,
                                      signpoint1d *inters,
                                      char dim, const double *vert,
                                      double u, char sa )
{
  point2d p;
  point3d q;

  if ( dim == 2 )
    InterPoint2d ( (point2d*)vert, (point2d*)&vert[2], u, &p );
  else {
    InterPoint3d ( (point3d*)vert, (point3d*)&vert[3], u, &q );
    Point3to2d ( &q, &p );
  }
  inters->t = _get_lparamd ( p0, t0, p1, t1, &p );
  inters->sign1 = sa;
} /*_mbs_AddLineLineIntersd*/

static double _mbs_FindLinePolylineIntersd ( const point2d *p0, double t0,
                                             const point2d *p1, double t1,
                                             const vector2d *normal,
                                             char dim, short nlines,
                                             const double *vert,
                                             double a,
                                             signpoint1d *inters, int *ninters )
{
  int    i, ni, maxinters;
  double b;
  char   sa, sb;

  maxinters = *ninters;
  *ninters = ni = 0;
  sa = pkv_signd ( a );
  for ( i = 1; i <= nlines; i++ ) {
    if ( dim == 2 )
      b = _get_a2 ( p0, normal, (point2d*)(&vert[2*i]) );
    else
      b = _get_a3 ( p0, normal, (point3d*)(&vert[3*i]) );
    sb = pkv_signd ( b );
    if ( sb != sa ) {   /* there is perhaps an intersection */
      if ( ni >= maxinters )
        goto out;
      _mbs_AddLineLineIntersd ( p0, t0, p1, t1, &inters[ni],
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
} /*_mbs_FindLinePolylineIntersd*/

static int   _deg;   /* polynomial degree, for function _mbs_apoly et al. */
static double *_ac;  /* polynomial coefficients */
static int   _sp;    /* stack pointer, 0 if stack empty */
static double *_sd;  /* stack data pointer */

static double _mbs_apoly ( double t )
{
  double value;

  mbs_BCHornerC1d ( _deg, _ac, t, &value );
  return value;
} /*_mbs_apoly*/

static boolean _mbs_varsign ( double *ac )
{
  int  i;
  char s, z;

  s = pkv_signd ( ac[0] );
  for ( i = 1; i <= _deg; i++ ) {
    z = pkv_signd ( ac[i] );
    if ( z != s )
      return true;
  }
  return false;
} /*_mbs_varsign*/

static boolean _mbs_onesignch ( double *ac, boolean *nonzero )
{
  int  i, k;
  char sa, sb;

  k = 0;
  sa = pkv_signd ( ac[0] );
  *nonzero = (boolean)(sa != 0);
  for ( i = 1; i <= _deg; i++ ) {
    sb = pkv_signd ( ac[i] );
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

static void _mbs_pushp ( double t0, double t1, double *ac )
{
  _sd = pkv_GetScratchMemd ( _deg+3 );
  if ( !_sd ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
    exit ( 1 );
  }
  _sd[0] = t0;
  _sd[1] = t1;
  memcpy ( &_sd[2], ac, (_deg+1)*sizeof(double) );
  _sp++;
} /*_mbs_pushp*/

static void _mbs_popp ( double *t0, double *t1, double *ac )
{
  if ( !_sp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_6, ERRMSG_6 );
    exit ( 1 );
  }
  _sp--;
  *t0 = _sd[0];
  *t1 = _sd[1];
  memcpy ( ac, &_sd[2], (_deg+1)*sizeof(double) );
  _sd -= _deg+3;
  pkv_FreeScratchMemd ( _deg+3 );
} /*_mbs_popp*/

static void _mbs_AddLineBezcIntersd ( const point2d *p0, double t0,
                                      const point2d *p1, double t1,
                                      signpoint1d *inters,
                                      char dim, int deg, const double *cpts,
                                      double u, char as )
{
  point2d q;

  if ( dim == 2 )
    mbs_BCHornerC2d ( deg, (point2d*)cpts, u, &q );
  else
    mbs_BCHornerC2Rd ( deg, (point3d*)cpts, u, &q );
  inters->t = _get_lparamd ( p0, t0, p1, t1, &q );
  inters->sign1 = as;
} /*_mbs_AddLineBezcIntersd*/

static double _mbs_FindLineBezcIntersd ( const point2d *p0, double t0,
                                         const point2d *p1, double t1,
                                         const vector2d *normal,
                                         char dim, short deg,
                                         const double *cpts,
                                         double a,
                                         signpoint1d *inters, int *ninters )
{
#define eps 1.0e-6
  void    *ssp;
  int     i, ni, maxinters;
  double  *ac, *bc;
  double  u0, u1, u;
  boolean nonzero, error;

  ssp = pkv_GetScratchMemTop ();

                        /* compute the polynomial coefficients */
  ac = pkv_GetScratchMemd ( 3*(deg+1) );
  if ( !ac )
    goto out;
  _ac = &ac[deg+1];
  bc  = &_ac[deg+1];

  ac[0] = a;
  if ( dim == 2 ) {
    for ( i = 1; i <= deg; i++ )
      ac[i] = _get_a2 ( p0, normal, (point2d*)&cpts[2*i] );
  }
  else { /* dim == 3 */
    for ( i = 1; i <= deg; i++ )
      ac[i] = _get_a3 ( p0, normal, (point3d*)&cpts[3*i]);
  }
  a = ac[deg];

                        /* find the polynomial zeros */
  maxinters = *ninters;
  *ninters = ni = 0;
  _deg = deg;
  _sp = 0;
  if ( _mbs_varsign ( ac ) ) {
    _mbs_pushp ( 0.0, 1.0, ac );
    do {
      _mbs_popp ( &u0, &u1, _ac );
      if ( _mbs_onesignch ( _ac, &nonzero ) ) {
        if ( nonzero ) {
          u = pkn_Illinoisd ( _mbs_apoly, 0.0, 1.0, eps, &error );
          u = u0 + u*(u1-u0);

          if ( ni >= maxinters )
            goto out;
          _mbs_AddLineBezcIntersd ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, u, pkv_signd(_ac[0]) );
          ni++;
        }
        else {
          if ( ni >= maxinters-1 )
            goto out;
          _mbs_AddLineBezcIntersd ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, u0, 0 );
          _mbs_AddLineBezcIntersd ( p0, t0, p1, t1, &inters[ni+1],
                                    dim, deg, cpts, u1, 0 );
          ni += 2;
        }
      }
      else {
        if ( u1-u0 > eps ) {
          mbs_BisectBC1d ( deg, _ac, bc );
          u = 0.5*(u0+u1);
          if ( _mbs_varsign ( _ac ) )
            _mbs_pushp ( u, u1, _ac );
          if ( _mbs_varsign ( bc ) )
            _mbs_pushp ( u0, u, bc );
        }
        else {
          if ( ni >= maxinters )
            goto out;
          _mbs_AddLineBezcIntersd ( p0, t0, p1, t1, &inters[ni],
                                    dim, deg, cpts, 0.5*(u0+u1),
                                    pkv_signd(_ac[0]) );
          ni++;
        }
      }
    } while ( _sp );
  }

  pkv_SetScratchMemTop ( ssp );
  *ninters = ni;
  return a;

out:
  pkv_SetScratchMemTop ( ssp );
  *ninters = -1;
  return 0.0;
#undef eps
} /*_mbs_FindLineBezcIntersd*/

static boolean _mbs_GetBoundPointd ( const char *bufp, short k,
                                     const point2d *p0, const vector2d *normal,
                                     point3d *point, double *a )
{
  point2d *q;

/* it computes the number a, in order to avoid */
/* the nasty consequences of rounding errors.  */

  if ( bufp[1] == 2 ) {
    q = ((point2d*)&bufp[2+sizeof(short)+2*k*sizeof(double)]);
    SetPoint3d ( point, q->x, q->y, 1.0 );
  }
  else if ( bufp[1] == 3 )
    *point = *(point3d*)&bufp[2+sizeof(short)+3*k*sizeof(double)];
  else
    return false;
  *a = _get_a3 ( p0, normal, point );
  return true;
} /*_mbs_GetBoundPointd*/

static void _mbs_CleanLoopd ( signpoint1d *inters, int *ninters )
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
} /*_mbs_CleanLoopd*/

void mbs_FindBoundLineIntersectionsd ( const void *bound,
                                       const point2d *p0, double t0,
                                       const point2d *p1, double t1,
                                       signpoint1d *inters, int *ninters )
{
  int      maxinters, ni, lni;
  short    N;
  vector2d normal;
  point3d  slp[3];
#define startp slp[0]
#define lastp  slp[1]
#define nextp  slp[2]
  char     *bufp, *bp, dim;
  double    *points, a, b = 0.0;

  SubtractPoints2d ( p1, p0, &normal );
  SetVector2d ( &normal, -normal.y, normal.x );

  maxinters = *ninters;
  *ninters = lni = 0;

  bufp = (char*)bound;
  if ( !_mbs_GetBoundPointd ( bufp, 0, p0, &normal, &startp, &a ) )
    goto out;

  while ( bufp[0] != 4 ) {
    dim    = bufp[1];
    N      = *((short*)(&bufp[2]));
    if ( dim < 2 || dim > 3 || N < 1 )
      goto out;
    points = ((double*)(&bufp[2+sizeof(short)]));

    ni = maxinters-*ninters;
    switch ( bufp[0] ) {
 case 0:
 case 1:
      a = _mbs_FindLinePolylineIntersd ( p0, t0, p1, t1, &normal, dim, N,
                                         points, a, &inters[*ninters], &ni );
      break;

 case 2:
 case 3:
      a = _mbs_FindLineBezcIntersd ( p0, t0, p1, t1, &normal, dim, N,
                                     points, a, &inters[*ninters], &ni );
      break;

 default:
      goto out;          /* signal invalid data */
    }
    if ( ni < 0 )
      goto out;
    *ninters += ni;

                         /* now deal with the joining line */
    if ( !_mbs_GetBoundPointd ( bufp, N, p0, &normal, &lastp, &a ) )
      goto out;
    bp = &bufp[2+sizeof(short)+(N+1)*dim*sizeof(double)];
    if ( bufp[0] == 1 || bufp[0] == 3 || bp[0] == 4 ) {
      nextp = startp;
      if ( bp[0] != 4 ) {
        if ( !_mbs_GetBoundPointd ( bp, 0, p0, &normal, &startp, &b ) )
          goto out;
      }
    }
    else {
      if ( !_mbs_GetBoundPointd ( bp, 0, p0, &normal, &nextp, &b ) )
        goto out;
    }

    ni = maxinters-*ninters;
    a = _mbs_FindLinePolylineIntersd ( p0, t0, p1, t1, &normal, 3, 1,
             (double*)&lastp, a, &inters[*ninters], &ni );
    if ( ni < 0 )
      goto out;
    *ninters += ni;
                         /* cleanup the closed outline */
                         /* it is necessary to review the signs */
    if ( bufp[0] == 1 || bufp[0] == 3 || bp[0] == 4 ) {
      ni = *ninters - lni;
      _mbs_CleanLoopd ( &inters[lni], &ni );
      *ninters = (lni += ni);
      a = b;
    }
                         /* and move to the next item */
    bufp = bp;
  }
                         /* finally, sort the intersection points */

  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(signpoint1d),
                      0, *ninters, inters ) != SORT_OK ) {
out:
    *ninters = -1;
  }

#undef startp
#undef lastp
#undef nextp
} /*mbs_FindBoundLineIntersectionsd*/

