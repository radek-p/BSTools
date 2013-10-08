
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean DegreeElevationU ( void )
{
  void  *sp;
  int   ku, na, Na, r;
  double *ua, *cpa;
  int   id;

  sp = pkv_GetScratchMemTop ();
  sw_bind_blending = blending_mat_valid = false;
  if ( degree_u >= MAX_DEGREE )
    goto failure;

  r = mbs_KnotMultiplicityd ( min(2*degree_u+1,lastknot_u-1), &knots_u[1],
                              knots_u[degree_u] );
  ku = mbs_NumKnotIntervalsd ( degree_u, lastknot_u, knots_u );
  ua = pkv_GetScratchMemd ( lastknot_u+2+ku+r );
  if ( bind_spr ) {
    cpa = pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
    if ( !ua || !cpa )
      goto failure;
    if ( kwind.closed_u ) {
      if ( !mbs_BSDegElevClosedC3d ( degree_u, lastknot_u, knots_u,
                                     equator_cpoints, 1, &na, &Na, ua, cpa ) )
        goto failure;
    }
    else {
      if ( !mbs_BSDegElevC3d ( degree_u, lastknot_u, knots_u,
                               equator_cpoints, 1, &na, &Na, ua, cpa, true ) )
        goto failure;
    }
    memcpy ( equator_knots, ua, (Na+1)*sizeof(double) );
    memcpy ( equator_cpoints, cpa, (Na-na)*sizeof(point3d) );
    memset ( equator_mkpoints, 0, (Na-na)*sizeof(boolean) );
    eq_ckwind.degree = na;
    eq_ckwind.lastknot = Na;
    if ( kwind.closed_u ) {
      eq_ckwind.clcK = Na-2*na;
      eq_ckwind.clcT = equator_knots[na+eq_ckwind.clcK]-equator_knots[na];
    }
    BindSphericalProduct ();
  }
  else {
    cpa = pkv_GetScratchMem (
            (lastknot_u-degree_u+ku+r)*(lastknot_v-degree_v)*sizeof(point4d) );
    if ( !ua || !cpa )
      goto failure;

    if ( kwind.closed_u ) {
      if ( !mbs_multiBSDegElevClosedd ( 1, (lastknot_v-degree_v)*4, degree_u,
               lastknot_u, knots_u, 0, (double*)cpoints, 1, &na, &Na, ua,
               0, cpa ) )
        goto failure;
      kwind.clcKu = Na-2*na;
      kwind.clcTu = ua[na+kwind.clcKu]-ua[na];
    }
    else {
      if ( !mbs_multiBSDegElevd ( 1, (lastknot_v-degree_v)*4, degree_u,
               lastknot_u, knots_u, 0, (double*)cpoints, 1, &na, &Na, ua,
               0, cpa, false ) )
        goto failure;
    }
    degree_u = na;
    lastknot_u = Na;
    memcpy ( knots_u, ua, (lastknot_u+1)*sizeof(double) );
    memcpy ( cpoints, cpa,
             (lastknot_u-degree_u)*(lastknot_v-degree_v)*sizeof(point4d) );
    for ( id = 0; id < 4; id++ )
      ProjectSurface ( id );
  }
  pkv_SetScratchMemTop ( sp );

  ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v),
                      mkpoints );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeElevationU*/

boolean DegreeElevationV ( void )
{
  void  *sp;
  int   kv, ma, Ma, pitch1, pitch2, r, s, t;
  double *va, *cpa;
  int   id;

  sp = pkv_GetScratchMemTop ();
  sw_bind_blending = blending_mat_valid = false;
  if ( degree_v >= MAX_DEGREE )
    goto failure;
  r = mbs_KnotMultiplicityd ( min(2*degree_v+1,lastknot_v-1), &knots_v[1],
                              knots_v[degree_v] );
  for ( s = 0; knots_v[lastknot_v-degree_v-s-1] == knots_v[lastknot_v-degree_v]; s++ )
    ;
  for ( t = 0; knots_v[degree_v+t+1] == knots_v[degree_v]; t++ )
    ;
  kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  va = pkv_GetScratchMemd ( lastknot_v+2+kv+r );
  if ( bind_spr ) {
    cpa = pkv_GetScratchMem ( (lastknot_v-degree_v+kv+r)*sizeof(point3d) );
    if ( !va || !cpa )
      goto failure;
    if ( kwind.closed_v ) {
      if ( !mbs_BSDegElevClosedC3d ( degree_v, lastknot_v, knots_v,
                                     meridian_cpoints, 1, &ma, &Ma, va, cpa ) )
        goto failure;
    }
    else {
      if ( !mbs_BSDegElevC3d ( degree_v, lastknot_v, knots_v,
                               meridian_cpoints, 1, &ma, &Ma, va, cpa, true ) )
        goto failure;
    }
    memcpy ( meridian_knots, va, (Ma+1)*sizeof(double) );
    memcpy ( meridian_cpoints, cpa, (Ma-ma)*sizeof(point3d) );
    memset ( meridian_mkpoints, 0, (Ma-ma)*sizeof(boolean) );
    mer_ckwind.degree = ma;
    mer_ckwind.lastknot = Ma;
    if ( kwind.closed_v ) {
      mer_ckwind.clcK = Ma-2*ma;
      mer_ckwind.lastknot = meridian_knots[Ma-ma]-meridian_knots[ma];
    }
    BindSphericalProduct ();
  }
  else {
    cpa = pkv_GetScratchMem (
            (lastknot_u-degree_u)*(lastknot_v-degree_v+kv+r)*sizeof(point4d) );
    if ( !va || !cpa )
      goto failure;
    pitch1 = (lastknot_v-degree_v)*4;
    if ( kwind.closed_v ) {
      pitch2 = (lastknot_v-degree_v+kv+r-s-t)*4;
      if ( !mbs_multiBSDegElevClosedd ( lastknot_u-degree_u, 4, degree_v,
               lastknot_v, knots_v, pitch1,
               (double*)cpoints, 1, &ma, &Ma, va, pitch2, cpa ) )
        goto failure;
      kwind.clcKv = Ma-2*ma;
      kwind.clcTv = va[ma+kwind.clcKv]-va[ma];
    }
    else {
      pitch2 = (lastknot_v-degree_v+kv-s-t)*4;
      if ( !mbs_multiBSDegElevd ( lastknot_u-degree_u, 4, degree_v,
               lastknot_v, knots_v, pitch1, (double*)cpoints, 1,
               &ma, &Ma, va, pitch2, cpa, false ) )
        goto failure;
    }
    degree_v = ma;
    lastknot_v = Ma;
    memcpy ( knots_v, va, (lastknot_v+1)*sizeof(double) );
    memcpy ( cpoints, cpa,
             (lastknot_u-degree_u)*(lastknot_v-degree_v)*sizeof(point4d) );
    for ( id = 0; id < 4; id++ )
      ProjectSurface ( id );
  }
  pkv_SetScratchMemTop ( sp );

  ClearPointMarking (
      (lastknot_u-degree_u)*(lastknot_v-degree_v),
      mkpoints );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeElevationV*/

boolean DegreeReductionU ( void )
{
  void   *sp;
  double *nkn;
  double *ncp;
  int    deg, lkn, id;

  sp = pkv_GetScratchMemTop ();
  sw_bind_blending = blending_mat_valid = false;
  if ( degree_u > 1 ) {
    nkn = (double*)pkv_GetScratchMemd ( lastknot_u+1 );
    if ( bind_spr ) {
      ncp = pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
      if ( !nkn || !ncp )
        goto failure;
      if ( kwind.closed_u ) {
        if ( !mbs_BSDegRedClosedC3d ( eq_ckwind.degree, eq_ckwind.lastknot,
                                      equator_knots, equator_cpoints, 1,
                                      &deg, &lkn, nkn, ncp ) )
          goto failure;
      }
      else {
        if ( !mbs_BSDegRedC3d ( eq_ckwind.degree, eq_ckwind.lastknot,
                                equator_knots, equator_cpoints, 1,
                                &deg, &lkn, nkn, ncp ) )
          goto failure;
      }
      eq_ckwind.degree = deg;
      eq_ckwind.lastknot = lkn;
      memcpy ( equator_knots, nkn, (lkn+1)*sizeof(double) );
      memcpy ( equator_cpoints, ncp, (lkn-deg)*sizeof(point3d) );
      memset ( equator_mkpoints, 0, (lkn-deg)*sizeof(boolean) );
      if ( kwind.closed_u ) {
        eq_ckwind.clcK = lkn - 2*deg;
        eq_ckwind.clcT = equator_knots[lkn-deg] - equator_knots[deg];
      }
      BindSphericalProduct ();
    }
    else {
      ncp = pkv_GetScratchMem (
              (lastknot_u-degree_u)*(lastknot_v-degree_v)*sizeof(point4d) );
      if ( !nkn || !ncp )
        goto failure;
      if ( kwind.closed_u ) {
        if ( !mbs_multiBSDegRedClosedd ( 1, 4*(lastknot_v-degree_v),
                 degree_u, lastknot_u, knots_u, 0, (double*)cpoints,
                 1, &deg, &lkn, nkn, 0, (double*)ncp ) )
          goto failure;
        kwind.clcKu = lkn-2*deg;
      }
      else {
        if ( !mbs_multiBSDegRedd ( 1, 4*(lastknot_v-degree_v),
                 degree_u, lastknot_u, knots_u, 0, (double*)cpoints,
                 1, &deg, &lkn, nkn, 0, (double*)ncp ) )
          goto failure;
      }
      lastknot_u = lkn;
      degree_u = deg;
      memcpy ( knots_u, nkn, (lkn+1)*sizeof(double) );
      memcpy ( cpoints, ncp, (lkn+1)*(lastknot_v-degree_v)*sizeof(point4d) );
      for ( id = 0; id < 4; id++ )
        ProjectSurface ( id );
    }
    ClearPointMarking (
        (lastknot_u-degree_u)*(lastknot_v-degree_v),
        mkpoints );
    pkv_SetScratchMemTop ( sp );
    return true;
  }

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeReductionU*/

boolean DegreeReductionV ( void )
{
  void    *sp;
  double   *nkn;
  point4d *ncp;
  int     deg, lkn, id, ncurv, pitch;

  sp = pkv_GetScratchMemTop ();
  sw_bind_blending = blending_mat_valid = false;
  if ( degree_v > 1 ) {
    nkn = (double*)pkv_GetScratchMemd ( lastknot_v+1 );
    if ( bind_spr ) {
      ncp = pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
      if ( !nkn || !ncp )
        goto failure;
      if ( kwind.closed_v ) {
        if ( !mbs_BSDegRedClosedC3d ( mer_ckwind.degree, mer_ckwind.lastknot,
                                      meridian_knots, meridian_cpoints, 1,
                                      &deg, &lkn, nkn, ncp ) )
          goto failure;
      }
      else {
        if ( !mbs_BSDegRedC3d ( mer_ckwind.degree, mer_ckwind.lastknot,
                                meridian_knots, meridian_cpoints, 1,
                                &deg, &lkn, nkn, ncp ) )
          goto failure;
      }
      mer_ckwind.degree= deg;
      mer_ckwind.lastknot = lkn;
      memcpy ( meridian_knots, nkn, (lkn+1)*sizeof(double) );
      memcpy ( meridian_cpoints, ncp, (lkn-deg)*sizeof(point3d) );
      memset ( meridian_mkpoints, 0, (lkn-deg)*sizeof(boolean) );
      if ( kwind.closed_u ) {
        mer_ckwind.clcK = lkn-2*deg;
        mer_ckwind.clcT = meridian_knots[lkn-deg]-meridian_knots[deg];
      }
      BindSphericalProduct ();
    } 
    else {
      ncurv = lastknot_u-degree_u;
      pitch = 4*(lastknot_v-degree_v);
      ncp = (point4d*)pkv_GetScratchMemd ( ncurv*pitch );
      if ( !nkn || !ncp )
        goto failure;
      if ( kwind.closed_v ) {
        if ( !mbs_multiBSDegRedClosedd ( ncurv, 4, degree_v, lastknot_v, knots_v,
                pitch, (double*)cpoints, 1, &deg, &lkn, nkn, pitch,
                (double*)ncp ) )
          goto failure;
        kwind.clcKv = lkn-2*deg;
      }
      else {
        if ( !mbs_multiBSDegRedd ( ncurv, 4, degree_v, lastknot_v, knots_v,
                pitch, (double*)cpoints, 1, &deg, &lkn, nkn, pitch,
                (double*)ncp ) )
          goto failure;
      }
      lastknot_v = lkn;
      degree_v = deg;
      memcpy ( knots_v, nkn, (lkn+1)*sizeof(double) );
      pkv_Selectd ( ncurv, 4*(lkn-deg), pitch, 4*(lkn-deg), ncp, cpoints );
      ClearPointMarking (
          (lastknot_u-degree_u)*(lastknot_v-degree_v),
          mkpoints );
      for ( id = 0; id < 4; id++ )
        ProjectSurface ( id );
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeReductionV*/

