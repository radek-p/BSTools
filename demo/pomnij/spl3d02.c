
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
boolean InsertKnotU ( double newknot )
{
  if ( bind_spr ) {
    if ( kwind.closed_u )
      mbs_KnotInsClosedC3d ( degree_u, &lastknot_u, equator_knots,
                             equator_cpoints, newknot );
    else
      mbs_KnotInsC3d ( degree_u, &lastknot_u, equator_knots,
                       equator_cpoints, newknot );
    if ( kwind.closed_u ) {
      eq_ckwind.clcK = kwind.clcKu = lastknot_u - 2*degree_u;
      eq_ckwind.clcT = kwind.clcTu =
          equator_knots[lastknot_u-degree_u]-equator_knots[degree_u];
    }
    BindSphericalProduct ();
  }
  else {
    if ( kwind.closed_u ) {
      mbs_multiKnotInsClosedd ( degree_u, &lastknot_u, knots_u, 1,
                                (lastknot_v-degree_v)*4, 0, 0,
                                (double*)cpoints, newknot );
      kwind.clcKu = lastknot_u-2*degree_u;
    }
    else {
      mbs_multiKnotInsd ( degree_u, &lastknot_u, knots_u, 1,
                          (lastknot_v-degree_v)*4, 0, 0,
                          (double*)cpoints, newknot );
    }
    kwind.lastknot_u = lastknot_u;
    if ( sw_bind_blending ) {
      SetEquidistantU ();
      blending_mat_valid = false;
      if ( sw_blending_g1 )
        ConstructG1BlendingSurface ();
      else
        ConstructG2BlendingSurface ();
    }
    ResizeObject ();
  }
  ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v), mkpoints );
  return true;
} /*InsertKnotU*/

boolean InsertKnotV ( double newknot )
{
  if ( bind_spr ) {
    if ( kwind.closed_v )
      mbs_KnotInsClosedC3d ( degree_v, &lastknot_v, meridian_knots,
                             meridian_cpoints, newknot );
    else
      mbs_KnotInsC3d ( degree_v, &lastknot_v, meridian_knots,
                       meridian_cpoints, newknot );
    if ( kwind.closed_v ) {
      mer_ckwind.clcK = kwind.clcKv = lastknot_v - 2*degree_v;
      mer_ckwind.clcT = kwind.clcTv =
          meridian_knots[lastknot_v-degree_v]-meridian_knots[degree_v];
    }
    BindSphericalProduct ();
  }
  else {
    if ( kwind.closed_v ) {
      mbs_multiKnotInsClosedd ( degree_v, &lastknot_v, knots_v,
               lastknot_u-degree_u, 4,
               (lastknot_v-degree_v)*4,
               (lastknot_v-degree_v+1)*4,
               (double*)cpoints, newknot );
      kwind.clcKv = lastknot_v-2*degree_v;
    }
    else {
      mbs_multiKnotInsd ( degree_v, &lastknot_v, knots_v,
                          lastknot_u-degree_u, 4,
                          (lastknot_v-degree_v)*4,
                          (lastknot_v-degree_v+1)*4,
                          (double*)cpoints, newknot );
    }
    kwind.lastknot_v = lastknot_v;
    if ( sw_bind_blending ) {
      SetEquidistantV ();
      blending_mat_valid = false;
      if ( sw_blending_g1 )
        ConstructG1BlendingSurface ();
      else
        ConstructG2BlendingSurface ();
    }
    ResizeObject ();
  }
  ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v), mkpoints );
  return true;
} /*InsertKnotV*/

boolean RemoveKnotU ( int knotnum )
{
  if ( bind_spr ) {
    if ( kwind.closed_u )
      mbs_KnotRemoveClosedC3d ( degree_u, &lastknot_u, equator_knots,
                                equator_cpoints, knotnum );
    else
      mbs_KnotRemoveC3d ( degree_u, &lastknot_u, equator_knots,
                          equator_cpoints, knotnum );
    if ( kwind.closed_u ) {
      eq_ckwind.clcK = kwind.clcKu = lastknot_u - 2*degree_u;
      eq_ckwind.clcT = kwind.clcTu =
          equator_knots[lastknot_u-degree_u]-equator_knots[degree_u];
    }
    BindSphericalProduct ();
  }
  else {
    if ( sw_bind_blending ) {
      if ( (sw_blending_g1 && lastknot_u <= 7) || 
           (sw_blending_g2 && lastknot_u <= 10) )
      return false;
    }
    if ( kwind.closed_u ) {
      mbs_multiKnotRemoveClosedd ( degree_u, &lastknot_u, knots_u,
                                   1, (lastknot_v-degree_v)*4, 0, 0,
                                   (double*)cpoints, knotnum );
      kwind.clcKu = lastknot_u-2*degree_u;
    }
    else {
      mbs_multiKnotRemoved ( degree_u, &lastknot_u, knots_u,  
                             1, (lastknot_v-degree_v)*4, 0, 0,  
                             (double*)cpoints, knotnum );
    }
    kwind.lastknot_u = lastknot_u;
    if ( sw_bind_blending ) {
      SetEquidistantU ();
      blending_mat_valid = false;
      if ( sw_blending_g1 )
        ConstructG1BlendingSurface ();
      else
        ConstructG2BlendingSurface ();
    }
    ResizeObject ();
  }
  ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v), mkpoints );
  return true;
} /*RemoveKnotU*/

boolean RemoveKnotV ( int knotnum )
{
  if ( bind_spr ) {
    if ( kwind.closed_v )
      mbs_KnotRemoveClosedC3d ( degree_v, &lastknot_v, meridian_knots,
                                meridian_cpoints, knotnum );
    else
      mbs_KnotRemoveC3d ( degree_v, &lastknot_v, meridian_knots,
                          meridian_cpoints, knotnum );
    if ( kwind.closed_v ) {
      mer_ckwind.clcK = kwind.clcKv = lastknot_v - 2*degree_v;
      mer_ckwind.clcT = kwind.clcTv =
          equator_knots[lastknot_v-degree_v]-equator_knots[degree_v];
    }
    BindSphericalProduct ();
  }
  else {
    if ( sw_bind_blending ) {
      if ( (sw_blending_g1 && lastknot_v <= 7) || 
           (sw_blending_g2 && lastknot_v <= 10) )
      return false;
    }
    if ( kwind.closed_v ) {
      mbs_multiKnotRemoveClosedd ( degree_v, &lastknot_v, knots_v,
                                   lastknot_u-degree_u, 4,
                                   (lastknot_v-degree_v)*4,
                                   (lastknot_v-degree_v-1)*4,
                                   (double*)cpoints, knotnum );
      kwind.clcKv = lastknot_u-2*degree_u;
    }
    else {
      mbs_multiKnotRemoved ( degree_v, &lastknot_v, knots_v,
                             lastknot_u-degree_u, 4,
                             (lastknot_v-degree_v)*4,
                             (lastknot_v-degree_v-1)*4,
                             (double*)cpoints, knotnum );
    }
    kwind.lastknot_v = lastknot_v;
    if ( sw_bind_blending ) {
      blending_mat_valid = false;
      SetEquidistantV ();
      if ( sw_blending_g1 )
        ConstructG1BlendingSurface ();
      else
        ConstructG2BlendingSurface ();
    }
    ResizeObject ();
  }
  ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v), mkpoints );
  return true;
} /*RemoveKnotV*/

void SetEquidistantU ( void )
{
  int i;

  for ( i = 0; i <= lastknot_u; i++ )
    knots_u[i] = (double)i;
  xge_T2KnotWindFindMapping ( &kwind );
  if ( bind_spr )
    memcpy ( equator_knots, knots_u, (lastknot_u+1)*sizeof(double) );
  if ( kwind.closed_u ) {
    kwind.clcTu = lastknot_u-2*degree_u;
    if ( bind_spr )
      eq_ckwind.clcT = kwind.clcTu;
  }
  ResizeObject ();
} /*SetEquidistantU*/

void SetEquidistantV ( void )
{
  int i;

  for ( i = 0; i <= lastknot_v; i++ )
    knots_v[i] = (double)i;
  xge_T2KnotWindFindMapping ( &kwind );
  if ( bind_spr )
    memcpy ( meridian_knots, knots_v, (lastknot_v+1)*sizeof(double) );
  if ( kwind.closed_v ) {
    kwind.clcTv = lastknot_v-2*degree_v;
    if ( bind_spr )
      mer_ckwind.clcT = kwind.clcTv;
  }
  ResizeObject ();
} /*SetEquidistantV*/

