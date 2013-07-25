
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
boolean FindNearestCPoint ( int id, short x, short y )
{
  int     ncp, i, j;
  short   dist;
  boolean result;

  id &= 0x03;
  dist = xge_MINDIST;
  swind.current_tab = 0;
  switch ( win1_contents ) {
case WIN1_BLENDING:
    result = false;
    if ( display_control_net ) {
      ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
      result = FindNearestPoint ( ncp, rpoints[id], x, y, &dist,
                                  &swind.current_point );
    }
    if ( display_constr_poly ) {
      ncp = (lastknot_v-3*degree_v);
      for ( i = 1; i <= n_blending_constraints; i++ ) {
        j = (i-1)*(lastknot_v-degree_v) + degree_v;
        if ( FindNearestPoint ( ncp, &blending_constr_rp[id][j], x, y, &dist,
                                &swind.current_point ) ) {
          swind.current_tab = i;
          result = true;
        }
      }
    }
    return result;
default:
    if ( !display_control_net )
      return false;
    ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
    return FindNearestPoint ( ncp, rpoints[id], x, y, &dist,
                              &swind.current_point );
  }
} /*FindNearestCPoint*/

static void AuxCopyCPoint ( int i, int j )
{
  int id;

  cpoints[j] = cpoints[i];
  for ( id = 0; id < 4; id++ )
    rpoints[id][j] = rpoints[id][i];
} /*AuxCopyCPoint*/

void SetCPoint ( int id, short x, short y )
{
  point3d q, r;
  int     pitch, i, j, k, l;

  id &= 0x03;
  if ( swind.current_tab == 0 ) {  /* a surface patch control point */
    Point4to3d ( &cpoints[swind.current_point], &q );
    CameraProjectPoint3d ( &swind.CPos[id], &q, &r );
    r.x = (double)x;  r.y = (double)y;
    CameraUnProjectPoint3d ( &swind.CPos[id], &r, &q );
    Point3to4d ( &q, cpoints[swind.current_point].w, &cpoints[swind.current_point] );
    for ( id = 0; id < 4; id++ )
      clpoints[id][swind.current_point] =
        RzutujPunkt ( id, &cpoints[swind.current_point],
                      &rpoints[id][swind.current_point] );

    if ( kwind.closed_u || kwind.closed_v ) {
      pitch = lastknot_v-degree_v;
      i = swind.current_point / pitch;
      j = swind.current_point % pitch;
      k = l = swind.current_point;
      if ( kwind.closed_u ) {
        if ( i < degree_u ) k = swind.current_point+kwind.clcKu*pitch;
        else if ( i >= kwind.clcKu ) k = swind.current_point-kwind.clcKu*pitch;
        if ( k != swind.current_point )
          AuxCopyCPoint ( swind.current_point, k );
      }
      if ( kwind.closed_v ) {
        if ( j < degree_v ) l = swind.current_point+kwind.clcKv;
        else if ( j >= kwind.clcKv ) l = swind.current_point-kwind.clcKv;
        if ( l != swind.current_point ) {
          AuxCopyCPoint ( swind.current_point, l );
          if ( k != swind.current_point )
            AuxCopyCPoint ( swind.current_point, k+l-swind.current_point );
        }
      }
    }
  }
  else {  /* a blending surface constraint curve control point */
    j = (swind.current_tab-1)*(lastknot_v-degree_v) + degree_v + swind.current_point;
    Point4to3d ( &blending_constr_cp[j], &q );
    CameraProjectPoint3d ( &swind.CPos[id], &q, &r );
    r.x = (double)x;  r.y = (double)y;
    CameraUnProjectPoint3d ( &swind.CPos[id], &r, &q );
    Point3to4d ( &q, blending_constr_cp[j].w, &blending_constr_cp[j] );
    for ( id = 0; id < 4; id++ )
      clblending_constr[id][j] = RzutujPunkt ( id, &blending_constr_cp[j],
                                               &blending_constr_rp[id][j] );
  }
  if ( sw_bind_blending ) {
    if ( sw_blending_g1 )
      sw_bind_blending = ConstructG1BlendingSurface ();
    else
      sw_bind_blending = ConstructG2BlendingSurface ();
  }
} /*SetCPoint*/

static void MkPoint1 ( int i, int j, boolean mark )
{
  int pitch;

  pitch = lastknot_v-degree_v;
  if ( mark )
    mkpoints[i*pitch+j] |= 0x01;
  else
    mkpoints[i*pitch+j] &= 0xFE;
} /*MkPoint1*/

void MkPoint ( int i, int j, boolean mark )
{
  int k, l;

  MkPoint1 ( i, j, mark );
  if ( kwind.closed_u ) {
    if ( i < degree_u )    k = i+kwind.clcKu;
    else if ( i >= kwind.clcKu ) k = i-kwind.clcKu;
    else k = i;
    if ( k != i ) MkPoint1 ( k, j, mark );
    if ( kwind.closed_v ) {
      if ( j < degree_v )    l = j+kwind.clcKv;
      else if ( j >= kwind.clcKv ) l = j-kwind.clcKv;
      else l = j;
      if ( l != j ) {
        MkPoint1 ( i, l, mark );
        if ( k != i ) MkPoint1 ( k, l, mark );
      }
    }
  }
  else if ( kwind.closed_v ) {
    if ( j < degree_v )    l = j+kwind.clcKv;
    else if ( j >= kwind.clcKv ) l = j-kwind.clcKv;
    else l = j;
    if ( l != j ) MkPoint1 ( i, l, mark );
  }
} /*MkPoint*/

void MarkBlConstrPoint ( int i, boolean mark )
{
  if ( mark )
    mkblending_cp[i] |= 0x01;
  else
    mkblending_cp[i] &= 0xFE;
} /*MarkBlConstrPoint*/

void SelectPoints ( int id, const Box2s *sbox, boolean mk )
{
  int i, j, k, np, pitch;

  id &= 0x03;
  pitch = lastknot_v-degree_v;
  if ( sbox->x1-sbox->x0 < 2 && sbox->y1-sbox->y0 < 2 ) {
    if ( FindNearestCPoint ( id, sbox->x0, sbox->y0 ) ) {
      if ( swind.current_tab == 0 ) {
        i = swind.current_point / pitch;
        j = swind.current_point % pitch;
        MkPoint ( i, j, mk );
      }
      else
        MarkBlConstrPoint ( swind.current_point, mk );
    }
  }
  else {
    if ( display_control_net ) {
      np = (lastknot_u-degree_u)*pitch;
      for ( k = 0; k < np; k++ )
        if ( clpoints[id][k] &&
             rpoints[id][k].x >= sbox->x0 && rpoints[id][k].x <= sbox->x1 &&
             rpoints[id][k].y >= sbox->y0 && rpoints[id][k].y <= sbox->y1 ) {
          i = k / pitch;
          j = k % pitch;
          MkPoint ( i, j, mk );
        }
    }
    if ( display_constr_poly ) {
      for ( i = 1; i <= n_blending_constraints; i++ )
        if ( blending_constr_poly_valid[i] ) {
          k = (i-1)*pitch;
          for ( j = 3; j < pitch-3; j++ ) {
            if ( clblending_constr[id][k+j] &&
                 blending_constr_rp[id][k+j].x >= sbox->x0 &&
                 blending_constr_rp[id][k+j].x <= sbox->x1 &&
                 blending_constr_rp[id][k+j].y >= sbox->y0 &&
                 blending_constr_rp[id][k+j].y <= sbox->y1 )
            MarkBlConstrPoint ( k+j, mk );
          }
        }
    }
  }
} /*SelectPoints*/

