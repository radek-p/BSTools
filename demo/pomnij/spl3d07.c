
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
static void GetClipLines ( CameraRecd *CPos, vector3d cliplines[4] )
{
  SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-CPos->xmin) );
  SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-CPos->ymin) );
  SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(CPos->xmin+CPos->width) );
  SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(CPos->ymin+CPos->height) );
} /*GetClipLines*/

static void RzutujPK ( CameraRecd *CPos, point4d *p, point3d *q )
{
  point4d r;

  Trans3Point4d ( &CPos->CTr, p, &r );
  if ( CPos->parallel )
    SetPoint3d ( q, r.x, r.y, r.w );
  else
    SetPoint3d ( q, r.x+CPos->vd.persp.xi0*r.z,
                 r.y+CPos->vd.persp.eta0*r.z, r.z );
} /*RzutujPK*/

void DisplayBezCurve4 ( CameraRecd *CPos, int degree, point4d *cp )
{
  void     *sp;
  int      i;  
  point3d  *cq;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
  cq = (point3d*)pkv_GetScratchMem ( (degree+1)*sizeof(point3d) );
  if ( !cp )
    { printf ( "*" );  exit ( 1 ); }
  GetClipLines ( CPos, cliplines );   
  for ( i = 0; i <= degree; i++ )
    RzutujPK ( CPos, &cp[i], &cq[i] );
  mbs_ClipBC2Rd ( 4, cliplines, degree, cq, xge_DrawBC2Rd );
  pkv_SetScratchMemTop ( sp );
} /*DisplayBezCurve4*/

void DisplayBezPatch4 ( CameraRecd *CPos,
                        int degree_u, int degree_v, point4d *cp,
                        boolean u_edge, boolean v_edge )
{
  void    *sp;
  point4d *c;
  int     i, k, pitch;
  double   t;

  sp = pkv_GetScratchMemTop ();
  c = pkv_GetScratchMem ( (max(degree_u,degree_v)+1)*sizeof(point4d) );
  if ( c ) {
    pitch = (degree_v+1)*4;
    if ( u_edge )
      k = display_bez_dens_u;
    else
      k = display_bez_dens_u-1;
    for ( i = 0; i <= k; i++ ) {
      if ( i == 0 || i == display_bez_dens_u )
        xgeSetForeground ( xgec_White );
      else if ( i == 1 )
        xgeSetForeground ( xgec_Grey3 );
      t = (double)i/(double)display_bez_dens_u;
      mbs_multiBCHornerd ( degree_u, 1, pitch, 0, (double*)cp, t, (double*)c );
      DisplayBezCurve4 ( CPos, degree_v, c );
    }
    if ( v_edge )
      k = display_bez_dens_v;
    else
      k = display_bez_dens_v-1;
    for ( i = 0; i <= k; i++ ) {
      if ( i == 0 || i == display_bez_dens_v )
        xgeSetForeground ( xgec_White );
      else if ( i == 1 )
        xgeSetForeground ( xgec_Grey3 );
      t = (double)i/(double)display_bez_dens_v;
      mbs_multiBCHornerd ( degree_v, degree_u+1, 4, pitch, (double*)cp, t,
                           (double*)c );
      DisplayBezCurve4 ( CPos, degree_u, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayBezPatch4*/

void DisplaySurface ( int id )
{
  void   *sp;
  int    ku, kv, pitch1, pitch2, pitch3;
  double  *b, *c;
  int    i, j, start;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  ku = mbs_NumKnotIntervalsd ( degree_u, lastknot_u, knots_u );
  kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  pitch1 = (lastknot_v-degree_v)*4;
  pitch2 = (degree_v+1)*4*kv;
  pitch3 = (degree_v+1)*4;
  b = pkv_GetScratchMemd ( pitch2*ku*(degree_u+1)*4 );
  c = pkv_GetScratchMemd ( (degree_u+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezd ( 4, degree_u, lastknot_u, knots_u,
                        degree_v, lastknot_v, knots_v,
                        pitch1, (double*)cpoints,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(degree_u+1)*pitch2+pitch3*j;
        pkv_Selectd ( degree_u+1, pitch3, pitch2, pitch3, &b[start], c );
        DisplayBezPatch4 ( &swind.CPos[id], degree_u, degree_v, (point4d*)c,
                           (boolean)(i == ku-1), (boolean)(j == kv-1) );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplaySurface*/

void DisplayConstraintCurves ( int id )
{
  void   *sp;
  int    i, j, kv;
  double *ccp, *ccb, t;

  if ( n_blending_constraints <= 0 )
    return;
  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  ccp = pkv_GetScratchMemd ( (lastknot_v-degree_v)*4 );
  ccb = pkv_GetScratchMemd ( kv*(degree_v+1)*4 );
  if ( ccp && ccb ) {
    xgeSetForeground ( xgec_Orchid1 );
    for ( i = 1; i <= n_blending_constraints; i++ ) {
      t = blending_constr_knots[i];
      if ( t > knots_u[degree_u] && t < knots_u[lastknot_u-degree_u] ) {
        mbs_multideBoord ( degree_u, lastknot_u, knots_u, 1,
                           4*(lastknot_v-degree_v), 4*(lastknot_v-degree_v),
                           (double*)cpoints, t, ccp );
        mbs_BSToBezC4d ( degree_v, lastknot_v, knots_v, ccp, NULL, NULL, NULL, ccb );
        for ( j = 0; j < kv; j++ )
          DisplayBezCurve4 ( &swind.CPos[id], degree_v,
                             (point4d*)&ccb[4*j*(degree_v+1)] );
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayConstraintCurves*/

void DisplayLine3 ( int id, point3d *p0, point3d *p1 )
{
  point3d r0, r1;

  id &= 0x03;
  if ( CameraClipLine3d ( &swind.CPos[id], p0, 0.0, p1, 1.0, &r0, &r1 ) ) {
    xgeDrawLine ( (int)(r0.x+0.5), (int)(r0.y+0.5),  
                  (int)(r1.x+0.5), (int)(r1.y+0.5) );
  }
} /*DisplayLine4*/

void DisplayLine4 ( int id, point4d *p0, point4d *p1 )
{
  point3d q0, q1, r0, r1;

  id &= 0x03;
  Point4to3d ( p0, &q0 );
  Point4to3d ( p1, &q1 );
  if ( CameraClipLine3d ( &swind.CPos[id], &q0, 0.0, &q1, 1.0, &r0, &r1 ) ) {
    xgeDrawLine ( (int)(r0.x+0.5), (int)(r0.y+0.5),  
                  (int)(r1.x+0.5), (int)(r1.y+0.5) );
  }
} /*DisplayLine4*/

void DisplayControlNet ( int id )
{
  int i, j, k, pitch;

  id &= 0x03;
  xgeSetForeground ( xgec_Green );
  pitch = lastknot_v-degree_v;
  for ( i = k = 0; i < lastknot_u-degree_u; i++ ) {
    for ( j = 0; j < lastknot_v-degree_v-1; j++, k++ ) {
      if ( clpoints[id][k] && clpoints[id][k+1] )
        xgeDrawLine ( (int)(rpoints[id][k].x+0.5), (int)(rpoints[id][k].y+0.5),
              (int)(rpoints[id][k+1].x+0.5), (int)(rpoints[id][k+1].y+0.5) );
      else
        DisplayLine4 ( id, &cpoints[k], &cpoints[k+1] );
    }
    k++;
  }
  for ( i = k = 0; i < lastknot_u-degree_u-1; i++ ) {
    for ( j = 0; j < lastknot_v-degree_v; j++, k++ ) {
      if ( clpoints[id][k] && clpoints[id][k+pitch] )
        xgeDrawLine ( (int)(rpoints[id][k].x+0.5), (int)(rpoints[id][k].y+0.5),
            (int)(rpoints[id][k+pitch].x+0.5),
            (int)(rpoints[id][k+pitch].y+0.5) );
      else
        DisplayLine4 ( id, &cpoints[k], &cpoints[k+pitch] );
    }
  }
} /*DisplayControlNet*/

static void DisplayBezNet4 ( int id, int n, int m, point4d *cp )
{
  int i, j;

  id &= 0x03;
  xgeSetForeground ( xgec_SeaGreen3 );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j < m; j++ )
      DisplayLine4 ( id, &cp[(m+1)*i+j], &cp[(m+1)*i+j+1] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      DisplayLine4 ( id, &cp[(m+1)*i+j], &cp[(m+1)*(i+1)+j] );
} /*DisplayBezNet4*/

void DisplayBezierNets ( int id )
{
  void   *sp;
  int    ku, kv, pitch1, pitch2, pitch3;
  double  *b, *c;
  int    i, j, start;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  ku = mbs_NumKnotIntervalsd ( degree_u, lastknot_u, knots_u );
  kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  pitch1 = (lastknot_v-degree_v)*4;
  pitch2 = (degree_v+1)*4*kv;
  pitch3 = (degree_v+1)*4;
  b = pkv_GetScratchMemd ( pitch2*ku*(degree_u+1)*4 );
  c = pkv_GetScratchMemd ( (degree_u+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezd ( 4, degree_u, lastknot_u, knots_u,
                        degree_v, lastknot_v, knots_v,
                        pitch1, (double*)cpoints,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(degree_u+1)*pitch2+pitch3*j;
        pkv_Selectd ( degree_u+1, pitch3, pitch2, pitch3, &b[start], c );
        DisplayBezNet4 ( id, degree_u, degree_v, (point4d*)c );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayBezierNets*/

void DisplayControlPoints ( int id )
{
  int i, j, k;

  id &= 0x03;
  xgeSetForeground ( xgec_Yellow );
  if ( win1_contents == WIN1_BLENDING ) {
    for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
      for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
        if ( clpoints[id][k] && !(mkpoints[k] & 0x01) )
          xgeFillRectangle ( 3, 3,
              (int)(rpoints[id][k].x-0.5), (int)(rpoints[id][k].y-0.5) );
      }
    xgeSetForeground ( xgec_OrangeRed );
    for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
      for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
        if ( clpoints[id][k] && (mkpoints[k] & 0x01) )
          xgeFillRectangle ( 3, 3,
              (int)(rpoints[id][k].x-0.5), (int)(rpoints[id][k].y-0.5) );
      }
    xgeSetForeground ( xgec_Green );
    if ( kwind.closed_u ) {
      for ( i = k = 0; i < lastknot_u-2*degree_u;  i++ )
        for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
          if ( j >= blending_opt_part[2] && j <= blending_opt_part[3] &&
               clpoints[id][k] ) {
            if ( blending_opt_part[0] <= blending_opt_part[1] ) {
              if ( i >= blending_opt_part[0] && i <= blending_opt_part[1] )
                xgeFillRectangle ( 3, 3,
                    (int)(rpoints[id][k].x-0.5), (int)(rpoints[id][k].y-0.5) );
            }
            else {
              if ( i <= blending_opt_part[1] || i >= blending_opt_part[0] )
                xgeFillRectangle ( 3, 3,
                    (int)(rpoints[id][k].x-0.5), (int)(rpoints[id][k].y-0.5) );
            }
          }
        }
    }
    else {
      for ( i = k = 0; i < lastknot_u-degree_u;  i++ )
        for ( j = 0;  j < lastknot_v-degree_v;  j++, k++ ) {
          if ( i >= blending_opt_part[0] && i <= blending_opt_part[1] &&
               j >= blending_opt_part[2] && j <= blending_opt_part[3] &&
               clpoints[id][k] )
            xgeFillRectangle ( 3, 3,
                (int)(rpoints[id][k].x-0.5), (int)(rpoints[id][k].y-0.5) );
        }
    }
  }
  else {
    for ( i = 0;
          i < (lastknot_u-degree_u)*(lastknot_v-degree_v);
          i++ )
      if ( clpoints[id][i] && !(mkpoints[i] & 0x01) )
        xgeFillRectangle ( 3, 3,
            (int)(rpoints[id][i].x-0.5), (int)(rpoints[id][i].y-0.5) );
    xgeSetForeground ( xgec_OrangeRed );
    for ( i = 0;
          i < (lastknot_u-degree_u)*(lastknot_v-degree_v);
          i++ )
      if ( clpoints[id][i] && (mkpoints[i] & 0x01) )
        xgeFillRectangle ( 3, 3,
            (int)(rpoints[id][i].x-0.5), (int)(rpoints[id][i].y-0.5) );
  }
} /*DisplayControlPoints*/

/* this procedure displays the control net subject to the transformation */
/* applied before the blending surface nonlinear optimization */
void DisplayPreTransControlNet ( int id )
{
  void    *sp;
  point3d *tcp, q0, q1;
  int     ncp, i, j, k, pitch, apitch;
  int     i0, i1, j0, j1;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  pitch = lastknot_v-degree_v;
  ncp = (lastknot_u-degree_u)*pitch;
  tcp = pkv_GetScratchMem ( ncp*sizeof(point3d) );
  if ( tcp ) {
    j0 = max ( 0, blending_opt_part[2]-degree_v );
    j1 = min ( lastknot_v-degree_v, blending_opt_part[3]+degree_v+1 );
    apitch = j1-j0;
    if ( kwind.closed_u && blending_opt_part[0] > blending_opt_part[1] ) {
      i0 = blending_opt_part[0]-degree_u;
      i1 = blending_opt_part[1]+degree_u+1;
      for ( i = i0, k = 0;  i < lastknot_u-2*degree_u;  i++ )
        for ( j = j0;  j < j1;  j++, k++ ) {
          Point4to3d ( &cpoints[i*pitch+j], &tcp[k] );
          TransPoint3d ( &blending_opt_transform, &tcp[k], &tcp[k] );
        }
      for ( i = 0; i < i1; i++ )
        for ( j = j0;  j < j1;  j++, k++ ) {
          Point4to3d ( &cpoints[i*pitch+j], &tcp[k] );
          TransPoint3d ( &blending_opt_transform, &tcp[k], &tcp[k] );
        }
      i1 = lastknot_u-2*degree_u+i1-i0;
    }
    else {
      i0 = max ( 0, blending_opt_part[0]-degree_u );
      i1 = min ( lastknot_u-degree_u, blending_opt_part[1]+degree_u+1 );
      for ( i = i0, k = 0;  i < i1;  i++ )
        for ( j = j0;  j < j1;  j++, k++ ) {
          Point4to3d ( &cpoints[i*pitch+j], &tcp[k] );
          TransPoint3d ( &blending_opt_transform, &tcp[k], &tcp[k] );
        }
      i1 -= i0;
    }
    xgeSetForeground ( xgec_NavajoWhite3 );
    for ( i = 0; i < i1; i++ )
      for ( j = 0; j < apitch-1; j++ ) {
        CameraClipLine3d ( &swind.CPos[id], &tcp[i*apitch+j], 0.0,
                           &tcp[i*apitch+j+1], 1.0, &q0, &q1 );
        xgeDrawLine ( (int)(q0.x+0.5), (int)(q0.y+0.5),
                      (int)(q1.x+0.5), (int)(q1.y+0.5) );
      }
    for ( i = 0; i < i1-1; i++ )
      for ( j = 0; j < apitch; j++ ) {
        CameraClipLine3d ( &swind.CPos[id], &tcp[i*apitch+j], 0.0,
                           &tcp[(i+1)*apitch+j], 1.0, &q0, &q1 );
        xgeDrawLine ( (int)(q0.x+0.5), (int)(q0.y+0.5),
                      (int)(q1.x+0.5), (int)(q1.y+0.5) );
      }
    for ( i = 0; i < k; i++ )
      if ( CameraClipPoint3d ( &swind.CPos[id], &tcp[i], &q0 ) )
        xgeFillRectangle ( 3, 3, (int)(q0.x-0.5), (int)(q0.y-0.5) );
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayPreTransControlNet*/

