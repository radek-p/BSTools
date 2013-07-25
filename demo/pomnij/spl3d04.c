
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
boolean RzutujPunkt3 ( int id, point3d *p, point2d *q )
{
  point3d b;
  boolean clp;

  id &= 0x03;
  if ( id < 3 ) {
    CameraProjectPoint3d ( &swind.CPos[id], p, &b );
    clp = true;  /* for now */
  }
  else {
    clp = CameraClipPoint3d ( &swind.CPos[3], p, &b );
    if ( !clp )
      CameraProjectPoint3d ( &swind.CPos[3], p, &b );  /* a horrible solution */
  }
  SetPoint2d ( q, b.x, b.y );
  return clp;
} /*RzutujPunkt3*/

boolean RzutujPunkt ( int id, point4d *p, point2d *q )
{
  point3d a;

  id &= 0x03;
  Point4to3d ( p, &a );
  return RzutujPunkt3 ( id, &a, q );
} /*RzutujPunkt*/

void ProjectSurface ( int id )
{
  int i;

  id &= 0x03;
  for ( i = 0; i < (lastknot_u-degree_u)*(lastknot_v-degree_v); i++ )
    clpoints[id][i] = RzutujPunkt ( id, &cpoints[i], &rpoints[id][i] );
} /*ProjectSurface*/

void ResizeObject ( void )
{
  int id;

  for ( id = 0; id < 4; id++ )
    ProjectSurface ( id );
} /*ResizeObject*/

void ResetObject ( void )
{
  pkv_SetErrorHandler ( ErrorHandler );
  kwind.closed_u = kwind.closed_v = false;
  xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
  degree_u = degree_v = kwind.degree_u = kwind.degree_v = 1;
  kwind.altknu = kwind.altknv = false;
  knots_u[0] = knots_u[1] = 0.0;
  knots_u[2] = knots_u[3] = 1.0;
  knots_v[0] = knots_v[1] = 0.0;
  knots_v[2] = knots_v[3] = 1.0;
  lastknot_u = lastknot_v = kwind.lastknot_u = kwind.lastknot_v = 3;
  kwind.altlastknot_u = kwind.altlastknot_v = 0;
  kwind.altdeg_u = kwind.altdeg_v = 0;
  SetPoint4d ( &cpoints[0], -1.0, -1.0, 0.0, 1.0 );
  SetPoint4d ( &cpoints[1], -1.0,  1.0, 0.0, 1.0 );
  SetPoint4d ( &cpoints[2],  1.0, -1.0, 0.0, 1.0 );
  SetPoint4d ( &cpoints[3],  1.0,  1.0, 0.0, 1.0 );
  xge_T2KnotWindFindMapping ( &kwind );
  ClearPointMarking ( 4, mkpoints );
  ResetEquator ();
  ResetMeridian ();
  SelectEquator ();
  bind_spr = false;
  sw_bind_blending = sw_nonlin_blending = false;
  blending_mat_valid = false;
  n_blending_constraints = 0;
  FindBoundingBox ( &swind.RefBBox );
  swind.PerspBBox = swind.RefBBox;
  xge_3DwindSetupParProj ( &swind, &swind.RefBBox );
  xge_3DwindSetupPerspProj ( &swind, false );
  ResizeObject ();
  IdentTrans3d ( &blending_opt_transform );
} /*ResetObject*/

void FlipPatch ( void )
{
#define SWAP(a,b,c) { c = a;  a = b;  b = c; }
  void    *sp;
  int     ncp, s;
  double  *buf;
  boolean cl;

        /* swap the "u" and "v" parameters */
  sp = pkv_GetScratchMemTop ();
  ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
  buf = pkv_GetScratchMemd ( max ( 4*ncp, (MAX_DEGREE+1)*MAX_KNOTS ) );
  if ( buf ) {
        /* Flip the patch */
    pkv_TransposeMatrixc ( lastknot_u-degree_u, lastknot_v-degree_v,
                           sizeof(point4d),
                           (lastknot_v-degree_v)*sizeof(point4d), (char*)cpoints,
                           (lastknot_u-degree_u)*sizeof(point4d), (char*)buf );
    memcpy ( cpoints, buf, ncp*sizeof(point4d) );
    memcpy ( buf, knots_u, (lastknot_u+1)*sizeof(double) );
    memcpy ( knots_u, knots_v, (lastknot_v+1)*sizeof(double) );
    memcpy ( knots_v, buf, (lastknot_u+1)*sizeof(double) );
    SWAP ( lastknot_u, lastknot_v, s )
    SWAP ( degree_u, degree_v, s );
        /* setup it in the knot window */
    xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
    kwind.lastknot_u = lastknot_u;
    kwind.degree_u = degree_u;
    kwind.lastknot_v = lastknot_v;
    kwind.degree_v = degree_v;
    SWAP ( kwind.closed_u, kwind.closed_v, cl );
    xge_T2KnotWindFindMapping ( &kwind );
        /* turn off the switches involved */
    bind_spr = false;

    sw_bind_blending = blending_mat_valid = sw_triharmonic_blending =
    sw_nonlin_blending = false;
        /* swap the optimization range bounds */
    SWAP ( blending_opt_part[0], blending_opt_part[2], s );
    SWAP ( blending_opt_part[1], blending_opt_part[3], s );
  }
  pkv_SetScratchMemTop ( sp );
#undef SWAP
} /*FlipPatch*/

void FindBoundingBox ( Box3d *box )
{
#define EPS 1.0e-5
  int i;
  double xyz, c;

  box->x0 = box->x1 = cpoints[0].x/cpoints[0].w;
  box->y0 = box->y1 = cpoints[0].y/cpoints[0].w;
  box->z0 = box->z1 = cpoints[0].z/cpoints[0].w;
  for ( i = 1;
        i < (lastknot_u-degree_u)*(lastknot_v-degree_v);
        i++ ) {
    xyz = cpoints[i].x/cpoints[i].w;
    if ( xyz < box->x0 )      box->x0 = xyz;
    else if ( xyz > box->x1 ) box->x1 = xyz;
    xyz = cpoints[i].y/cpoints[i].w;
    if ( xyz < box->y0 )      box->y0 = xyz;
    else if ( xyz > box->y1 ) box->y1 = xyz;
    xyz = cpoints[i].z/cpoints[i].w;
    if ( xyz < box->z0 )      box->z0 = xyz;
    else if ( xyz > box->z1 ) box->z1 = xyz;
  }
      /* now extend it by 10% and by EPS */
  c = (double)(0.5*(box->x0+box->x1));  xyz = (double)(0.55*(box->x1-box->x0))+EPS;
  box->x0 = c-xyz;  box->x1 = c+xyz;
  c = (double)(0.5*(box->y0+box->y1));  xyz = (double)(0.55*(box->y1-box->y0))+EPS;
  box->y0 = c-xyz;  box->y1 = c+xyz;
  c = (double)(0.5*(box->z0+box->z1));  xyz = (double)(0.55*(box->z1-box->z0))+EPS;
  box->z0 = c-xyz;  box->z1 = c+xyz;
#undef EPS
} /*FindBoundingBox*/

