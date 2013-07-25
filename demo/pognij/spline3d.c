
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

#include "pkvaria.h"
#include "pkgeom.h"
#include "convh.h"
#include "multibs.h"
#include "camera.h"
#include "xgedit.h"

#include "spline3d.h"
#include "ed3dspl.h"

xge_3Dwind   cwind;
xge_KnotWind kwind;

boolean curve   = true;
boolean lamana  = true;
boolean bezpoly = false;
boolean ticks   = false;
boolean convh   = false;
boolean nurbs   = false;
boolean polarf  = false;
boolean curvgr  = false;
boolean torsgr  = false;
boolean uniform = false;
boolean display_coord   = false;
boolean move_many_knots = false;
double  curvgraphscale = 0.0;  /* from 0 to 1 */
double  torsgraphscale = 0.0;
int     curv_graph_dens = 16;

int     npoints;

point4d cpoints[8*MAX_KNOTS];      /* curve control points */
point4d savedcpoints[8*MAX_KNOTS];
byte    cpmark[8*MAX_KNOTS];
point2d rpoints[4][8*MAX_KNOTS];   /* projections */
boolean clpoints[4][8*MAX_KNOTS];  /* true if fits into the frame */
double   knots[8*MAX_KNOTS+1];

/* ///////////////////////////////////////////////////////////////////////// */
static void ErrorHandler ( int module, int errno, const char *errstr )
{
  FILE *f;
  int  i;

  printf ( "Error %d in %d: %s,\n", errno, module, errstr );
  printf ( "  writing file 'pognij.dat'\n" );
  f = fopen ( "pognij.dat", "w+" );
  fprintf ( f, "degree = %d\n", kwind.degree );
  fprintf ( f, "lastknot = %d\n", kwind.lastknot );
  fprintf ( f, "knots:\n" );
  for ( i = 0; i <= kwind.lastknot; i++ )
    fprintf ( f, "%f,\n", knots[i] );
  fprintf ( f, "control points:\n" );
  for ( i = 0; i < kwind.lastknot-kwind.degree; i++ )
    fprintf ( f, "%f,%f,%f,%f\n", cpoints[i].x, cpoints[i].y,
              cpoints[i].z, cpoints[i].w );
  fclose ( f );
  exit ( 1 );
} /*ErrorHandler*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean RzutujPunkt3 ( int id, point3d *p, point2d *q )
{
  point3d b;
  boolean clp;

  id &= 0x03;
  if ( id < 3 ) {
    CameraProjectPoint3d ( &cwind.CPos[id], p, &b );
    clp = true;  /* for now */
  }
  else {
    clp = CameraClipPoint3d ( &cwind.CPos[3], p, &b );
    if ( !clp )
      CameraProjectPoint3d ( &cwind.CPos[3], p, &b );  /* a horrible solution */
  }
  q->x = b.x;  q->y = b.y;
  return clp;
} /*RzutujPunkt3*/

static boolean RzutujPunkt ( int id, point4d *p, point2d *q )
{
  point3d a;

  id &= 0x03;
  Point4to3d ( p, &a );
  return RzutujPunkt3 ( id, &a, q );
} /*RzutujPunkt*/

void ProjectCurve ( int id )
{
  int j;

  id &= 0x03;
  for ( j = 0; j < npoints-1; j++ )
    MidPoint4d ( &cpoints[j], &cpoints[j+1], &cpoints[j+npoints] );
  for ( j = 0; j < 2*npoints-1; j++ )
    clpoints[id][j] = RzutujPunkt ( id, &cpoints[j], &rpoints[id][j] );
} /*ProjectCurve*/

void ResetObject ( void )
{
  int     i, j;
  point3d q, r;

  pkv_SetErrorHandler ( ErrorHandler );

  curve   = true;
  lamana  = true;
  ticks   = false;
  convh   = false;
  nurbs   = false;
  kwind.closed  = false;
  polarf  = false;
  curvgr  = false;

  kwind.lastknot = 3;
  knots[0] = 0.0;
  knots[1] = 0.0;
  knots[2] = 1.0;
  knots[3] = 1.0;
  kwind.degree  = 1;
  SetPoint4d ( &cpoints[0], -1.0, 0.0, 0.0, 1.0 );
  SetPoint4d ( &cpoints[1],  1.0, 0.0, 0.0, 1.0 );
  npoints = kwind.lastknot-kwind.degree;
  ClearPointMarking ();
  for ( i = 0; i < 2; i++ ) {
    Point4to3d ( &cpoints[i], &q );
    for ( j = 0; j < 3; j++ ) {
      CameraProjectPoint3d ( &cwind.CPos[j], &q, &r );
      SetPoint3d ( &r, (double)((int)(r.x+0.5)), (double)((int)(r.y+0.5)), r.z );
      CameraUnProjectPoint3d ( &cwind.CPos[j], &r, &q );
    }
    Point3to4d ( &q, cpoints[i].w, &cpoints[i] );
  }

  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );

  cwind.current_point = 0;
} /*ResetObject*/

void ClearPointMarking ( void )
{
  memset ( cpmark, 0, npoints*sizeof(boolean) );
} /*ClearPointMarking*/

void ResizeObject ( void )
{
  int i;

  xge_KnotWindInitMapping ( &kwind, kwind.umin, kwind.umax );
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
} /*ResizeObject*/

boolean FindNearestPoint ( int id, int x, int y, int mindist )
{
  int   i, m;
  double d, e;

  id &= 0x03;
  if ( nurbs )
    m = 2*npoints-1;
  else m = npoints;
  e = (double)mindist;
  for ( i = 0; i < m; i++ ) {
    if ( clpoints[id][i] ) {
      d = (double)(fabs((double)x-rpoints[id][i].x)+fabs((double)y-rpoints[id][i].y));
      if ( d < e ) {
        cwind.current_point = i;
        e = d;
      }
    }
  }
  return (boolean)(e < mindist);
} /*FindNearestPoint*/

void SetCPoint ( int id, int x, int y )
{
  int      i;
  point3d  q, r;
  vector2d v1, v2;
  double   x1, y1, x2, y2, s, t, w;

  id &= 0x03;
  i = cwind.current_point;
  if ( kwind.closed ) {
    if ( (i < npoints && i >= kwind.clcK) || (nurbs && i >= npoints+kwind.clcK) )
      i -= kwind.clcK;
  }

  if ( nurbs )
    w = cpoints[i].w;
    else w = 1.0;

  if ( i < npoints ) {
    Point4to3d ( &cpoints[i], &q );
    CameraProjectPoint3d ( &cwind.CPos[id], &q, &r );
    r.x = (double)x;  r.y = (double)y;
    CameraUnProjectPoint3d ( &cwind.CPos[id], &r, &q );
    Point3to4d ( &q, cpoints[i].w, &cpoints[i] );
    if ( kwind.closed ) {
      if ( i < kwind.degree )
        cpoints[i+kwind.clcK] = cpoints[i];
      else if ( i >= kwind.clcK )
        cpoints[i-kwind.clcK] = cpoints[i];
    }
  }
  else if ( nurbs ) {
    x1 = rpoints[id][i-npoints].x;    y1 = rpoints[id][i-npoints].y;
    x2 = rpoints[id][i-npoints+1].x;  y2 = rpoints[id][i-npoints+1].y;
    SetVector2d ( &v1, x2-x1, y2-y1 );
    s = (double)DotProduct2d ( &v1, &v1 );
    if ( s >= 1.0 ) {
      SetVector2d ( &v2, (double)x-x1, (double)y-y1 );
      t = (double)DotProduct2d ( &v1, &v2 )/s;
      t = (double)min( t, 0.99 );
      t = (double)max( t, 0.01 );
      i -= npoints;
      w = (double)(cpoints[i].w*(t/(1.0-t)));
      Point4to3d ( &cpoints[i+1], &r );
      Point3to4d ( &r, w, &cpoints[i+1] );
      if ( kwind.closed ) {
        if ( i < kwind.degree )
          cpoints[i+kwind.clcK] = cpoints[i];
        else if ( i >= kwind.clcK )
          cpoints[i-kwind.clcK] = cpoints[i];
      }
    }
  }
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
} /*SetCPoint*/

void SelectCPoints ( int id, short x0, short x1, short y0, short y1 )
{
  int i, k;

  id &= 0x03;
  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < npoints; i++ )
      if ( rpoints[id][i].x >= x0 && rpoints[id][i].x <= x1 &&
           rpoints[id][i].y >= y0 && rpoints[id][i].y <= y1 )
        cpmark[i] = true;
  }
  else {
    if ( FindNearestPoint ( id, x0, y0, xge_MINDIST ) )
      if ( cwind.current_point < npoints ) {
        cpmark[cwind.current_point] = true;
        if ( kwind.closed ) {
          if ( cwind.current_point < kwind.degree )
            k = cwind.current_point+kwind.clcK;
          else if ( cwind.current_point >= kwind.clcK )
            k = cwind.current_point-kwind.clcK;
          else return;
            cpmark[k] = true;
        }
      }
  }
} /*SelectCPoints*/

void UnSelectCPoints ( int id, short x0, short x1, short y0, short y1 )
{
  int i, k;

  id &= 0x03;
  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < npoints; i++ )
      if ( rpoints[id][i].x >= x0 && rpoints[id][i].x <= x1 &&
           rpoints[id][i].y >= y0 && rpoints[id][i].y <= y1 )
        cpmark[i] = false;
  }
  else {
    if ( FindNearestPoint ( id, x0, y0, xge_MINDIST ) )
      if ( cwind.current_point < npoints ) {
        cpmark[cwind.current_point] = false;
        if ( kwind.closed ) {
          if ( cwind.current_point < kwind.degree )
            k = cwind.current_point+kwind.clcK;
          else if ( cwind.current_point >= kwind.clcK )
            k = cwind.current_point-kwind.clcK;
          else return;
            cpmark[k] = false;
        }
      }
  }
} /*UnSelectCPoints*/

void TransformCPoints ( trans3d *tr )
{
  int i;

  for ( i = 0; i < npoints; i++ )
    if ( cpmark[i] )
      Trans3Point4d ( tr, &savedcpoints[i], &cpoints[i] );
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
} /*TransformCPoints*/

void SaveCPoints ( void )
{
  memcpy ( savedcpoints, cpoints, npoints*sizeof(point4d) );
} /*SaveCPoints*/

boolean DegreeElevation ( void )
{
  void    *sp;
  int     k;
  double   *nkn;
  point4d *ncp;

  if ( kwind.degree < MAX_DEGREE ) {
                             /* compute the final number of knots */
    k = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
    if ( k*kwind.degree+2 < MAX_KNOTS ) {
      sp = pkv_GetScratchMemTop ();
      nkn = (double*)pkv_GetScratchMemd ( MAX_KNOTS );
      ncp = (point4d*)pkv_GetScratchMem ( MAX_KNOTS*sizeof(point4d) );
      if ( !nkn || !ncp ) {
        printf ( "Not enough memory\n" );
        exit ( 1 );
      }
      if ( kwind.closed )
        mbs_BSDegElevClosedC4d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                                 &kwind.degree, &kwind.lastknot, nkn, ncp );
      else
        mbs_BSDegElevC4d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                           &kwind.degree, &kwind.lastknot, nkn, ncp, true );
      npoints = kwind.lastknot-kwind.degree;
      memcpy ( knots, nkn, (kwind.lastknot+1)*sizeof(double) );
      memcpy ( cpoints, ncp, npoints*sizeof(point4d) );
      if ( kwind.closed ) {
        kwind.clcK = kwind.lastknot-2*kwind.degree;
        kwind.clcT = knots[kwind.degree+kwind.clcK]-knots[kwind.degree];
      }
      pkv_SetScratchMemTop ( sp );
      memset ( cpmark, 0, npoints*sizeof(boolean) );
      uniform = false;
      for ( k = 0; k < 4; k++ )
        ProjectCurve ( k );
      return true;
    }
    else
      return false;
  }
  else {
    xge_DisplayErrorMessage ( "Error: Cannot raise degree above the limit.", -1 );
    return false;
  }
} /*DegreeElevation*/

boolean DegreeReduction ( void )
{
  void    *sp;
  double   *nkn;
  point4d *ncp;
  int     id;

  sp = pkv_GetScratchMemTop ();
  if ( kwind.degree > 1 ) {
    nkn = pkv_GetScratchMemd ( MAX_KNOTS );
    ncp = pkv_GetScratchMem ( MAX_KNOTS*sizeof(point4d) );
    if ( !nkn || !ncp ) {
      printf ( "Not enough memory\n" );
      exit ( 1 );
    }
    if ( kwind.closed ) {
      if ( !mbs_BSDegRedClosedC4d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                                    &kwind.degree, &kwind.lastknot, nkn, ncp ) )
        goto failure;
      kwind.clcK = kwind.lastknot-2*kwind.degree;
    }
    else {
      if ( !mbs_BSDegRedC4d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                              &kwind.degree, &kwind.lastknot, nkn, ncp ) )
        goto failure;
    }
    npoints = kwind.lastknot-kwind.degree;
    memcpy ( knots, nkn, (kwind.lastknot+1)*sizeof(double) );
    memcpy ( cpoints, ncp, npoints*sizeof(point4d) );
    memset ( cpmark, 0, npoints*sizeof(boolean) );
    pkv_SetScratchMemTop ( sp );
    uniform = false;
    for ( id = 0; id < 4; id++ )
      ProjectCurve ( id );
    return true;
  }

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeReduction*/

void InsertKnot ( void )
{
  int i;

  if ( kwind.closed ) {
    mbs_KnotInsClosedC4d ( kwind.degree, &kwind.lastknot, knots,
                           cpoints, kwind.newknot );
    kwind.clcK = kwind.lastknot - 2*kwind.degree;
  }
  else
    mbs_KnotInsC4d ( kwind.degree, &kwind.lastknot, knots,
                     cpoints, kwind.newknot );
  npoints = kwind.lastknot - kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  uniform = false;
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
} /*InsertKnot*/

void RemoveKnot ( void )
{
  int i;

  if ( kwind.closed ) {
    mbs_KnotRemoveClosedC4d ( kwind.degree, &kwind.lastknot, knots,
                              cpoints, kwind.current_knot );
    kwind.clcK = kwind.lastknot - 2*kwind.degree;
  }
  else
    mbs_KnotRemoveC4d ( kwind.degree, &kwind.lastknot, knots,
                        cpoints, kwind.current_knot );
  npoints = kwind.lastknot - kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  uniform = false;
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
} /*RemoveKnot*/

void SetUniformKnots ( void )
{
  double du, delta;
  int   i, N;

  N = kwind.lastknot;
  du = knots[N-1]-knots[1];
  delta = du/(double)(N-2);
  for ( i = 2; i < N-1; i++ )
    knots[i] = knots[1]+(double)(i-1)*delta;
  if ( kwind.closed ) {
    kwind.clcK = kwind.lastknot-2*kwind.degree;
    kwind.clcT = knots[kwind.clcK+kwind.degree]-knots[kwind.degree];
  }
  uniform = true;
} /*SetUniformKnots*/

static void RzutujPK ( int id, point4d *p, point3d *q )
{
  point4d r;

  id &= 0x03;
  Trans3Point4d ( &cwind.CPos[id].CTr, p, &r );
  if ( id < 3 )
    SetPoint3d ( q, r.x, r.y, r.w );
  else
    SetPoint3d ( q, r.x+cwind.CPos[3].vd.persp.xi0*r.z,
                 r.y+cwind.CPos[3].vd.persp.eta0*r.z, r.z );
} /*RzutujPK*/

/* TODO: clipping of all items to be displayed */
static void GetClipLines ( int id, vector3d cliplines[4] )
{
  id &= 0x03;
  if ( id >= 0 && id < 3 ) {
    SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-cwind.CPos[id].xmin) );
    SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-cwind.CPos[id].ymin) );
    SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(cwind.CPos[id].xmin+cwind.CPos[id].width) );
    SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(cwind.CPos[id].ymin+cwind.CPos[id].height) );
  }
  else {
    SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-cwind.CPos[3].xmin) );
    SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-cwind.CPos[3].ymin) );
    SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(cwind.CPos[3].xmin+cwind.CPos[3].width) );
    SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(cwind.CPos[3].ymin+cwind.CPos[3].height) );
  }
} /*GetClipLines*/

static void DisplayBez ( int id, int n, point4d *p )
{
  void     *sp;
  int      i;
  point3d  *cp;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
  cp = (point3d*)pkv_GetScratchMem ( (n+1)*sizeof(point3d) );
  if ( !cp )
    { printf ( "*" );  exit ( 1 ); }
  GetClipLines ( id, cliplines );
  for ( i = 0; i <= n; i++ )
    RzutujPK ( id, &p[i], &cp[i] );
  mbs_ClipBC2Rd ( 4, cliplines, n, cp, xge_DrawBC2Rd );
  pkv_SetScratchMemTop ( sp );
} /*DisplayBez*/

static void DisplayTick ( int id, point4d *p1, point4d *p2 )
{
  int      x1, y1, x2, y2;
  point2d  a, b;
  vector2d v;

  if ( RzutujPunkt ( id, p1, &a ) && RzutujPunkt ( id, p2, &b ) ) {
    SubtractPoints2d ( &a, &b, &v );
    NormalizeVector2d ( &v );
    MultVector2d ( 4.5, &v, &v );
    x1 = (int)(a.x+v.y+0.5);
    y1 = (int)(a.y-v.x+0.5);
    x2 = (int)(a.x-v.y+0.5);
    y2 = (int)(a.y+v.x+0.5);
    xgeDrawLine ( x1, y1, x2, y2 );
  }
} /*DisplayTick*/

void DisplayCurve ( int id )
{
  int     i, j, kpcs;
  point4d *ncp;
  int     scratchsize;

  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point4d*)pkv_GetScratchMem (
                      scratchsize = (kwind.degree+1)*kpcs*sizeof(point4d) );
  mbs_BSToBezC4d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs, NULL, NULL, ncp );

  xgeSetForeground ( xgec_White );
  for ( i = 0; i < kpcs; i++ )
    DisplayBez ( id, kwind.degree, &ncp[i*(kwind.degree+1)] );
  if ( ticks ) {
    xgeSetLineAttributes ( 2, LineSolid, CapButt, JoinMiter );
    xgeSetForeground ( xgec_DeepPink );
    for ( i = 0; i < kpcs; i++ ) {
      j = i*(kwind.degree+1);
      DisplayTick ( id, &ncp[j], &ncp[j+1] );
    }
    j = (kwind.degree+1)*kpcs-1;
    DisplayTick ( id, &ncp[j], &ncp[j-1] );
    xgeSetLineAttributes ( 1, LineSolid, CapButt, JoinMiter );
  }

  pkv_FreeScratchMem ( scratchsize );
} /*DisplayCurve*/

static void DisplayLine3 ( int id, point3d *p0, point3d *p1 )
{
  point3d r0, r1;

  id &= 0x03;
  if ( id == 3 ) {
    if ( CameraClipLine3d ( &cwind.CPos[3], p0, 0.0, p1, 1.0, &r0, &r1 ) ) {
      xgeDrawLine ( (int)(r0.x+0.5), (int)(r0.y+0.5),
                    (int)(r1.x+0.5), (int)(r1.y+0.5) );
    }
  }
  else {
  }
} /*DisplayLine3*/

static void DisplayPolyline4 ( int id, int np, point4d *p )
{
  void    *sp;
  int     i;
  point2d q;
  point3d p0, p1, q0, q1;
  boolean *clp;
  XPoint  *poly;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  clp = (boolean*)pkv_GetScratchMem ( np );
  poly = (XPoint*)pkv_GetScratchMem ( np*sizeof(XPoint) );
  if ( clp && poly ) {
    for ( i = 0; i < np; i++ ) {
      if ( (clp[i] = RzutujPunkt ( id, &p[i], &q )) ) {
        poly[i].x = (short)(q.x+0.5);
        poly[i].y = (short)(q.y+0.5);
      }
    }
    for ( i = 0; i < np-1; i++ ) {
      if ( clp[i] && clp[i+1] )
        xgeDrawLine ( poly[i].x, poly[i].y, poly[i+1].x, poly[i+1].y );
      else {
        if ( id == 3 ) {
          Point4to3d ( &p[i], &p0 );
          Point4to3d ( &p[i+1], &p1 );
          if ( CameraClipLine3d ( &cwind.CPos[3], &p0, 0.0, &p1, 1.0, &q0, &q1 ) )
            xgeDrawLine ( (int)(q0.x+0.5), (int)(q0.y+0.5),
                          (int)(q1.x+0.5), (int)(q1.y+0.5) );
        }
        else {
        }
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayPolyline4*/

void DisplayCPolygon ( int id )
{
  id &= 0x03;
  xgeSetForeground ( xgec_Green );
  DisplayPolyline4 ( id, npoints, cpoints );
} /*DisplayCPolygon*/

void DisplayCPoints ( int id ) 
{
  int i;
  int x, y;

  id &= 0x03;
  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < npoints; i++ )
    if ( clpoints[id][i] && !cpmark[i] ) {
      x = (int)(rpoints[id][i].x+0.5);
      y = (int)(rpoints[id][i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < npoints; i++ )
    if ( clpoints[id][i] && cpmark[i] ) {
      x = (int)(rpoints[id][i].x+0.5);
      y = (int)(rpoints[id][i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
} /*DisplayCPoints*/

void DisplayBezPoly ( int id )
{
  void    *sp;
  int     i, kpcs;
  point4d *ncp;
  boolean *clp;
  point2d q;
  XPoint  *poly;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );

  ncp = (point4d*)pkv_GetScratchMem ( (kwind.degree+1)*kpcs*sizeof(point4d) );
  clp = (boolean*)pkv_GetScratchMem ( kwind.degree+1 );
  poly = (XPoint*)pkv_GetScratchMem ( (kwind.degree+1)*sizeof(XPoint) );

  if ( ncp && clp && poly ) {
    mbs_BSToBezC4d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs, NULL, NULL, ncp );

    xgeSetForeground ( xgec_Green );
    for ( i = 0; i < kpcs; i++ )
      DisplayPolyline4 ( id, kwind.degree+1, &ncp[i*(kwind.degree+1)] );

    xgeSetForeground ( xgec_CornflowerBlue );
    for ( i = 0; i < kpcs*(kwind.degree+1); i++ ) {
      if ( RzutujPunkt ( id, &ncp[i], &q ) )
        xgeFillRectangle ( 3, 3, (int)(q.x+0.5)-1, (int)(q.y+0.5)-1 );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayBezPoly*/

void DisplayAuxPoints ( int id )
{
  int i;
  int x, y;

  id &= 0x03;
  xgeSetForeground ( xgec_HotPink );
  for ( i = npoints; i < 2*npoints-1; i++ )
    if ( clpoints[id][i] ) {
      x = (int)(rpoints[id][i].x+0.5);
      y = (int)(rpoints[id][i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
} /*DisplayAuxPoints*/

void DisplayConvh ( int id )
{
  int     i, j, k, kpcs;
  point4d *ncp;
  int     scratchsize;
  point2d chf[MAX_DEGREE+2];
  XPoint  ch[MAX_DEGREE+3];

  id &= 0x03;
  xgeSetForeground ( xgec_DimGrey );
  if ( bezpoly ) {
    kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
    ncp = (point4d*)pkv_GetScratchMem (
                      scratchsize = (kwind.degree+1)*kpcs*sizeof(point4d) );
    mbs_BSToBezC4d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs, NULL, NULL, ncp );
    for ( i = 0; i < kpcs; i++ ) {
      for ( j = 0; j <= kwind.degree; j++ )
        RzutujPunkt ( id, &ncp[i*(kwind.degree+1)+j], &chf[j] );
      j = kwind.degree+1;
      FindConvexHull2d ( &j, chf );
      if ( j > 2 ) {
        for ( k = 0; k < j; k++ ) {
          ch[k].x = (short)(chf[k].x+0.5);
          ch[k].y = (short)(chf[k].y+0.5);
        }
        xgeFillPolygon ( Convex, j, ch );
      }
    }
    xgeSetForeground ( xgec_Grey );
    for ( i = 0; i < kpcs; i++ ) {
      for ( j = 0; j <= kwind.degree; j++ )
        RzutujPunkt ( id, &ncp[i*(kwind.degree+1)+j], &chf[j] );
      j = kwind.degree+1;
      FindConvexHull2d ( &j, chf );
      if ( j > 1 ) {
        for ( k = 0; k < j; k++ ) {
          ch[k].x = (short)(chf[k].x+0.5);
          ch[k].y = (short)(chf[k].y+0.5);
        }
        ch[j] = ch[0];
        xgeDrawLines ( j+1, ch );
      }
    }
    pkv_FreeScratchMem ( scratchsize );
  }
  else {
    for ( i = 0; i < kwind.lastknot-2*kwind.degree; i++ ) {
      if ( knots[i+kwind.degree] < knots[i+kwind.degree+1] ) {
        for ( j = 0; j <= kwind.degree; j++ )
          RzutujPunkt ( id, &cpoints[i+j], &chf[j] );
        j = kwind.degree+1;
        FindConvexHull2d ( &j, chf );
        if ( j > 2 ) {
          for ( k = 0; k < j; k++ ) {
            ch[k].x = (short)(chf[k].x+0.5);
            ch[k].y = (short)(chf[k].y+0.5);
          }
          xgeFillPolygon ( Convex, j, ch );
        }
      }
    }
    xgeSetForeground ( xgec_Grey );
    for ( i = 0; i <= kwind.lastknot-2*kwind.degree-1; i++ ) {
      if ( knots[i+kwind.degree] < knots[i+kwind.degree+1] ) {
        for ( j = 0; j <= kwind.degree; j++ )
          RzutujPunkt ( id, &cpoints[i+j], &chf[j] );
        j = kwind.degree+1;
        FindConvexHull2d ( &j, chf );
        if ( j > 2 ) {
          for ( k = 0; k < j; k++ ) {
            ch[k].x = (short)(chf[k].x+0.5);
            ch[k].y = (short)(chf[k].y+0.5);
          }
          ch[j] = ch[0];
          xgeDrawLines ( j+1, ch );
        }
      }
    }
  }
} /*DisplayConvh*/

static void GetPFKnots ( int i, int *rr, double *uu, double *pfknots )
{
  int j, r, dr;
  double u;

  *uu = u = knots[i];
  r = 0;
  while ( knots[i+r] <= u )
    r++;
  *rr = r;
  dr = kwind.degree-r;
  for ( j = 0; j <= dr; j++ ) {
    pfknots[j] = knots[i-(dr+1)+j];
    pfknots[dr+1+j] = knots[i+r+j];
  }
} /*GetPFKnots*/

void DisplayPolarf ( int id )
{
  void     *sp;
  int      i, r, dr, nf, kpcs;
  double    u;
  double    *pfknots;
  point4d  *bcpf;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  if ( !(bcpf = pkv_GetScratchMem (3*(kwind.degree+1)*sizeof(point4d)) ) )
    goto failure;
  if ( !(pfknots = pkv_GetScratchMemd ( 2*(kwind.degree+1) )) )
    goto failure;

  xgeSetForeground ( xgec_Cyan );
  if ( kwind.closed )
    nf = kwind.lastknot-kwind.degree+1;
  else
    nf = kwind.lastknot-kwind.degree;
  for ( i = kwind.degree+1; i < nf; i += r ) {
    GetPFKnots ( i, &r, &u, pfknots );
    if ( r < kwind.degree ) {
      dr = kwind.degree-r;
      mbs_BSToBezC4d ( dr, 2*dr+1, pfknots, &cpoints[i-dr-1],
                       &kpcs, NULL, NULL, bcpf );
      DisplayBez ( id, dr, bcpf );
    }
  }
failure:
  pkv_SetScratchMemTop ( sp );
} /*DisplayPolarf*/

static void DisplayBezHedgehog ( int id, int n, point4d *cp )
{
  double    scale_c, scale_t;
  int      i;
  double    t;
  double    curvature[2];
  point3d  cpoint, a, b;
  vector3d fframe[3];

  id &= 0x03;
  if ( curvgr ) {
    xgeSetForeground ( xgec_CornflowerBlue );
    scale_c = (double)exp(-4.6051702+6.9077553*curvgraphscale); /* from 0.01 to 10.0 */
    for ( i = 0; i <=curv_graph_dens; i++ ) {
      t = (double)i/(double)curv_graph_dens;
      mbs_BCFrenetC3Rd ( n, cp, t, &cpoint, fframe, curvature );
      AddVector3Md ( &cpoint, &fframe[1], scale_c*curvature[0], &a );
      DisplayLine3 ( id, &cpoint, &a );
    }
  }
  if ( torsgr ) {
    xgeSetForeground ( xgec_SpringGreen2 );
    scale_t = (double)exp(-4.6051702+6.9077553*torsgraphscale);
    for ( i = 0; i <= curv_graph_dens; i++ ) {
      t = (double)i/(double)curv_graph_dens;
      mbs_BCFrenetC3Rd ( n, cp, t, &cpoint, fframe, curvature );
      AddVector3Md ( &cpoint, &fframe[2], scale_t*curvature[1], &b );
      DisplayLine3 ( id, &cpoint, &b );
    }
  }
} /*DisplayBezHedgehog*/

void DisplayCurvGraph ( int id )
{
  int i, kpcs;   
  point4d *ncp;
  int scratchsize;  

  id &= 0x03;
  if ( kwind.degree < 2 ) 
    return;
  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point4d*)pkv_GetScratchMem (
                    scratchsize = (kwind.degree+1)*kpcs*sizeof(point4d) );
  mbs_BSToBezC4d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs, NULL, NULL, ncp );


  for ( i = 0; i < kpcs; i++ )
    DisplayBezHedgehog ( id, kwind.degree, &ncp[i*(kwind.degree+1)] );

  pkv_FreeScratchMem ( scratchsize );
} /*DisplayCurvGraph*/

void FindBoundingBox ( Box3d *box )
{
  int i;
  double xyz;

  box->x0 = box->x1 = cpoints[0].x/cpoints[0].w;
  box->y0 = box->y1 = cpoints[0].y/cpoints[0].w;
  box->z0 = box->z1 = cpoints[0].z/cpoints[0].w;
  for ( i = 1; i < kwind.lastknot-kwind.degree; i++ ) {
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
} /*FindBoundingBox*/

void UstawNURBS ( void )
{
  int i;
  point3d p;

  if ( !nurbs ) {
    for ( i = 0; i < npoints; i++ ) {
      Point4to3d ( &cpoints[i], &p );
      Point3to4d ( &p, 1.0, &cpoints[i] );
    }
    for ( i = 0; i < 4; i++ )
      ProjectCurve ( i );         
  }
} /*UstawNURBS*/

void UstawPolarf ( void )
{
} /*UstawPolarf*/

void UstawCurvGraph ( void )
{
} /*UstawCurvGraph*/

void UstawTorsGraph ( void )
{
} /*UstawTorsGraph*/

static void SetClosedCurve ( void )
{
  int i;

  kwind.clcK = kwind.lastknot-2*kwind.degree;
  kwind.clcT = knots[kwind.clcK+kwind.degree]-knots[kwind.degree];
  for ( i = 1; i < 2*kwind.degree; i++ )
    knots[i] = (double)(0.5*(knots[i]+knots[kwind.clcK+i]-kwind.clcT));
  for ( i = 1; i < 2*kwind.degree; i++ )
    knots[kwind.clcK+i] = knots[i]+kwind.clcT;
  knots[0] = min (knots[1], knots[kwind.clcK]-kwind.clcT );
  knots[kwind.lastknot] = max ( knots[kwind.lastknot-1],
                                knots[kwind.lastknot-kwind.clcK]+kwind.clcT );
  for ( i = 0; i < kwind.degree; i++ ) {
    MidPoint4d ( &cpoints[kwind.clcK+i], &cpoints[i], &cpoints[i] );
    cpoints[kwind.clcK+i] = cpoints[i];
  }
} /*SetClosedCurve*/

boolean UstawZamknieta ( void )
{
  int i;

  if ( kwind.closed ) {
    if ( kwind.lastknot <= 3*kwind.degree ) {
      xge_DisplayErrorMessage ( "Error: Not enough knots to close the curve.", -1 );
      kwind.closed = false;
      return false;
    }
    SetClosedCurve ();
  }
  xge_KnotWindInitMapping ( &kwind, (double)min(0.0, kwind.knots[1]),
                            (double)max(1.0, knots[kwind.lastknot-1]) );
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
  return true;
} /*UstawZamknieta*/

boolean RefineUniform ( void )
{
  void   *sp;
  double *acp;
  int    lkn, i;

  if ( !uniform || 2*(kwind.lastknot-kwind.degree) >= MAX_KNOTS )
    return false;
  sp = pkv_GetScratchMemTop ();
  acp = pkv_GetScratchMem ( (2*kwind.lastknot-3*kwind.degree)*sizeof(point4d) );
  if ( !acp )
    goto failure;
  if ( !mbs_LaneRiesenfeldC4d ( kwind.degree, kwind.lastknot, cpoints,
                                &lkn, acp ) )
    goto failure;
  memcpy ( cpoints, acp, (lkn-kwind.degree)*sizeof(point4d) );
  knots[lkn] = knots[lkn-1] = knots[kwind.lastknot-1];
  kwind.lastknot = lkn;
  SetUniformKnots ();
  npoints = kwind.lastknot-kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  for ( i = 0; i < 4; i++ )
    ProjectCurve ( i );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RefineUniform*/

