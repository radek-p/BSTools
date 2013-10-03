
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
#include "spline2d.h"

#define KNOT_EPS 1.0e-4 /* the smallest distance between neighbouring knots */

boolean curve   = true;
boolean lamana  = true;
boolean bezpoly = false;
boolean funkcja = false;
boolean uniform = false;
boolean ticks   = false;
boolean convh   = false;
boolean nurbs   = false;
boolean baza    = false;
boolean polarf  = false;
boolean curvgr  = false;
double  curvgraphscale = 0.0;  /* from 0 to 1 */
int     curv_graph_dens = 16;

xge_2Dwind   cwind;
xge_KnotWind kwind;

int npoints;

point3d cpoints[MAX_DEGREE*MAX_KNOTS];
point3d savedcpoints[MAX_DEGREE*MAX_KNOTS];
byte    cpmark[MAX_DEGREE*MAX_KNOTS];
point2d rpoints[MAX_DEGREE*MAX_KNOTS];
double  knots[MAX_DEGREE*MAX_KNOTS+1];

/* ///////////////////////////////////////////////////////////////////////// */
void DumpData ( void )
{
  FILE *f;
  int  i;

  printf ( "  writing file 'pokrzyw.dat'\n" );
  f = fopen ( "pokrzyw.dat", "w+" );
  fprintf ( f, "degree = %d\n", kwind.degree );
  fprintf ( f, "lastknot = %d\n", kwind.lastknot );
  fprintf ( f, "knots:\n" );
  for ( i = 0; i <= kwind.lastknot; i++ )
    fprintf ( f, "%f,\n", knots[i] );
  fprintf ( f, "control points:\n\n" );
  for ( i = 0; i < kwind.lastknot-kwind.degree; i++ )
    fprintf ( f, "%f,%f,%f\n", cpoints[i].x, cpoints[i].y, cpoints[i].z );
  fclose ( f );
} /*DumpData*/

static void ErrorHandler ( int module, const char *file, int line,
                           int errcode, const char *errstr )
{
  printf ( "Error in module %d, file %s, line %d: %s\n",
           module, file, line, errstr );
  DumpData ();
  exit ( 1 );
} /*ErrorHandler*/

/* ///////////////////////////////////////////////////////////////////////// */
static void RzutujPunkt ( point3d *p, point2d *q )
{
  point2d a ;

  Point3to2d ( p, &a );
  CameraProjectPoint2d ( &cwind.CPos, &a, q );
} /*RzutujPunkt*/

static void RzutujRPunkt ( point3d *p, point3d *q )
{
  point2d ap;

  RzutujPunkt ( p, &ap );
  Point2to3d ( &ap, p->z, q );
} /*RzutujRPunkt*/

void ProjectCurve ( void )
{
  int i;

  for ( i = 0; i < npoints; i++ )
    RzutujPunkt ( &cpoints[i], &rpoints[i] );
  for ( i = 0; i < npoints-1; i++ ) {
    MidPoint3d ( &cpoints[i], &cpoints[i+1], &cpoints[i+npoints] );
    RzutujPunkt ( &cpoints[i+npoints], &rpoints[i+npoints] );
  }
} /*ProjectCurve*/

boolean FindNearestPoint ( int x, int y, int mindist )
{
  int i, m;
  double d, e;

  if ( nurbs )
    m = 2*npoints-1;
  else m = npoints;
  e = (double)mindist;
  for ( i = 0; i < m; i++ ) {
    d = (double)(fabs((double)x-rpoints[i].x) + fabs((double)y-rpoints[i].y));
    if ( d < e ) {
      cwind.current_point = i;
      e = d;
    }
  }
  return (boolean)(e < mindist);
} /*FindNearestPoint*/

void SetCPoint ( int x, int y )
{
  int i;
  double w, r, t;
  double x1, y1, x2, y2;
  vector2d v1, v2;
  point3d q, s;

  i = cwind.current_point;
  if ( kwind.closed ) {
    if ( (i < npoints && i >= kwind.clcK) || (nurbs && i >= npoints+kwind.clcK) )
      i -= kwind.clcK;
  }

  if ( nurbs )
    w = cpoints[i].z;
    else w = 1.0;
  SetPoint3d ( &q, (double)x, (double)y, 1.0 );
  CameraUnProjectPoint3d ( &cwind.CPos, &q, &s );
  MultVector3d ( w, &s, &s );

  if ( funkcja ) {
    cpoints[i].y = s.y;
    if ( kwind.closed ) {
      if ( i < kwind.degree )
        cpoints[i+kwind.clcK].y = cpoints[i].y;
      else if ( i >= kwind.clcK )
        cpoints[i-kwind.clcK].y = cpoints[i].y;
    }
  }
  else if ( i < npoints ) {
    cpoints[i] = s;
    if ( kwind.closed ) {
      if ( i < kwind.degree )
        cpoints[i+kwind.clcK] = cpoints[i];
      else if ( i >= kwind.clcK )
        cpoints[i-kwind.clcK] = cpoints[i];
    }
  }
  else { /* moving a Farin point to change weight */
    x1 = rpoints[i-npoints].x;    y1 = rpoints[i-npoints].y;
    x2 = rpoints[i-npoints+1].x;  y2 = rpoints[i-npoints+1].y;
    SetVector2d ( &v1, x2-x1, y2-y1 );
    r = (double)DotProduct2d ( &v1, &v1 );
    if ( r >= 1 ) {
      SetVector2d ( &v2, (double)x-x1, (double)y-y1 );
      t = (double)DotProduct2d ( &v1, &v2 )/r;
      t = (double)max ( t, 0.01 );  t = (double)min ( t, 0.99 );
      i -= npoints;
      w = (double)(cpoints[i].z*(t/(1.0-t)));   /* the new weight */
      w = (double)min ( w, 10000.0 );  w = (double)max ( w, 0.0001 );
      i++;
      Point2to3d ( &rpoints[i], 1.0, &q );
      CameraUnProjectPoint3d ( &cwind.CPos, &q, &s );
      MultVector3d ( w, &s, &cpoints[i] );
      if ( kwind.closed ) {
        if ( i < kwind.degree )
          cpoints[i+kwind.clcK] = cpoints[i];
        else if ( i >= kwind.clcK )
          cpoints[i-kwind.clcK] = cpoints[i];
      }
    }
  }
  ProjectCurve ();
} /*SetCPoint*/

void SelectCPoints ( short x0, short x1, short y0, short y1 )
{
  int i, k;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < npoints; i++ )
      if ( rpoints[i].x >= x0 && rpoints[i].x <= x1 &&
           rpoints[i].y >= y0 && rpoints[i].y <= y1 )
        cpmark[i] = true;
  }
  else {
    if ( FindNearestPoint ( x0, y0, xge_MINDIST ) )
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

void UnSelectCPoints ( short x0, short x1, short y0, short y1 )
{
  int i, k;

  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 0; i < npoints; i++ )
      if ( rpoints[i].x >= x0 && rpoints[i].x <= x1 &&
           rpoints[i].y >= y0 && rpoints[i].y <= y1 )
        cpmark[i] = false;
  }
  else {
    if ( FindNearestPoint ( x0, y0, xge_MINDIST ) )
      if ( cwind.current_point < npoints )
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
} /*UnSelectCPoints*/

void TransformCPoints ( trans2d *tr )
{
  int i;

  for ( i = 0; i < npoints; i++ )
    if ( cpmark[i] )
      Trans2Point3d ( tr, &savedcpoints[i], &cpoints[i] );
  if ( funkcja )
    UstawFunkcje ();
  ProjectCurve ();
} /*TransformCPoints*/

void SaveCPoints ( void )
{
  memcpy ( savedcpoints, cpoints, npoints*sizeof(point3d) );
} /*SaveCPoints*/

static void DisplayTick ( point3d *p1, point3d *p2 )
{
  vector2d v;
  double x, y;
  int x1, y1, x2, y2;

  if ( funkcja ) {
    x1 = x2 = (int)(p1->x/p1->z);
    y1 = (int)(p1->y/p1->z)-4;
    y2 = y1+9;
  }
  else {
    x = p1->x/p1->z;
    y = p1->y/p1->z;
    SetVector2d ( &v, p2->x/p2->z-x, p2->y/p2->z-y );
    NormalizeVector2d ( &v );
    MultVector2d ( 4.5, &v, &v );
    x1 = (int)(x+v.y);
    y1 = (int)(y-v.x);
    x2 = (int)(x-v.y);
    y2 = (int)(y+v.x);
  }
  xgeDrawLine ( x1, y1, x2, y2 );
} /*DisplayTick*/

static void GetClipLines ( vector3d cliplines[4] )
{
  SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-cwind.CPos.xmin) );
  SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-cwind.CPos.ymin) );
  SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(cwind.CPos.xmin+cwind.CPos.width) );
  SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(cwind.CPos.ymin+cwind.CPos.height) );
} /*GetClipLines*/

void DisplayCurve ( void )
{
  void     *sp;
  int      i, j, kpcs;
  point3d  *ncp;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
    /* Get the clipping line equation coefficients */
  GetClipLines ( cliplines );
    /* divide the spline curve into rational Bezier arcs */
  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point3d*)pkv_GetScratchMem ( (kwind.degree+1)*kpcs*sizeof(point3d) );
  if ( !ncp )
    goto failure;
  mbs_BSToBezC3d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs, NULL, NULL, ncp );
  for ( i = 0; i < kpcs*(kwind.degree+1); i++ )
    RzutujRPunkt ( &ncp[i], &ncp[i] );

    /* clip and draw the arcs */
  xgeSetForeground ( xgec_White );
  for ( i = 0; i < kpcs; i++ ) {
    mbs_ClipBC2Rd ( 4, cliplines, kwind.degree, &ncp[i*(kwind.degree+1)],
                    xge_DrawBC2Rd );
  }

  if ( ticks ) {
    xgeSetLineAttributes ( 2, LineSolid, CapButt, JoinMiter );
    xgeSetForeground ( xgec_DeepPink );
    for ( i = 0; i < kpcs; i++ ) {
      j = i*(kwind.degree+1);
      DisplayTick ( &ncp[j], &ncp[j+1] );
    }
    j = (kwind.degree+1)*kpcs-1;
    DisplayTick ( &ncp[j], &ncp[j-1] );
    xgeSetLineAttributes ( 1, LineSolid, CapButt, JoinMiter );
  }

failure:
  pkv_SetScratchMemTop ( sp );
} /*DisplayCurve*/

void DisplayCPolygon ( void )
{
  int i;
  XPoint poly[MAX_KNOTS];
/*  char s[6]; */

  for ( i = 0; i < npoints; i++ ) {
    poly[i].x = (short)(rpoints[i].x+0.5);
    poly[i].y = (short)(rpoints[i].y+0.5);
  }
  xgeSetForeground ( xgec_Green2 );
  xgeDrawLines ( npoints, poly );
/*
  xgeSetForeground ( xgec_Green );
  for ( i = 0; i < npoints; i++ ) {
    sprintf ( s, "%d", i );
    xgeDrawString ( s, poly[i].x, poly[i].y );
  }
*/
} /*DisplayCPolygon*/

void DisplayCPoints ( void ) 
{
  int i;
  int x, y;

  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < npoints; i++ )
    if ( !cpmark[i] ) {
      x = (int)(rpoints[i].x+0.5);
      y = (int)(rpoints[i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < npoints; i++ )
    if ( cpmark[i] ) {
      x = (int)(rpoints[i].x+0.5);
      y = (int)(rpoints[i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
} /*DisplayCPoints*/

void DisplayBezPoly ( void )
{
  int     i, j, kpcs;
  point3d *ncp;
  point2d q;
  int     scratchsize;
  XPoint  poly[MAX_KNOTS];

  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point3d*)pkv_GetScratchMem (
                    scratchsize = (kwind.degree+1)*kpcs*sizeof(point3d) );
  mbs_BSToBezC3d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs,
                   NULL, NULL, ncp );

  for ( i = 0; i < kpcs; i++ ) {
    for ( j = 0; j <= kwind.degree; j++ ) {
      RzutujPunkt ( &ncp[i*(kwind.degree+1)+j], &q );
      poly[j].x = (short)(q.x+0.5);
      poly[j].y = (short)(q.y+0.5);
    }
    xgeSetForeground ( xgec_LawnGreen );
    xgeDrawLines ( kwind.degree+1, poly );
    xgeSetForeground ( xgec_DodgerBlue );
    for ( j = 0; j <= kwind.degree; j++ )
      xgeFillRectangle ( 3, 3, poly[j].x-1, poly[j].y-1 );
  }

  pkv_FreeScratchMem ( scratchsize );
} /*DisplayBezPoly*/

void DisplayAuxPoints ( void )
{
  int i;
  int x, y;

  xgeSetForeground ( xgec_HotPink );
  for ( i = npoints; i < 2*npoints-1; i++ ) {
    x = (int)(rpoints[i].x);
    y = (int)(rpoints[i].y);
    xgeFillRectangle ( 3, 3, x-1, y-1 );
  }
} /*DisplayAuxPoints*/

void DisplayConvh ( void )
{
  void    *sp;
  int     i, j, k, kpcs;
  point3d *ncp;
  point2d q;
  int     scratchsize;
  point2d chf[MAX_DEGREE+2];
  XPoint  ch[MAX_DEGREE+3];

  sp = pkv_GetScratchMemTop ();
  xgeSetForeground ( xgec_Grey5 );
  if ( bezpoly ) {
    kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
    ncp = (point3d*)pkv_GetScratchMem (
                      scratchsize = (kwind.degree+1)*kpcs*sizeof(point3d) );
    if ( !ncp )
      goto way_out;
    mbs_BSToBezC3d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs,
                     NULL, NULL, ncp );
    for ( i = 0; i < kpcs; i++ ) {
      for ( j = 0; j <= kwind.degree; j++ ) {
        RzutujPunkt ( &ncp[i*(kwind.degree+1)+j], &q );
        chf[j].x = q.x;
        chf[j].y = q.y;
      }
      j = kwind.degree+1;
      if ( !FindConvexHull2d ( &j, chf ) )
        goto way_out;
      if ( j > 2 ) {
        for ( k = 0; k < j; k++ ) {
          ch[k].x = (short)(chf[k].x+0.5);
          ch[k].y = (short)(chf[k].y+0.5);
        }
        xgeFillPolygon ( Convex, j, ch );
      }
    }
    xgeSetForeground ( xgec_Grey3 );
    for ( i = 0; i < kpcs; i++ ) {
      for ( j = 0; j <= kwind.degree; j++ ) {
        RzutujPunkt ( &ncp[i*(kwind.degree+1)+j], &q );
        chf[j].x = (short)(q.x+0.5);
        chf[j].y = (short)(q.y+0.5);
      }
      j = kwind.degree+1;
      if ( !FindConvexHull2d ( &j, chf ) )
        goto way_out;
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
          chf[j] = rpoints[i+j];
        j = kwind.degree+1;
        if ( !FindConvexHull2d ( &j, chf ) )
          goto way_out;
        if ( j > 2 ) {
          for ( k = 0; k < j; k++ ) {
            ch[k].x = (short)(chf[k].x+0.5);
            ch[k].y = (short)(chf[k].y+0.5);
          }
          xgeFillPolygon ( Convex, j, ch );
        }
      }
    }
    xgeSetForeground ( xgec_Grey3 );
    for ( i = 0; i <= kwind.lastknot-2*kwind.degree-1; i++ ) {
      if ( knots[i+kwind.degree] < knots[i+kwind.degree+1] ) {
        for ( j = 0; j <= kwind.degree; j++ )
          chf[j] = rpoints[i+j];
        j = kwind.degree+1;
        if ( !FindConvexHull2d ( &j, chf ) )
          goto way_out;
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
way_out:
  pkv_SetScratchMemTop ( sp );
} /*DisplayConvh*/

void DisplayBasis ( int y, int yu )
{
#define dd 300
  point2d *fgraph;
  point2d cpoint;
  XPoint *graph;
  int scratchsize;
  int i, j;
  double xi, s, t, xn, xNn;

  scratchsize = (kwind.lastknot-kwind.degree+1)*sizeof(point2d) +
                (dd+1)*sizeof(XPoint);
  fgraph = (point2d*)pkv_GetScratchMem ( scratchsize );
  if ( !fgraph )
    return;
  graph = (XPoint*)&fgraph[kwind.lastknot+1-kwind.degree];
  for ( i = 0; i < kwind.lastknot-kwind.degree; i++ ) {
    xi = knots[i+1];
    for ( j = i+2; j <= i+kwind.degree; j++ )
      xi += knots[j];
    xi /= (double)kwind.degree;
    fgraph[i].x = xi;
    fgraph[i].y = 0.0;
  }
  xn  = (double)xge_KnotWindMapKnot ( &kwind, knots[kwind.degree] );
  xNn = (double)xge_KnotWindMapKnot ( &kwind, knots[kwind.lastknot-kwind.degree] );
  xgeSetForeground ( xgec_Grey );
  for ( i = 0; i < kwind.lastknot-kwind.degree; i++ ) {
    if ( i > 0 )
      fgraph[i-1].y = 0.0;
    fgraph[i].y = 1.0;
    for ( j = 0; j <= dd; j++ ) {
      s = ((double)j/(double)dd);
      t = knots[kwind.degree] +
          s*(knots[kwind.lastknot-kwind.degree]-knots[kwind.degree]);
      mbs_deBoorC2d ( kwind.degree, kwind.lastknot, knots, fgraph, t, &cpoint );
      graph[j].x = (short)(xn + s*(double)(xNn-xn)+0.5);
      graph[j].y = (short)(y + (short)(cpoint.y*(double)yu+0.5));
    }
    xgeDrawLines ( dd+1, graph );
  }

  pkv_FreeScratchMem ( scratchsize );
#undef dd
} /*DisplayBasis*/

static double GrevilleAbscissa ( int i )
{
  double xi;
  int   j;

  xi = knots[i+1];
  for ( j = 2; j <= kwind.degree; j++ )
    xi += knots[i+j];
  return xi/(double)kwind.degree;
} /*GrevilleAbscissa*/

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
/*
printf ( "i = %d,  r = %d,  dr = %d,  u = %6.3f\n", i, r, dr, u );
*/
  for ( j = 0; j <= dr; j++ ) {
    pfknots[j] = knots[i-(dr+1)+j];
    pfknots[dr+1+j] = knots[i+r+j];
  }
/*
WriteArrayd ( "knots", 2*(dr+1), pfknots );
printf ( "\n" );
*/
} /*GetPFKnots*/

void DisplayPolarf ( void )
{
  void     *sp;
  int      i, j, r, dr, nf, kpcs;
  double    u, xi;
  double    *pfknots;
  point3d  *bcpf, *ccpf;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
  if ( !(bcpf = pkv_GetScratchMem ( 3*(kwind.degree+1)*sizeof(point3d) )) )
    goto failure;
  if ( !(pfknots = pkv_GetScratchMemd ( 2*(kwind.degree+1) )) )
    goto failure;

  GetClipLines ( cliplines );
  xgeSetForeground ( xgec_DeepSkyBlue );
  if ( kwind.closed )
    nf = kwind.lastknot+1-kwind.degree;
  else
    nf = kwind.lastknot-kwind.degree;
  if ( funkcja ) {
    ccpf = pkv_GetScratchMem ( (kwind.degree+1)*sizeof(point3d) );
    if ( !ccpf )
      goto failure;
    for ( i = kwind.degree+1; i < nf; i+= r ) {
      GetPFKnots ( i, &r, &u, pfknots );
      if ( r < kwind.degree ) {
        dr = kwind.degree-r;
        for ( j = 0; j <= dr; j++ ) {
          RzutujRPunkt ( &cpoints[i-dr-1+j], &ccpf[j] );
          xi = GrevilleAbscissa ( i-dr-1+j );
          ccpf[j].x = (double)xge_KnotWindMapKnot ( &kwind,
                         u + (xi-u)*((double)kwind.degree/(double)dr) );
        }
        mbs_BSToBezC3d ( dr, 2*dr+1, pfknots, ccpf, &kpcs, NULL, NULL, bcpf );
        mbs_ClipBC2Rd ( 4, cliplines, dr, bcpf, xge_DrawBC2Rd );
      }
    }
  }
  else {
    for ( i = kwind.degree+1; i < nf; i += r ) {
      GetPFKnots ( i, &r, &u, pfknots );
      if ( r < kwind.degree ) {
        dr = kwind.degree-r;
        mbs_BSToBezC3d ( dr, 2*dr+1, pfknots, &cpoints[i-dr-1],
                         &kpcs, NULL, NULL, bcpf );
        for ( j = 0; j <= dr; j++ )
          RzutujRPunkt ( &bcpf[j], &bcpf[j] );
        mbs_ClipBC2Rd ( 4, cliplines, dr, bcpf, xge_DrawBC2Rd );
      }
    }
  }
failure:
  pkv_SetScratchMemTop ( sp );
} /*DisplayPolarf*/

static void DisplayBezHedgehog ( int n, point3d *cp )
{
  double    scale;
  int      i;
  double    t, curvature;
  point2d  cpoint, hpoint;
  vector2d fframe[2];

  scale = xge_LogSlidebarValued ( 0.002, 2.0, curvgraphscale );
  for ( i = 0; i <= curv_graph_dens; i++ ) {
    t = (double)i/(double)curv_graph_dens;
    mbs_BCFrenetC2Rd ( n, cp, t, &cpoint, fframe, &curvature );
    AddVector2Md ( &cpoint, &fframe[1], scale*curvature, &hpoint );
    CameraProjectPoint2d ( &cwind.CPos, &cpoint, &cpoint );
    CameraProjectPoint2d ( &cwind.CPos, &hpoint, &hpoint );
    xgeDrawLine ( (int)cpoint.x, (int)cpoint.y, (int)hpoint.x, (int)hpoint.y );
  }
} /*DisplayBezHedgehog*/

void DisplayCurvGraph ( void )
{
  int i, kpcs;
  point3d *ncp;
  int scratchsize;

  if ( kwind.degree < 2 )
    return;
  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point3d*)pkv_GetScratchMem (
                    scratchsize = (kwind.degree+1)*kpcs*sizeof(point3d) );
  mbs_BSToBezC3d ( kwind.degree, kwind.lastknot, knots, cpoints, &kpcs,
                   NULL, NULL, ncp );

  xgeSetForeground ( xgec_DodgerBlue );

  for ( i = 0; i < kpcs; i++ )
    DisplayBezHedgehog ( kwind.degree, &ncp[i*(kwind.degree+1)] );

  pkv_FreeScratchMem ( scratchsize );
} /*DisplayCurvGraph*/

/* ///////////////////////////////////////////////////////////////////////// */
void FindBoundingBox ( Box2d *box )
{
  int     i;
  point2d p;

  Point3to2d ( &cpoints[0], &p );
  box->x0 = box->x1 = p.x;
  box->y0 = box->y1 = p.y;
  for ( i = 1; i < kwind.lastknot-kwind.degree; i++ ) {    
    Point3to2d ( &cpoints[i], &p );   
    if ( p.x < box->x0 )      box->x0 = p.x;
    else if ( p.x > box->x1 ) box->x1 = p.x;
    if ( p.y < box->y0 )      box->y0 = p.y;
    else if ( p.y > box->y1 ) box->y1 = p.y;
  }
  /* extend it by 5% */
  p.x = (double)(1.025*box->x1-0.025*box->x0);
  p.y = (double)(1.025*box->y1-0.025*box->y0);
  box->x0 = (double)(1.025*box->x0-0.025*box->x1);
  box->y0 = (double)(1.025*box->y0-0.025*box->y1);
  box->x1 = p.x;
  box->y1 = p.y;
 } /*FindBoundingBox*/

void ResetObject ( void )
{
  pkv_SetErrorHandler ( ErrorHandler );

  curve   = true;
  lamana  = true;
  funkcja = false;
  ticks   = false;
  convh   = false;
  nurbs   = false;
  kwind.closed  = false;
  baza    = false;
  polarf  = false;
  curvgr  = false;
  uniform = false;

  kwind.lastknot = 3;
  knots[0] = 0.0;
  knots[1] = 0.0;
  knots[2] = 1.0;
  knots[3] = 1.0;
  kwind.degree  = 1;
  SetPoint3d ( &cpoints[0], -1.0, 0.0, 1.0 );
  SetPoint3d ( &cpoints[1], +1.0, 0.0, 1.0 );
  npoints = kwind.lastknot-kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );

  xge_KnotWindInitMapping ( &kwind, 0.0, 1.0 );
  ProjectCurve ();

  cwind.current_point = 0;
  kwind.current_knot = 1;
} /*ResetObject*/

void ClearPointMarking ( void )
{
  memset ( cpmark, 0, npoints*sizeof(boolean) );
} /*ClearPointMarking*/

void ResizeObject ( void )
{
  xge_KnotWindInitMapping ( &kwind, kwind.umin, kwind.umax );
  UstawFunkcje ();
  ProjectCurve ();
} /*ResizeObject*/
                    
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
  if ( funkcja ) {
    for ( i = 0; i < kwind.degree; i++ )
      cpoints[kwind.clcK+i].y = cpoints[i].y =
        (double)(0.5*(cpoints[kwind.clcK+i].y+cpoints[i].y));
  }
  else {
    for ( i = 0; i < kwind.degree; i++ ) {
      MidPoint3d ( &cpoints[kwind.clcK+i], &cpoints[i], &cpoints[i] );
      cpoints[kwind.clcK+i] = cpoints[i];
    }
  }
  memset ( cpmark, 0, npoints*sizeof(boolean) );
} /*SetClosedCurve*/

boolean DegreeElevation ( void )
{
  void    *sp;
  int     k;
  double   *nkn;
  point3d *ncp;
/*
int Nt, r, s, t, pt;
*/

  if ( kwind.degree < MAX_DEGREE ) {
                             /* compute the final number of knots */
    k = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
    if ( k*kwind.degree+2 < MAX_KNOTS ) {
      sp = pkv_GetScratchMemTop ();
      nkn = (double*)pkv_GetScratchMemd ( MAX_KNOTS );
      ncp = (point3d*)pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
      if ( !nkn || !ncp ) {
        printf ( "Not enough memory\n" );
        exit ( 1 );
      }
      if ( kwind.closed )
        mbs_BSDegElevClosedC3d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                                 &kwind.degree, &kwind.lastknot, nkn, ncp );
      else
        mbs_BSDegElevC3d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                           &kwind.degree, &kwind.lastknot, nkn, ncp, true );
      npoints = kwind.lastknot-kwind.degree;
      memcpy ( knots, nkn, (kwind.lastknot+1)*sizeof(double) );
      memcpy ( cpoints, ncp, npoints*sizeof(point3d) );
      memset ( cpmark, 0, npoints*sizeof(boolean) );
      if ( kwind.closed ) {
        kwind.clcK = kwind.lastknot-2*kwind.degree;
        kwind.clcT = knots[kwind.degree+kwind.clcK]-knots[kwind.degree];
      }
      pkv_SetScratchMemTop ( sp );
      uniform = false;
      if ( funkcja )
        UstawFunkcje ();
      ProjectCurve ();
      return true;
    }
    else
      return false;
  }
  else
    return false;
} /*DegreeElevation*/

boolean DegreeReduction ( void )
{
  void    *sp;
  double   *nkn;
  point3d *ncp;

  sp = pkv_GetScratchMemTop ();
  if ( kwind.degree > 1 ) {
    nkn = (double*)pkv_GetScratchMemd ( MAX_KNOTS );
    ncp = (point3d*)pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
    if ( !nkn || !ncp ) {
      printf ( "Not enough memory\n" );
      exit ( 1 );
    }
    if ( kwind.closed ) {
      if ( !mbs_BSDegRedClosedC3d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                                    &kwind.degree, &kwind.lastknot, nkn, ncp ) )
        goto failure;
      kwind.clcK = kwind.lastknot-2*kwind.degree;
    }
    else {
      if ( !mbs_BSDegRedC3d ( kwind.degree, kwind.lastknot, knots, cpoints, 1,
                              &kwind.degree, &kwind.lastknot, nkn, ncp ) )
        goto failure;
    }
    npoints = kwind.lastknot-kwind.degree;
    memcpy ( knots, nkn, (kwind.lastknot+1)*sizeof(double) );
    memcpy ( cpoints, ncp, npoints*sizeof(point3d) );
    memset ( cpmark, 0, npoints*sizeof(boolean) );
    pkv_SetScratchMemTop ( sp );
    uniform = false;
    if ( funkcja )
      UstawFunkcje ();
    ProjectCurve ();
    return true;
  }

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*DegreeReduction*/

void InsertKnot ( void )
{
  if ( kwind.closed ) {
    mbs_KnotInsClosedC3d ( kwind.degree, &kwind.lastknot, knots,
                           cpoints, kwind.newknot );
    kwind.clcK = kwind.lastknot - 2*kwind.degree;
  }
  else {
    mbs_KnotInsC3d ( kwind.degree, &kwind.lastknot, knots,
                     cpoints, kwind.newknot );
  }
  npoints = kwind.lastknot - kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  uniform = false;
  ProjectCurve ();
} /*InsertKnot*/

void RemoveKnot ( void )
{
  if ( kwind.closed ) {
    mbs_KnotRemoveClosedC3d ( kwind.degree, &kwind.lastknot, knots,
                              cpoints, kwind.current_knot );
    kwind.clcK = kwind.lastknot - 2*kwind.degree;
  }
  else {
    mbs_KnotRemoveC3d ( kwind.degree, &kwind.lastknot, knots,
                        cpoints, kwind.current_knot );
  }
  npoints = kwind.lastknot-kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  uniform = false;
  ProjectCurve ();
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

boolean RefineUniform ( void )
{
  void   *sp;
  double *acp;
  int    lkn;

  if ( !uniform || 2*(kwind.lastknot-kwind.degree) >= MAX_KNOTS )
    return false;
  sp = pkv_GetScratchMemTop ();
  acp = pkv_GetScratchMem ( (2*kwind.lastknot-3*kwind.degree)*sizeof(point3d) );
  if ( !acp )
    goto failure;
  if ( !mbs_LaneRiesenfeldC3d ( kwind.degree, kwind.lastknot, cpoints,
                                &lkn, acp ) )
    goto failure;
  memcpy ( cpoints, acp, (lkn-kwind.degree)*sizeof(point3d) );
  knots[lkn] = knots[lkn-1] = knots[kwind.lastknot-1];
  kwind.lastknot = lkn;
  SetUniformKnots ();
  npoints = kwind.lastknot-kwind.degree;
  memset ( cpmark, 0, npoints*sizeof(boolean) );
  ProjectCurve ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*RefineUniform*/

void UstawFunkcje ( void )
{
  int     i;
  point2d p;

  if ( funkcja ) {
    curvgr = 0;
    if ( nurbs ) {
      nurbs = 0;
      UstawNURBS ();
    }
    ProjectCurve ();
    for ( i = 0; i < kwind.lastknot-kwind.degree; i++ ) {
      rpoints[i].x = (double)xge_KnotWindMapKnot ( &kwind, GrevilleAbscissa ( i ) );
      CameraUnProjectPoint2d ( &cwind.CPos, &rpoints[i], &p );
      Point2to3d ( &p, cpoints[i].z, &cpoints[i] );
    }
  }
  else {
    if ( kwind.closed )
      SetClosedCurve ();
  }
} /*UstawFunkcje*/

void UstawNURBS ( void )
{
  int i;

  if ( !nurbs ) {
    for ( i = 0; i < npoints; i++ ) {
      SetPoint3d ( &cpoints[i], cpoints[i].x/cpoints[i].z,
                   cpoints[i].y/cpoints[i].z, 1.0 );
    }
    ProjectCurve ();
  }
  else
    funkcja = 0;
} /*UstawNURBS*/

void UstawPolarf ( void )
{
  UstawNURBS ();
} /*UstawPolarf*/

void UstawCurvGraph ( void )
{
  if ( curvgr && funkcja ) {
    funkcja = 0;
    UstawFunkcje ();
  }
} /*UstawCurvGraph*/

boolean UstawZamknieta ( void )
{
  if ( kwind.closed ) {
    if ( kwind.lastknot <= 3*kwind.degree ) {
      kwind.closed = false;
      return false;
    }
    SetClosedCurve ();
  }
  xge_KnotWindInitMapping ( &kwind, (double)min(0.0, kwind.knots[1]),
                            (double)max(1.0, knots[kwind.lastknot-1]) );
  ProjectCurve ();
  return true;
} /*UstawZamknieta*/

void ExportPovRay ( void )
{
  void    *sp;
  FILE    *f; 
  int     kpcs, i, j, k;
  point3d *ncp;

  sp = pkv_GetScratchMemTop ();
  f = fopen ( "pokrzyw.pov", "w+" );
  kpcs = mbs_NumKnotIntervalsd ( kwind.degree, kwind.lastknot, knots );
  ncp = (point3d*)pkv_GetScratchMem ( (kwind.degree+1)*kpcs*sizeof(point3d) );
  if ( !ncp )
    exit ( 1 );
  mbs_BSToBezC3d ( kwind.degree, kwind.lastknot, knots, cpoints,
                   &kpcs, NULL, NULL, ncp );

  fprintf ( f, "degree = %d, arcs = %d\n", kwind.degree, kpcs );
  for ( i = k = 0; i < kpcs; i++ ) {
    for ( j = 0; j <= kwind.degree; j++, k++ ) {
      fprintf ( f, "<%f,%f>", ncp[k].x/ncp[k].z, ncp[k].y/ncp[k].z );
      if ( j < kwind.degree )
        fprintf ( f, "," );
    }
    fprintf ( f, "\n" );
  }
  fclose ( f );
  printf ( "%s\n", "pokrzyw.pov" );
  pkv_SetScratchMemTop ( sp );
} /*ExportPovRay*/

