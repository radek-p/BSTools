
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
#include "ed3dswidgets.h"


static void ProjectEqMerPoint ( point3d *p, point2d *q )
{
  point2d a;
  xge_2Dwind *cwind;

  cwind = eq_cwind.er->data0;
  Point3to2d ( p, &a );
  CameraProjectPoint2d ( &cwind->CPos, &a, q );
} /*ProjectEqMerPoint*/

static void ProjectEqMerRPoint ( point3d *p, point3d *q )
{
  point2d ap;

  ProjectEqMerPoint ( p, &ap );
  Point2to3d ( &ap, p->z, q );
} /*ProjectEqMerRPoint*/

void ProjectEqMerCurve ( void )
{
  int i;

  for ( i = 0; i < neqmerpoints; i++ )
    ProjectEqMerPoint ( &eqmer_cpoints[i], &eqmer_rpoints[i] );
  for ( i = 0; i < neqmerpoints-1; i++ ) {
    MidPoint3d ( &eqmer_cpoints[i], &eqmer_cpoints[i+1],
                 &eqmer_cpoints[i+neqmerpoints] );
    ProjectEqMerPoint ( &eqmer_cpoints[i+neqmerpoints],
                        &eqmer_rpoints[i+neqmerpoints] );
  }
} /*ProjectEqMerCurve*/

boolean SetClosedEqMer ( void )
{
  int          i, clcK, degree, lastknot;
  double        clcT, *knots;
  xge_KnotWind *ckwind;

  ckwind = eq_ckwind.er->data0;
  knots = ckwind->knots;
  degree = ckwind->degree;
  lastknot = ckwind->lastknot;
  clcK = 0;  clcT = 0.0;
  if ( (ckwind->closed = eqmer_closed) ) {
    if ( lastknot <= 3*degree ) {
      xge_DisplayErrorMessage ( ErrorMsgNotEnoughKnots, -1 );
      ckwind->closed = false;
      return false;
    }
    clcK = ckwind->clcK = lastknot-2*degree;
    clcT = ckwind->clcT = knots[clcK+degree]-knots[degree];
    for ( i = 1; i < 2*degree; i++ )
      knots[i] = (double)(0.5*(knots[i]+knots[i+clcK]-clcT));
    for ( i = 1; i < 2*degree; i++ )
      knots[clcK+i] = knots[i]+clcT;
    knots[0] = knots[1];
    knots[lastknot] = knots[lastknot-1];
    for ( i = 0; i < degree; i++ ) {
      MidPoint3d ( &eqmer_cpoints[i+clcK], &eqmer_cpoints[i],
                   &eqmer_cpoints[i] );
      eqmer_cpoints[i+clcK] = eqmer_cpoints[i];
    }
    ClearPointMarking ( neqmerpoints, eqmer_mkpoints );
  }
  ProjectEqMerCurve ();
  if ( equator ) {
    if ( (kwind.closed_u = ckwind->closed) ) {
      kwind.clcKu = clcK;
      kwind.clcTu = clcT;
    }
  }
  else {
    if ( (kwind.closed_v = ckwind->closed) ) {
      kwind.clcKv = clcK;
      kwind.clcTu = clcT;
    }
  }
  return true;
} /*SetClosedEqMer*/

boolean FindNearestEqMerCPoint ( int x, int y, int mindist )
{
  int        i, m;
  double      d, e;
  xge_2Dwind *cwind;

  cwind = eq_cwind.er->data0;
  if ( eqmer_nurbs )
    m = 2*neqmerpoints-1;
  else m = neqmerpoints;
  e = (double)mindist;
  for ( i = 0; i < m; i++ ) {
    d = (double)(fabs((double)x-eqmer_rpoints[i].x) +
                fabs((double)y-eqmer_rpoints[i].y));
    if ( d < e ) {
      cwind->current_point = i;
      e = d;
    }
  }
  return (boolean)(e < mindist);
} /*FindNearestEqMerCPoint*/

void DisplayEqMerControlPolygon ( void )
{
  void *sp;
  int i;
  XPoint *poly;

  sp = pkv_GetScratchMemTop ();
  poly = pkv_GetScratchMem ( MAX_KNOTS*sizeof(XPoint) );
  if ( poly ) {
    for ( i = 0; i < neqmerpoints; i++ ) {
      poly[i].x = (short)(eqmer_rpoints[i].x+0.5);
      poly[i].y = (short)(eqmer_rpoints[i].y+0.5);
    }
    xgeSetForeground ( xgec_Green2 );
    xgeDrawLines ( neqmerpoints, poly );
  }
  pkv_SetScratchMemTop ( sp );
} /*DisplayEqMerControlPolygon*/

static void GetClipLines ( vector3d cliplines[4] )
{
  xge_2Dwind *cwind;
  CameraRecd *CPos;

  cwind = eq_cwind.er->data0;
  CPos = &cwind->CPos;
  SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-CPos->xmin) );
  SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-CPos->ymin) );
  SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(CPos->xmin+CPos->width) ); 
  SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(CPos->ymin+CPos->height) );
} /*GetClipLines*/

static void DisplayTick ( point3d *p1, point3d *p2 )
{
  vector2d v;
  double x, y;
  int x1, y1, x2, y2;

  x = p1->x/p1->z;
  y = p1->y/p1->z;
  SetVector2d ( &v, p2->x/p2->z-x, p2->y/p2->z-y );
  NormalizeVector2d ( &v );
  MultVector2d ( 4.5, &v, &v );
  x1 = (int)(x+v.y);
  y1 = (int)(y-v.x);
  x2 = (int)(x-v.y);
  y2 = (int)(y+v.x);
  xgeDrawLine ( x1, y1, x2, y2 );
} /*DisplayTick*/

void DisplayEqMerCurve ( void )
{
  void         *sp;
  int          i, j, kpcs;
  point3d      *ncp;
  vector3d     cliplines[4];
  xge_KnotWind *ckwind;

  sp = pkv_GetScratchMemTop ();
  ckwind = eq_ckwind.er->data0;
    /* Get the clipping line equation coefficients */
  GetClipLines ( cliplines );
    /* divide the spline curve into rational Bezier arcs */
  kpcs = mbs_NumKnotIntervalsd ( ckwind->degree,
                                 ckwind->lastknot, ckwind->knots );
  ncp = (point3d*)pkv_GetScratchMem ( (ckwind->degree+1)*kpcs*sizeof(point3d) );
  if ( !ncp )
    goto failure;
  mbs_BSToBezC3d ( ckwind->degree, ckwind->lastknot,
                   ckwind->knots, eqmer_cpoints, &kpcs, NULL, NULL, ncp );
  for ( i = 0; i < kpcs*(ckwind->degree+1); i++ )
    ProjectEqMerRPoint ( &ncp[i], &ncp[i] );

    /* clip and draw the arcs */
  xgeSetForeground ( xgec_White );
  for ( i = 0; i < kpcs; i++ ) {
    mbs_ClipBC2Rd ( 4, cliplines, ckwind->degree, &ncp[i*(ckwind->degree+1)],
                    xge_DrawBC2Rd );
  }

  if ( display_eqmer_ticks ) {
    xgeSetLineAttributes ( 2, LineSolid, CapButt, JoinMiter );
    xgeSetForeground ( xgec_DeepPink );
    for ( i = 0; i < kpcs; i++ ) {
      j = i*(ckwind->degree+1);
      DisplayTick ( &ncp[j], &ncp[j+1] );
    }
    j = (ckwind->degree+1)*kpcs-1;
    DisplayTick ( &ncp[j], &ncp[j-1] );
    xgeSetLineAttributes ( 1, LineSolid, CapButt, JoinMiter );
  }

failure:
  pkv_SetScratchMemTop ( sp );
} /*DisplayEqMerCurve*/

void DisplayEqMerControlPoints ( void )
{
  int   i;
  short x, y;

  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < neqmerpoints; i++ )
    if ( !eqmer_mkpoints[i] ) {
      x = (short)(eqmer_rpoints[i].x+0.5);
      y = (short)(eqmer_rpoints[i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < neqmerpoints; i++ )
    if ( eqmer_mkpoints[i] ) {
      x = (short)(eqmer_rpoints[i].x+0.5);
      y = (short)(eqmer_rpoints[i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
  if ( eqmer_nurbs ) {
    xgeSetForeground ( xgec_HotPink );
    for ( i = neqmerpoints; i < 2*neqmerpoints; i++ ) {
      x = (short)(eqmer_rpoints[i].x+0.5);
      y = (short)(eqmer_rpoints[i].y+0.5);
      xgeFillRectangle ( 3, 3, x-1, y-1 );
    }
  }
} /*DisplayEqMerControlPoints*/

void DisplayEqMerAxes ( void )
{
  point2d    p, q;
  short      x, y;
  xge_2Dwind *cwind;

  cwind = eq_cwind.er->data0;
  SetPoint2d ( &p, 0.0, 0.0 );
  CameraProjectPoint2d ( &cwind->CPos, &p, &q );
  x = (short)(q.x+0.5);
  y = (short)(q.y+0.5);
  xgeSetForeground ( xgec_Grey3 );
  if ( x >= cwind->er->x && x < cwind->er->x+cwind->er->w )
    xgeDrawLine ( x, cwind->er->y, x, cwind->er->y+cwind->er->h );
  if ( y >= cwind->er->y && y < cwind->er->y+cwind->er->h )
    xgeDrawLine ( cwind->er->x, y, cwind->er->x+cwind->er->w, y );
} /*DisplayEqMerAxes*/

void DisplayEqMerBezierPolygons ( void )
{
} /*DisplayEqMerBezierPolygons*/

void InsertEqMerKnot ( void )
{
  xge_KnotWind *ckwind;

  ckwind = eq_ckwind.er->data0;
  if ( ckwind->closed ) {
    mbs_KnotInsClosedC3d ( ckwind->degree, &ckwind->lastknot, ckwind->knots,
                           eqmer_cpoints, ckwind->newknot );
    ckwind->clcK = ckwind->lastknot - 2*ckwind->degree;
  }
  else {
    mbs_KnotInsC3d ( ckwind->degree, &ckwind->lastknot, ckwind->knots,
                     eqmer_cpoints, ckwind->newknot );
  }
  neqmerpoints = ckwind->lastknot - ckwind->degree;
  memset ( eqmer_mkpoints, 0, neqmerpoints*sizeof(byte) );
  ProjectEqMerCurve ();
} /*InsertEqMerKnot*/

void RemoveEqMerKnot ( void )
{
  xge_KnotWind *ckwind;

  ckwind = eq_ckwind.er->data0;
  if ( ckwind->closed ) {
    mbs_KnotRemoveClosedC3d ( ckwind->degree, &ckwind->lastknot, ckwind->knots,
                              eqmer_cpoints, ckwind->current_knot );
    ckwind->clcK = ckwind->lastknot - 2*ckwind->degree;
  }
  else {
    mbs_KnotRemoveC3d ( ckwind->degree, &ckwind->lastknot, ckwind->knots,
                        eqmer_cpoints, ckwind->current_knot );
  }
  neqmerpoints = ckwind->lastknot - ckwind->degree;
  memset ( eqmer_mkpoints, 0, neqmerpoints*sizeof(byte) );
  ProjectEqMerCurve ();
} /*RemoveEqMerKnot*/

boolean EqMerDegreeElevation ( void )
{
  void         *sp;
  int          k;
  double        *nkn;
  point3d      *ncp;
  xge_KnotWind *ckwind;

  ckwind = eq_ckwind.er->data0;
  if ( ckwind->degree < MAX_DEGREE ) {
                             /* compute the final number of knots */
    k = mbs_NumKnotIntervalsd ( ckwind->degree, ckwind->lastknot, ckwind->knots );
    if ( k*ckwind->degree+2 < MAX_KNOTS ) {
      sp = pkv_GetScratchMemTop ();
      nkn = (double*)pkv_GetScratchMemd ( MAX_KNOTS );
      ncp = (point3d*)pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
      if ( !nkn || !ncp ) {
        printf ( "Not enough memory\n" );
        exit ( 1 );
      }
      if ( ckwind->closed )
        mbs_BSDegElevClosedC3d ( ckwind->degree, ckwind->lastknot,
                                 ckwind->knots, eqmer_cpoints, 1,
                                 &ckwind->degree, &ckwind->lastknot, nkn, ncp );
      else
        mbs_BSDegElevC3d ( ckwind->degree, ckwind->lastknot,
                           ckwind->knots, eqmer_cpoints, 1,
                           &ckwind->degree, &ckwind->lastknot, nkn, ncp, true );
      neqmerpoints = ckwind->lastknot-ckwind->degree;
      memcpy ( ckwind->knots, nkn, (ckwind->lastknot+1)*sizeof(double) );
      memcpy ( eqmer_cpoints, ncp, neqmerpoints*sizeof(point3d) );
      memset ( eqmer_mkpoints, 0, neqmerpoints*sizeof(boolean) );
      if ( ckwind->closed ) {
        ckwind->clcK = ckwind->lastknot-2*ckwind->degree;
        ckwind->clcT = ckwind->knots[ckwind->degree+ckwind->clcK]-
                       ckwind->knots[ckwind->degree];
      }
      pkv_SetScratchMemTop ( sp );
      ProjectEqMerCurve ();
      return true;
    }
    else
      return false;
  }
  else
    return false;
} /*EqMerDegreeElevation*/

boolean EqMerDegreeReduction ( void )
{
  void        *sp;
  double       *nkn;
  point3d     *ncp;
  xge_KnotWind *ckwind;

  ckwind = eq_ckwind.er->data0;
  sp = pkv_GetScratchMemTop ();
  if ( ckwind->degree > 1 ) {
    nkn = (double*)pkv_GetScratchMemd ( MAX_KNOTS );
    ncp = (point3d*)pkv_GetScratchMem ( MAX_KNOTS*sizeof(point3d) );
    if ( !nkn || !ncp ) {
      printf ( "Not enough memory\n" );
      exit ( 1 );
    }
    if ( ckwind->closed ) {
      if ( !mbs_BSDegRedClosedC3d ( ckwind->degree, ckwind->lastknot,
                                    ckwind->knots, eqmer_cpoints, 1,
                                    &ckwind->degree, &ckwind->lastknot, nkn, ncp ) )
        goto failure;
      ckwind->clcK = ckwind->lastknot-2*ckwind->degree;
    }
    else {
      if ( !mbs_BSDegRedC3d ( ckwind->degree, ckwind->lastknot,
                              ckwind->knots, eqmer_cpoints, 1,
                              &ckwind->degree, &ckwind->lastknot, nkn, ncp ) )
        goto failure;
    }
    neqmerpoints = ckwind->lastknot-ckwind->degree;
    memcpy ( ckwind->knots, nkn, (ckwind->lastknot+1)*sizeof(double) );
    memcpy ( eqmer_cpoints, ncp, neqmerpoints*sizeof(point3d) );
    memset ( eqmer_mkpoints, 0, neqmerpoints*sizeof(boolean) );
    pkv_SetScratchMemTop ( sp );
    ProjectEqMerCurve ();
    return true;
  }

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*EqMerDegreeReduction*/

void EqMerSetCPoint ( short x, short y )
{
  int          i;
  double        w, r, t;
  double        x1, y1, x2, y2;
  vector2d     v1, v2;
  point3d      q, s;
  xge_2Dwind   *cwind;
  xge_KnotWind *ckwind;

  cwind = eq_cwind.er->data0;
  ckwind = eq_ckwind.er->data0;
  i = cwind->current_point;
  if ( ckwind->closed ) {
    if ( (i < neqmerpoints && i >= ckwind->clcK) ||
         (eqmer_nurbs && i >= neqmerpoints+ckwind->clcK) )
      i -= ckwind->clcK;
  }

  if ( eqmer_nurbs )
    w = eqmer_cpoints[i].z;
    else w = 1.0;
  SetPoint3d ( &q, (double)x, (double)y, 1.0 );
  CameraUnProjectPoint3d ( &cwind->CPos, &q, &s );
  MultVector3d ( w, &s, &s );

  if ( i < neqmerpoints ) {
    eqmer_cpoints[i] = s;
    if ( ckwind->closed ) {
      if ( i < ckwind->degree )
        eqmer_cpoints[i+ckwind->clcK] = eqmer_cpoints[i];
      else if ( i >= ckwind->clcK )
        eqmer_cpoints[i-ckwind->clcK] = eqmer_cpoints[i];
    }
  }  
  else { /* moving a Farin point to change weight */
    x1 = eqmer_rpoints[i-neqmerpoints].x;
    y1 = eqmer_rpoints[i-neqmerpoints].y;
    x2 = eqmer_rpoints[i-neqmerpoints+1].x;
    y2 = eqmer_rpoints[i-neqmerpoints+1].y;
    SetVector2d ( &v1, x2-x1, y2-y1 );
    r = (double)DotProduct2d ( &v1, &v1 );
    if ( r >= 1 ) {
      SetVector2d ( &v2, (double)x-x1, (double)y-y1 );
      t = (double)DotProduct2d ( &v1, &v2 )/r;
      t = (double)max ( t, 0.01 );  t = (double)min ( t, 0.99 );
      i -= neqmerpoints;
      w = (double)(eqmer_cpoints[i].z*(t/(1.0-t)));   /* the new weight */
      w = (double)min ( w, 10000.0 );  w = (double)max ( w, 0.0001 );
      i++;
      Point2to3d ( &eqmer_rpoints[i], 1.0, &q );
      CameraUnProjectPoint3d ( &cwind->CPos, &q, &s );
      MultVector3d ( w, &s, &eqmer_cpoints[i] );
      if ( ckwind->closed ) {
        if ( i < ckwind->degree )
          eqmer_cpoints[i+ckwind->clcK] = eqmer_cpoints[i];
        else if ( i >= ckwind->clcK )
          eqmer_cpoints[i-ckwind->clcK] = eqmer_cpoints[i];
      }
    }
  }
  ProjectEqMerCurve ();
} /*EqMerSetCPoint*/

void SelectEqMerPoints ( const Box2s *sbox, boolean mk )
{
  MarkPoints ( neqmerpoints, eqmer_rpoints, eqmer_mkpoints, sbox, mk );
} /*SelectEqMerPoints*/

void ClearEqMerPointMarking ( void )
{
  ClearPointMarking ( neqmerpoints, eqmer_mkpoints );
} /*ClearEqMerPointMarking*/

void SaveEqMerControlPoints ( void )
{
  memcpy ( savedcpoints, eqmer_cpoints, neqmerpoints*sizeof(point3d) );
} /*SaveEqMerControlPoints*/

void TransformEqMerMarkedControlPoints ( trans2d *tr )
{
  int i;
  point3d *savedcp;

  savedcp = (void*)savedcpoints;
  for ( i = 0; i < neqmerpoints; i++ )
    if ( eqmer_mkpoints[i] )
      Trans2Point3d ( tr, &savedcp[i], &eqmer_cpoints[i] );
  ProjectEqMerCurve ();
} /*TransformEqMerMarkedControlPoints*/

void ResetEquator ( void )
{
  eq_ckwind.degree = 1;
  eq_ckwind.lastknot = 3;
  eq_ckwind.closed = false;
  equator_nurbs = false;
  if ( equator )
    eqmer_nurbs = false;
  equator_knots[0] = equator_knots[1] = 0.0;
  equator_knots[2] = equator_knots[3] = 1.0;
  SetPoint3d ( &equator_cpoints[0], -1.0, 0.0, 1.0 );
  SetPoint3d ( &equator_cpoints[1],  1.0, 0.0, 1.0 );
  neqmerpoints = 2;
} /*ResetEquator*/

void ResetMeridian ( void )
{
  mer_ckwind.degree = 1;
  mer_ckwind.lastknot = 3;
  mer_ckwind.closed = false;
  meridian_nurbs = false;
  if ( meridian )
    eqmer_nurbs = false;
  meridian_knots[0] = meridian_knots[1] = 0.0;
  meridian_knots[2] = meridian_knots[3] = 1.0;
  SetPoint3d ( &meridian_cpoints[0], 1.0, -1.0, 1.0 );
  SetPoint3d ( &meridian_cpoints[1], 1.0,  1.0, 1.0 );
  neqmerpoints = 2;
} /*ResetMeridian*/

void ResetEqMer ( void )
{
  if ( equator )
    ResetEquator ();
  else
    ResetMeridian ();
  ClearEqMerPointMarking ();
  ProjectEqMerCurve ();
  if ( bind_spr )
    BindSphericalProduct ();
} /*ResetEqMer*/

void FindEqMerRefBox ( Box2d *box )
{
  int     i;
  point2d p;

  /* the orogin of the coordinate system must also be within the box */
  box->x0 = box->x1 = 0.0;
  box->y0 = box->y1 = 0.0;
  for ( i = 0; i < neqmerpoints; i++ ) {
    Point3to2d ( &eqmer_cpoints[i], &p );
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
} /*FindEqMerRefBox*/

void SelectEquator ( void )
{
  eqdeg.er->data1 = &eq_ckwind.degree;
  equator = true;
  meridian = false;
  eq_cwind.er->data0 = &eq_cwind;
  eq_ckwind.er->data0 = &eq_ckwind;
  eqmer_cpoints   = &equator_cpoints[0];
  eqmer_rpoints   = &equator_rpoints[0];
  eqmer_mkpoints  = &equator_mkpoints[0];
  eqmer_closed    = eq_ckwind.closed;
  eqmer_nurbs     = equator_nurbs;
  neqmerpoints    = eq_ckwind.lastknot-eq_ckwind.degree;
  ProjectEqMerCurve ();
} /*SelectEquator*/

void SelectMeridian ( void )
{
  eqdeg.er->data1 = &mer_ckwind.degree;
  equator = false;
  meridian = true;
  eq_cwind.er->data0 = &mer_cwind;
  eq_ckwind.er->data0 = &mer_ckwind;
  eqmer_cpoints   = &meridian_cpoints[0];
  eqmer_rpoints   = &meridian_rpoints[0];
  eqmer_mkpoints  = &meridian_mkpoints[0];
  eqmer_closed    = mer_ckwind.closed;
  eqmer_nurbs     = meridian_nurbs;
  neqmerpoints    = mer_ckwind.lastknot-mer_ckwind.degree;
  ProjectEqMerCurve ();
} /*SelectMeridian*/

void SetupNURBSEqMer ( void )
{
  int i;

  if ( equator )
    equator_nurbs = eqmer_nurbs;
  else
    meridian_nurbs = eqmer_nurbs;
  if ( !eqmer_nurbs ) {
    for ( i = 0; i < neqmerpoints; i++ )
      SetPoint3d ( &eqmer_cpoints[i], eqmer_cpoints[i].x/eqmer_cpoints[i].z,
                   eqmer_cpoints[i].y/eqmer_cpoints[i].z, 1.0 );
  }
  ProjectEqMerCurve ();
} /*SetupNURBSEqMer*/

/* ////////////////////////////////////////////////////////////////////////// */
void BindSphericalProduct ( void )
{
  int i;

  degree_u   = eq_ckwind.degree;
  lastknot_u = eq_ckwind.lastknot;
  if ( (kwind.closed_u = eq_ckwind.closed) ) {
    kwind.clcKu = eq_ckwind.clcK;
    kwind.clcTu = eq_ckwind.clcT;
  }
  memcpy ( knots_u, equator_knots, (eq_ckwind.lastknot+1)*sizeof(double) );
  degree_v   = mer_ckwind.degree;
  lastknot_v = mer_ckwind.lastknot;
  if ( (kwind.closed_v = mer_ckwind.closed) ) {
    kwind.clcKv = mer_ckwind.clcK;
    kwind.clcTv = mer_ckwind.clcT;
  }
  memcpy ( knots_v, meridian_knots, (mer_ckwind.lastknot+1)*sizeof(double) );
  mbs_SphericalProductRd (
       eq_ckwind.degree, eq_ckwind.lastknot, equator_cpoints,
       mer_ckwind.degree, mer_ckwind.lastknot, meridian_cpoints,
       4*(mer_ckwind.lastknot-mer_ckwind.degree), cpoints );
  ClearPointMarking ( (eq_ckwind.lastknot-eq_ckwind.degree)*
                      (mer_ckwind.lastknot-mer_ckwind.degree), mkpoints );
  for ( i = 0; i < 4; i++ )
    ProjectSurface ( i );
} /*BindSphericalProduct*/

