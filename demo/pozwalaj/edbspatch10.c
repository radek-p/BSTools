
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>   
#include <stdio.h>
#include <math.h>  
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h" 
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "g2blendingd.h"
#include "egholed.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "render.h"
#include "editor.h"
#include "editor_bsp.h"


void GeomObjectBSplinePatchCopyTrimmedDomain ( GO_BSplinePatch *obj,
                                               GO_BSplinePatch *copy )
{
  BSP_TrimmedDomain *trpd;
  mbs_polycurved    *bound;
  int               i, nelem, dim, deg, lkn, ncp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH ||
       copy->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( !obj->trpd )
    return;
  nelem = obj->trpd->nelem;
  trpd = copy->trpd = malloc ( sizeof(BSP_TrimmedDomain) );
  if ( !trpd )
    return;
  memset ( trpd, 0, sizeof(BSP_TrimmedDomain) );
  bound = trpd->bound = malloc ( nelem*sizeof(mbs_polycurved) );
  if ( !bound ) {
    free ( trpd );
    return;
  }
  memset ( bound, 0, nelem*sizeof(mbs_polycurved) );
  trpd->nelem = nelem;
  for ( i = 0; i < nelem; i++ ) {
    bound[i].ident    = -1;
    bound[i].closing  = obj->trpd->bound[i].closing;
    bound[i].cpdimen  = dim = obj->trpd->bound[i].cpdimen;
    bound[i].degree   = deg = obj->trpd->bound[i].degree;
    bound[i].lastknot = lkn = obj->trpd->bound[i].lastknot;
    if ( obj->trpd->bound[i].knots ) {  /* a B-spline curve */
      bound[i].knots = malloc ( (lkn+1)*sizeof(double) );
      if ( bound[i].knots ) {
        memcpy ( bound[i].knots, obj->trpd->bound[i].knots,
                 (lkn+1)*sizeof(double) );
        ncp = lkn-deg;
      }
      else
        goto failure;
    }
    else if ( deg > 1 )  /* a Bezier curve */
      ncp = deg+1;
    else                 /* a polyline */
      ncp = lkn+1;
    bound[i].points   = malloc ( ncp*dim*sizeof(double) );
    if ( bound[i].points )
      memcpy ( bound[i].points, obj->trpd->bound[i].points,
               ncp*dim*sizeof(double) );
    else
      goto failure;
  }
  return;

failure:
  GeomObjectBSplinePatchDeleteTrimmedDomain ( copy );
  return;
} /*GeomObjectBSplinePatchCopyTrimmedDomain*/

void GeomObjectBSplinePatchDeleteTrimmedDomain ( GO_BSplinePatch *obj )
{
  int               i, nelem;
  BSP_TrimmedDomain *trpd;
  mbs_polycurved    *bound;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  trpd = obj->trpd;
  nelem = trpd->nelem;
  bound = trpd->bound;
    /* deallocate all elements of the domain boundary */
  for ( i = 0; i < nelem; i++ ) {
    if ( bound[i].knots )
      free ( bound[i].knots );
    if ( bound[i].points )
      free ( bound[i].points );
  }
  if ( trpd->bound )
    free ( trpd->bound );
  free ( trpd );
  obj->trpd = NULL;
  obj->trimmed = false;
} /*GeomObjectBSplinePatchDeleteTrimmedDomain*/

void GeomObjectBSplinePatchEnterTrimmedDomain ( GO_BSplinePatch *obj,
                                                int nelem, mbs_polycurved *elem )
{
  BSP_TrimmedDomain *trpd;
  int               i;

  trpd = NULL;
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    goto failure;
  if ( obj->trpd )
    GeomObjectBSplinePatchDeleteTrimmedDomain ( obj );
  trpd = malloc ( sizeof(BSP_TrimmedDomain) );
  if ( !trpd )
    goto failure;
  trpd->nelem = nelem;
  trpd->bound = elem;
  obj->trpd = trpd;
  obj->trimmed = true;
  return;

failure:
        /* deallocate all knot sequences and control points; */
        /* the elem array will be deallocated by the caller. */
  for ( i = 0; i < nelem; i++ ) {
    if ( elem[i].knots )  free ( elem[i].knots );
    if ( elem[i].points ) free ( elem[i].points );
  }
  free ( elem );
  if ( trpd ) free ( trpd );
  obj->trpd = NULL;
  obj->trimmed = false;
} /*GeomObjectBSplinePatchEnterTrimmedDomain*/

/* ////////////////////////////////////////////////////////////////////////// */
static void _GetClipLines ( CameraRecd *CPos, vector3d cliplines[4] )
{
  SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)(-CPos->xmin) );
  SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)(-CPos->ymin) );
  SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(CPos->xmin+CPos->width) );
  SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(CPos->ymin+CPos->height) );
} /*_GetClipLines*/
        
static void _DrawTrdBSCurve ( CameraRecd *CPos, int degree, int lastknot,
                              double *knots, int cpdim, double *cp )
{
  void *sp;
  double   *pcp, *bcp, w;
  vector3d cliplines[4];
  int      i, ncp, kpcs;
  point2d  q;

  sp = pkv_GetScratchMemTop ();
  kpcs = mbs_NumKnotIntervalsd ( degree, lastknot, knots );
  if ( kpcs < 1 )
    goto way_out;
  ncp = lastknot-degree;
  _GetClipLines ( CPos, cliplines );
  switch ( cpdim ) {
case 2:
    pcp = pkv_GetScratchMem ( ncp*sizeof(point2d) );
    bcp = pkv_GetScratchMem ( kpcs*(degree+1)*sizeof(point2d) );
    if ( !pcp || !bcp )
      goto way_out;
    for ( i = 0; i < ncp; i++ )
      CameraProjectPoint2d ( CPos, (point2d*)&cp[2*i], (point2d*)&pcp[2*i] );
    mbs_BSToBezC2d ( degree, lastknot, knots, pcp, &kpcs, NULL, NULL, bcp );
    for ( i = 0; i < kpcs; i++ )
      mbs_ClipBC2d ( 4, cliplines, degree, (point2d*)&bcp[i*2*(degree+1)],
                     xge_DrawBC2d );
    break;
case 3:
    pcp = pkv_GetScratchMem ( ncp*sizeof(point3d) );
    bcp = pkv_GetScratchMem ( kpcs*(degree+1)*sizeof(point3d) );
    if ( !pcp || !bcp )
      goto way_out;
    for ( i = 0; i < ncp; i++ ) {
      w = cp[3*i+2];
      SetPoint2d ( &q, cp[3*i]/w, cp[3*i+1]/w );
      CameraProjectPoint2d ( CPos, &q, &q );
      SetPoint3d ( (point3d*)&pcp[3*i], q.x*w, q.y*w, w );
    }
    mbs_BSToBezC3d ( degree, lastknot, knots, pcp, &kpcs, NULL, NULL, bcp );
    for ( i = 0; i < kpcs; i++ )
      mbs_ClipBC2Rd ( 4, cliplines, degree, (point3d*)&bcp[i*3*(degree+1)],
                      xge_DrawBC2Rd );
    break;
default:
    break;
  }

way_out:
  pkv_SetScratchMemTop ( sp );
} /*_DrawTrdBSCurve*/

static void _DrawTrdBezCurve ( CameraRecd *CPos, int degree,
                               int cpdim, double *cp )
{
  void     *sp;
  double   *pcp, w;
  vector3d cliplines[4];
  int      i;
  point2d  q;

  sp = pkv_GetScratchMemTop ();
  _GetClipLines ( CPos, cliplines );

  switch ( cpdim ) {
case 2:
    pcp = pkv_GetScratchMem ( (degree+1)*sizeof(point2d) );
    if ( !pcp )
      goto way_out;
    for ( i = 0; i <= degree; i++ )
      CameraProjectPoint2d ( CPos, (point2d*)&cp[2*i], (point2d*)&pcp[2*i] );
    mbs_ClipBC2d ( 4, cliplines, degree, (point2d*)pcp, xge_DrawBC2d );
    break;
case 3:
    pcp = pkv_GetScratchMem ( (degree+1)*sizeof(point3d) );
    if ( !pcp )
      goto way_out;
    for ( i = 0; i <= degree; i++ ) {
      w = cp[3*i+2];
      SetPoint2d ( &q, cp[3*i]/w, cp[3*i+1]/w );
      CameraProjectPoint2d ( CPos, &q, &q );
      SetPoint3d ( (point3d*)&pcp[3*i], q.x*w, q.y*w, w );
    }
    mbs_ClipBC2Rd ( 4, cliplines, degree, (point3d*)pcp, xge_DrawBC2Rd );
    break;
default:
    break;
  }

way_out:
  pkv_SetScratchMemTop ( sp );
} /*_DrawTrdBezCurve*/

static void _DrawTrdPolyline ( CameraRecd *CPos, int npoints,
                               int cpdim, double *cp )
{
  void    *sp;
  vector3d cliplines[4];
  point2d p, *q;
  int     i;

  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMem ( npoints*sizeof(point2d) );
  if ( !q )
    goto way_out;
  _GetClipLines ( CPos, cliplines );
  switch ( cpdim ) {
case 2:
    for ( i = 0; i < npoints; i++ )
      CameraProjectPoint2d ( CPos, (point2d*)&cp[2*i], &q[i] );
    goto draw_it;
case 3:
    for ( i = 0; i < npoints; i++ ) {
      Point3to2d ( (point3d*)&cp[3*i], &p );
      CameraProjectPoint2d ( CPos, &p, &q[i] );
    }
draw_it:
    for ( i = 0; i < npoints-1; i++ )
      mbs_ClipBC2d ( 4, cliplines, 1, &q[i], xge_DrawBC2d );
    break;
default:
    break;
  }
way_out:
  pkv_SetScratchMemTop ( sp );
} /*_DrawTrdPolyline*/

void GeomObjectBSplinePatchDrawTrimmedDomain ( GO_BSplinePatch *obj,  
                                               CameraRecd *CPos )
{
  BSP_TrimmedDomain *trpd;
  int               nelem, i;
  mbs_polycurved    *elem;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( !obj->trimmed || !obj->trpd )
    return;
  trpd  = obj->trpd;
  elem  = trpd->bound;
  nelem = trpd->nelem;
  xgeSetForeground ( xgec_White );
  for ( i = 0; i < nelem; i++ ) {
    if ( elem[i].knots )  /* draw a B-spline curve */
      _DrawTrdBSCurve ( CPos, elem[i].degree, elem[i].lastknot, 
                        elem[i].knots, elem[i].cpdimen, elem[i].points );
    else if ( elem[i].degree > 1 ) /* draw a Bezier curve */
      _DrawTrdBezCurve ( CPos, elem[i].degree,
                         elem[i].cpdimen, elem[i].points );
    else  /* draw a polyline */
      _DrawTrdPolyline ( CPos, elem[i].lastknot+1,
                         elem[i].cpdimen, elem[i].points );
  }
} /*GeomObjectBSplinePatchDrawTrimmedDomain*/

static void _FlipXYPoints ( int dim, int npoints, double *points )
{
  int    i, k;
  double s;

  for ( i = k = 0;  i < npoints;  i++, k += dim )
    { s = points[k];  points[k] = points[k+1];  points[k+1] = s; }
} /*_FlipXYPoints*/

void GeomObjectBSplinePatchFlipTrimmedDomain ( GO_BSplinePatch *obj )
  
{
  BSP_TrimmedDomain *trpd;
  int               nelem, i;
  mbs_polycurved    *elem;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( !obj->trimmed || !obj->trpd )
    return;
  trpd  = obj->trpd;
  elem  = trpd->bound;
  nelem = trpd->nelem;
  for ( i = 0; i < nelem; i++ ) {
    if ( elem[i].knots )  /* flip a B-spline curve */
      _FlipXYPoints ( elem[i].cpdimen, elem[i].lastknot-elem[i].degree,
                      elem[i].points );
    else if ( elem[i].degree > 1 ) /* flip a Bezier curve */
      _FlipXYPoints ( elem[i].cpdimen, elem[i].degree+1, elem[i].points );
    else /* flip a polyline */
      _FlipXYPoints ( elem[i].cpdimen, elem[i].lastknot+1, elem[i].points );
  }
} /*GeomObjectBSplinePatchFlipTrimmedDomain*/

