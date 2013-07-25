
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

#include "editor.h"
#include "editor_bsm.h"


boolean GeomObjectBSplineMeshSelectPoint ( GO_BSplineMesh *obj, CameraRecd *CPos,
                                           short x, short y )
{
  int dist, i;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  for ( i = BSM_NCURRENT_VERTICES-1; i > 0; i-- )
    obj->current_vertex[i] = obj->current_vertex[i-1];
  dist = MAXPIXDIST+1;
  if ( GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen, 
                         obj->nv, obj->meshvpc, obj->mkcp,
                         MASK_CP_MOVEABLE, CPos, x, y, &dist ) ) {
    return GeomObjectBSplineMeshSetCurrentVertex ( obj, current_point_ind, 0 );
  }
  else {
    GeomObjectBSplineMeshSetCurrentVertex ( obj, -1, 0 );
    return false;
  }
} /*GeomObjectBSplineMeshSelectPoint*/

boolean GeomObjectBSplineMeshSelectEdge ( GO_BSplineMesh *obj, CameraRecd *CPos,
                                          short x, short y )
{
  int         dist, dd, i, nhe, v0, v1, cpdimen, spdimen, che;
  point3d     q0, q1, qt, r;
  ray3d       ray;
  BSMhalfedge *mhe;
  double      *mvpc, s, t;

  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return false;
  CameraRayOfPixeld ( CPos, x, y, &ray );
  cpdimen = obj->me.cpdimen;
  spdimen = obj->me.spdimen;
  mvpc = obj->meshvpc;
  mhe = obj->meshhe;
  nhe = obj->nhe;
  dist = MAXPIXDIST+1;
  che = -1;
  for ( i = 0; i < nhe; i++ )
    if ( mhe[i].otherhalf < 0 || mhe[i].otherhalf > i ) {
      v0 = cpdimen*mhe[i].v0;
      v1 = cpdimen*mhe[i].v1;
        /* find the edge end points */
      switch ( cpdimen ) {
    case 2:
        SetPoint3d ( &q0, mvpc[v0], mvpc[v0+1], 0.0 );
        SetPoint3d ( &q1, mvpc[v1], mvpc[v1+1], 0.0 );
        break;
    case 3:
        if ( spdimen == 3 ) {
          memcpy ( &q0, &mvpc[v0], sizeof(point3d) );
          memcpy ( &q1, &mvpc[v1], sizeof(point3d) );
        }
        else {
          SetPoint3d ( &q0, mvpc[v0]/mvpc[v0+2], mvpc[v0+1]/mvpc[v0+2], 0.0 );
          SetPoint3d ( &q1, mvpc[v1]/mvpc[v1+2], mvpc[v1+1]/mvpc[v1+2], 0.0 );
        }
        break;
    case 4:
        Point4to3d ( (point4d*)&mvpc[v0], &q0 );
        Point4to3d ( (point4d*)&mvpc[v1], &q1 );
        break;
      }
        /* find the closest points of the ray and the edge */
      pkg_LineRayDist3d ( &ray, &q0, &q1, &s, &t, NULL, &qt );
      if ( t >= 0.0 && t <= 1.0 ) {
        CameraProjectPoint3d ( CPos, &qt, &r );
        dd = (int)(fabs(r.x-x)+fabs(r.y-y));
        if ( dd < dist ) {
          dist = dd;
          che = i;
        }
      }
    }
  for ( i = BSM_NCURRENT_EDGES-1; i > 0; i-- )
    obj->current_edge[i] = obj->current_edge[i-1];
  if ( che >= 0 ) {
    GeomObjectBSplineMeshSetCurrentEdge ( obj, che, 0 );
    i = mhe[che].otherhalf;
    if ( i >= 0 )
      GeomObjectBSplineMeshSetCurrentEdge ( obj, i, 1 );
  }
  return true;
} /*GeomObjectBSplineMeshSelectEdge*/

void GeomObjectBSplineMeshUnselectPoint ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type == GO_BSPLINE_MESH ) {
    GeomObjectBSplineMeshSetCurrentVertex ( obj, -1, 0 );
    GeomObjectBSplineMeshSetCurrentVertex ( obj, -1, 1 );
  }
} /*GeomObjectBSplineMeshUnselectPoint*/

void GeomObjectBSplineMeshUnselectEdge ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type == GO_BSPLINE_MESH ) {
    GeomObjectBSplineMeshSetCurrentEdge ( obj, -1, 0 );
    GeomObjectBSplineMeshSetCurrentEdge ( obj, -1, 1 );
  }
} /*GeomObjectBSplineMeshUnselectEdge*/

