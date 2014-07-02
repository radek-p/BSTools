
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
#include "mengerc.h"
#include "bsfile.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

static void GeomObjectGetCPoints ( geom_object *go, int *cpdimen, int *spdimen,
                                   int *ncp, double **cpoints, byte **mkcp )
{
  GO_BezierCurve  *bc;
  GO_BezierPatch  *bp;
  GO_BSplineCurve *bsc;
  GO_BSplinePatch *bsp;
  GO_BSplineMesh  *bsm;

  *cpdimen = go->cpdimen;
  *spdimen = go->spdimen;
  switch ( go->obj_type ) {
case GO_BEZIER_CURVE:
    bc = (GO_BezierCurve*)go;
    *ncp = bc->degree+1;
    *cpoints = bc->cpoints;
    *mkcp = bc->mkcp;
    break;
case GO_BEZIER_PATCH:
    bp = (GO_BezierPatch*)go;
    *ncp = (bp->degree_u+1)*(bp->degree_v+1);
    *cpoints = bp->cpoints;
    *mkcp = bp->mkcp;
    break;
case GO_BSPLINE_CURVE:
    bsc = (GO_BSplineCurve*)go;
    *ncp = bsc->lastknot - bsc->degree;
    *cpoints = bsc->cpoints;
    *mkcp = bsc->mkcp;
    break;
case GO_BSPLINE_PATCH:
    bsp = (GO_BSplinePatch*)go;
    *ncp = (bsp->lastknot_u - bsp->degree_u)*(bsp->lastknot_v - bsp->degree_v);
    *cpoints = bsp->cpoints;
    *mkcp = bsp->mkcp;
    break;
case GO_BSPLINE_MESH:
    bsm = (GO_BSplineMesh*)go;
    *ncp = bsm->nv;
    *cpoints = bsm->meshvpc;
    *mkcp = bsm->mkcp;
    break;
case GO_BSPLINE_HOLE:
default:
    *ncp = 0;
    *cpoints = NULL;
    *mkcp = NULL;
    break;
  }
} /*GeomObjectGetCPoints*/

void GeomObjectSpecial3DTranslate ( vector3d *v )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp;
  byte        *mkcp;
  int         i;
  point3d     *q, qq;
  point4d     *Q;

  GeomObjectSaveCPoints ( 3 );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              AddVector3d ( q, v, q );
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              AddVector3d ( &qq, v, &qq );
              Point3to4d ( &qq, Q->w, Q );
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectSpecial3DTranslate*/

void GeomObjectSpecial3DScale ( point3d *p, vector3d *s )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp;
  byte        *mkcp;
  int         i;
  point3d     *q, qq;
  point4d     *Q;

  GeomObjectSaveCPoints ( 3 );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              SubtractPoints3d ( q, p, q );
              q->x *= s->x;
              q->y *= s->y;
              q->z *= s->z;
              AddVector3d ( q, p, q );
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              SubtractPoints3d ( &qq, p, &qq );
              qq.x *= s->x;
              qq.y *= s->y;
              qq.z *= s->z;
              AddVector3d ( &qq, p, &qq );
              Point3to4d ( &qq, Q->w, Q );
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectSpecial3DScale*/

void GeomObjectSpecial3DRotate ( point3d *p, vector3d *v, double phi )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp;
  byte        *mkcp;
  int         i;
  point3d     *q;
  point4d     *Q;
  trans3d     tr;

  GeomObjectSaveCPoints ( 3 );
  IdentTrans3d ( &tr );
  ShiftTrans3d ( &tr, -p->x, -p->y, -p->z );
  RotVTrans3d ( &tr, v, phi );
  ShiftTrans3d ( &tr, p->x, p->y, p->z );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              TransPoint3d ( &tr, q, q );
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Trans3Point4d ( &tr, Q, Q );
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectSpecial3DRotate*/

void GeomObjectProject3DPointsOnALine ( point3d *p, vector3d *v )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp;
  byte        *mkcp;
  int         i;
  vector3d    vv, w;
  point3d     *q, qq;
  point4d     *Q;
  double      a;

  GeomObjectSaveCPoints ( 3 );
  vv = *v;
  NormalizeVector3d ( &vv );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              SubtractPoints3d ( q, p, &w );
              a = DotProduct3d ( &vv, &w );
              AddVector3Md ( p, &vv, a, q );
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              SubtractPoints3d ( &qq, p, &w );
              a = DotProduct3d ( &vv, &w );
              AddVector3Md ( p, &vv, a, &qq );
              Point3to4d ( &qq, Q->w, Q );
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectProject3DPointsOnALine*/

void GeomObjectProject3DPointsOnAPlane ( point3d *p, vector3d *v, char dir )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp;
  byte        *mkcp;
  int         i;
  vector3d    vv, w;
  point3d     *q, qq;
  point4d     *Q;
  double      a;
  boolean     doit;

  GeomObjectSaveCPoints ( 3 );
  vv = *v;
  NormalizeVector3d ( &vv );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              SubtractPoints3d ( q, p, &w );
              a = DotProduct3d ( &vv, &w );
              switch ( dir ) {
            case 1:  doit = a > 0;  break;
            case 2:  doit = a < 0;  break;
            default: doit = true;   break;
              }
              if ( doit )
                AddVector3Md ( q, &vv, -a, q );
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              SubtractPoints3d ( &qq, p, &w );
              a = DotProduct3d ( &vv, &w );
              switch ( dir ) {
            case 1:  doit = a > 0.0;  break;
            case 2:  doit = a < 0.0;  break;
            default: doit = true;     break;
              }
              if ( doit ) {
                AddVector3Md ( &qq, &vv, -a, &qq );
                Point3to4d ( &qq, Q->w, Q );
              }
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectProject3DPointsOnAPlane*/

void GeomObjectProject3DPointsOnASphere ( point3d *p, double r, char dir )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp, rr;
  byte        *mkcp;
  int         i;
  vector3d    w;
  point3d     *q, qq;
  point4d     *Q;
  boolean     doit;

  GeomObjectSaveCPoints ( 3 );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              SubtractPoints3d ( q, p, &w );
              rr = sqrt ( DotProduct3d ( &w, &w ) );
              switch ( dir ) {
            case 1:  doit = rr > r;  break;
            case 2:  doit = rr < r;  break;
            default: doit = true;    break;
              }
              if ( doit ) {
                NormalizeVector3d ( &w );
                AddVector3Md ( p, &w, r, q );
              }
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              SubtractPoints3d ( &qq, p, &w );
              rr = sqrt ( DotProduct3d ( &w, &w ) );
              switch ( dir ) {
            case 1:  doit = rr > r;  break;
            case 2:  doit = rr < r;  break;
            default: doit = true;    break;
              }
              if ( doit ) {
                NormalizeVector3d ( &w );
                AddVector3Md ( p, &w, r, &qq );
                Point3to4d ( &qq, Q->w, Q );
              }
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectProject3DPointsOnASphere*/

void GeomObjectProject3DPointsOnACylinder ( point3d *p, vector3d *v, double r,
                                            char dir )
{
  geom_object *go;
  int         cpdimen, spdimen, ncp;
  double      *cp, rr;
  byte        *mkcp;
  int         i;
  vector3d    vv, w, z;
  point3d     *q, qq;
  point4d     *Q;
  double      a;
  boolean     doit;

  GeomObjectSaveCPoints ( 3 );
  vv = *v;
  NormalizeVector3d ( &vv );
  for ( go = first_go; go; go = go->next )
    if ( go->active || go == current_go ) {
      GeomObjectGetCPoints ( go, &cpdimen, &spdimen, &ncp, &cp, &mkcp );
      if ( spdimen == 3 && ncp > 0 ) {
        for ( i = 0; i < ncp; i++ )
          if ( mkcp[i] & marking_mask ) {
            switch ( cpdimen ) {
          case 3:
              q = (point3d*)&cp[3*i];
              SubtractPoints3d ( q, p, &w );
              a = DotProduct3d ( &vv, &w );
              AddVector3Md ( &w, &vv, -a, &z );
              rr = sqrt ( DotProduct3d ( &z, &z ) );
              switch ( dir ) {
            case 1:  doit = rr > r;  break;
            case 2:  doit = rr < r;  break;
            default: doit = true;    break;
              }
              if ( doit ) {
                NormalizeVector3d ( &z );
                AddVector3Md ( p, &vv, a, q );
                AddVector3Md ( q, &z, r, q );
              }
              break;
          case 4:
              Q = (point4d*)&cp[4*i];
              Point4to3d ( Q, &qq );
              SubtractPoints3d ( &qq, p, &w );
              a = DotProduct3d ( &vv, &w );
              AddVector3Md ( &w, &vv, -a, &z );
              rr = sqrt ( DotProduct3d ( &z, &z ) );
              switch ( dir ) {
            case 1:  doit = rr > r;  break;
            case 2:  doit = rr < r;  break;
            default: doit = true;    break;
              }
              if ( doit ) {
                NormalizeVector3d ( &z );
                AddVector3Md ( p, &vv, a, &qq );
                AddVector3Md ( &qq, &z, r, &qq );
                Point3to4d ( &qq, Q->w, Q );
              }
              break;
          default:
              return;
            }
          }
        go->dlistmask = 0;
      }
    }
} /*GeomObjectProject3DPointsOnACylinder*/

