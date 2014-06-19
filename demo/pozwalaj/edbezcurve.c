
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
#include "editor_bezc.h"


boolean GeomObjectInitBezierCurve ( GO_BezierCurve *obj,
                                    char spdimen, boolean rational )
{
  static const double inicp[8] = {-1.0,0.0,0.0,1.0, 1.0,0.0,0.0,1.0};

  obj->me.obj_type = GO_BEZIER_CURVE;
  obj->me.ident = -1;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational ) {
    obj->me.cpdimen = spdimen+1;
    obj->weightpoints = malloc ( obj->me.cpdimen*sizeof(double) );
    if ( !obj->weightpoints )
      return false;
  }
  else {
    obj->me.cpdimen = spdimen;
    obj->weightpoints = NULL;
  }
  obj->me.active = false;
  obj->me.name[0] = 0;
  obj->savedsize = 0;
  obj->cpoints = malloc ( 2*obj->me.cpdimen*sizeof(double) );
  obj->mkcp = malloc ( 2 );
  if ( !obj->cpoints || !obj->mkcp ) {
    if ( obj->cpoints )      free ( obj->cpoints );
    if ( obj->weightpoints ) free ( obj->weightpoints );
    if ( obj->mkcp )         free ( obj->mkcp );
    return false;
  }
  GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                             2, inicp, obj->cpoints );
  if ( obj->rational )
    GeomObjectSetupWeightPoints ( obj->me.cpdimen, 2, obj->cpoints,
                                  obj->weightpoints );
  memset ( obj->mkcp, MASK_CP_MOVEABLE, 2 );
  obj->degree = 1;
  obj->view_curve = obj->view_cpoly = true;
  obj->view_curvature = obj->view_torsion = false;
  obj->graph_dens = 10;
  obj->me.displaylist = glGenLists ( BEZC_NDL );
  obj->me.dlistmask = 0;
  obj->me.colour[0] = obj->me.colour[1] = 1.0;
  obj->me.colour[2] = 0.0;
  obj->me.display_pretrans = false;
  IdentTrans3d ( &obj->me.pretrans );
  obj->pipe_diameter = 0.1;
  return true;
} /*GeomObjectInitBezierCurve*/

geom_object *GeomObjectAddBezierCurve ( const char *name, char spdimen,
                                        boolean rational )
{
  GO_BezierCurve *obj;

  obj = malloc ( sizeof(GO_BezierCurve) );
  if ( obj ) {
    memset ( obj, 0, sizeof(GO_BezierCurve) );
    if ( !GeomObjectInitBezierCurve ( obj, spdimen, rational ) ) {
      free ( obj );
      return NULL;
    }
    strncpy ( obj->me.name, name, MAX_NAME_LENGTH+1 );
    if ( !first_go )
      first_go = last_go = &obj->me;
    else {
      obj->me.prev = last_go;
      last_go->next = &obj->me;
      last_go = &obj->me;
    }
    current_go = &obj->me;
    return &obj->me;
  }
  else
    return NULL;
} /*GeomObjectAddBezierCurve*/

geom_object *GeomObjectCopyBezierCurve ( GO_BezierCurve *obj )
{
  GO_BezierCurve *copy;
  int            ncp;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return NULL;
  copy = malloc ( sizeof(GO_BezierCurve) );
  if ( copy ) {
    ncp = obj->degree+1;
    memset ( copy, 0, sizeof(GO_BezierCurve) );
    copy->me.obj_type = GO_BEZIER_CURVE;
    copy->me.ident = -1;
    copy->me.spdimen = obj->me.spdimen;
    copy->me.cpdimen = obj->me.cpdimen;
    copy->me.active = false;
    strcpy ( copy->me.name, obj->me.name );
    copy->me.display_pretrans = false;
    copy->me.pretrans = obj->me.pretrans;
    if ( obj->rational ) {
      copy->weightpoints = malloc ( (ncp-1)*obj->me.cpdimen*sizeof(double) );
      if ( !copy->weightpoints ) {
        free ( copy );
        return NULL;
      }
    }
    else
      copy->weightpoints = NULL;
    copy->cpoints = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    copy->mkcp = malloc ( ncp );
    if ( !copy->cpoints || !copy->mkcp ) {
      if ( copy->cpoints )      free ( copy->cpoints );
      if ( copy->mkcp )         free ( copy->mkcp );
      if ( copy->weightpoints ) free ( copy->weightpoints );
      free ( copy );
      return NULL;
    }
    copy->degree = obj->degree;
    memcpy ( copy->cpoints, obj->cpoints,
             ncp*obj->me.cpdimen*sizeof(double) );
    memset ( copy->mkcp, MASK_CP_MOVEABLE, ncp );
    copy->rational = obj->rational;
    if ( copy->rational )
      GeomObjectSetupWeightPoints ( copy->me.cpdimen, ncp,
                                    copy->cpoints, copy->weightpoints );
    copy->view_curve = copy->view_cpoly = true;
    copy->view_curvature = copy->view_torsion = false;
    copy->graph_dens = 10;
    copy->pipe_diameter = obj->pipe_diameter;
    copy->me.displaylist = glGenLists ( BEZC_NDL );
    copy->me.dlistmask = 0;
    return &copy->me;
  }
  else
    return NULL;
} /*GeomObjectCopyBezierCurve*/

void GeomObjectDeleteBezierCurve ( GO_BezierCurve *obj )
{
  glDeleteLists ( obj->me.displaylist, BEZC_NDL );
  if ( obj->cpoints )      free ( obj->cpoints );
  if ( obj->mkcp )         free ( obj->mkcp );
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  if ( obj->weightpoints ) free ( obj->weightpoints );
  free ( obj );
} /*GeomObjectDeleteBezierCurve*/

void GeomObjectDrawBezierCurve ( GO_BezierCurve *obj )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->me.dlistmask & BEZC_DLM_CURVE )
    glCallList ( obj->me.displaylist+1 );
  else {
    glNewList ( obj->me.displaylist+1, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_White );
    if ( obj->degree == 1 )
      DrawAPolyline ( obj->me.cpdimen, obj->me.spdimen, 2, obj->cpoints );
    else
      DrawBezierCurve ( obj->me.cpdimen, obj->me.spdimen, obj->degree, obj->cpoints,
                        8*obj->degree );
    glEndList ();
    obj->me.dlistmask |= BEZC_DLM_CURVE;
  }
} /*GeomObjectDrawBezierCurve*/

void GeomObjectDrawBezierCPoly ( GO_BezierCurve *obj )
{
  int ncp;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->me.dlistmask & BEZC_DLM_CPOLY )
    glCallList ( obj->me.displaylist );
  else {
    ncp = obj->degree+1;
    glNewList ( obj->me.displaylist, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_Green );
    if ( obj->rational && obj->weightpoints ) {
      DrawLineSegments ( obj->me.cpdimen, obj->me.spdimen, ncp-1,
                         obj->cpoints, obj->weightpoints );
      glColor3fv ( xglec_Green4 );
      DrawLineSegments ( obj->me.cpdimen, obj->me.spdimen, ncp-1,
                         obj->weightpoints, &obj->cpoints[(int)obj->me.cpdimen] );
    }
    else
      DrawAPolyline ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints );
    glColor3fv ( xglec_Yellow );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, 0, MASK_CP_MOVEABLE );
    glColor3fv ( xglec_OrangeRed );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
    if ( obj->rational ) {
      glColor3fv ( xglec_DarkOliveGreen2 );
      DrawPoints ( obj->me.cpdimen, obj->me.spdimen, ncp-1, obj->weightpoints );
    }
    glEndList ();
    obj->me.dlistmask |= BEZC_DLM_CPOLY;
  }
} /*GeomObjectDrawBezierCPoly*/

void GeomObjectDrawBezierCurvature ( GO_BezierCurve *obj )
{
  int      dim, deg, dens, j;
  double   t;
  point2d  p2, q2;
  vector2d f2[2];
  point3d  p3, q3;
  vector3d f3[3];
  double   curv[2], scf;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->me.dlistmask & BEZC_DLM_CURVATURE )
    glCallList ( obj->me.displaylist+2 );
  else {
    dim = obj->me.cpdimen;
    deg = obj->degree;
    dens = obj->graph_dens;
    scf = xge_LogSlidebarValued ( 0.01, 100.0, obj->curvature_scale );
    glNewList ( obj->me.displaylist+2, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_CornflowerBlue );
    glBegin ( GL_LINES );
    switch ( obj->me.spdimen ) {
  case 2:
      if ( dim == 2 ) {      /* non-rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC2d ( deg, (point2d*)obj->cpoints, t, &p2, f2, curv );
          AddVector2Md ( &p2, &f2[1], curv[0]*scf, &q2 );
          glVertex2dv ( &p2.x );
          glVertex2dv ( &q2.x );
        }
      }
      else if ( dim == 3 ) { /* rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC2Rd ( deg, (point3d*)obj->cpoints, t, &p2, f2, curv );
          AddVector2Md ( &p2, &f2[1], curv[0]*scf, &q2 );
          glVertex2dv ( &p2.x );
          glVertex2dv ( &q2.x );
        }
      }
      break;
  case 3:
      if ( dim == 3 ) {      /* non-rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC3d ( deg, (point3d*)obj->cpoints, t, &p3, f3, curv );
          AddVector3Md ( &p3, &f3[1], curv[0]*scf, &q3 );
          glVertex3dv ( &p3.x );
          glVertex3dv ( &q3.x );
        }
      }
      else if ( dim == 3 ) { /* rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC3Rd ( deg, (point4d*)obj->cpoints, t, &p3, f3, curv );
          AddVector3Md ( &p3, &f3[1], curv[0]*scf, &q3 );
          glVertex3dv ( &p3.x );
          glVertex3dv ( &q3.x );
        }
      }
      break;
  default:
      break;
    }
    glEnd ();
    glEndList ();
    obj->me.dlistmask |= BEZC_DLM_CURVATURE;
  }
} /*GeomObjectDrawBezierCurvature*/

void GeomObjectDrawBezierTorsion ( GO_BezierCurve *obj )
{
  int      dim, deg, j, dens;
  double   t;
  point3d  p3, q3;
  vector3d f3[3];
  double   curv[2], scf;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->me.dlistmask & BEZC_DLM_TORSION )
    glCallList ( obj->me.displaylist+3 );
  else {
    dim = obj->me.cpdimen;
    deg = obj->degree;
    dens = obj->graph_dens;
    scf = xge_LogSlidebarValued ( 0.01, 100.0, obj->torsion_scale );
    glNewList ( obj->me.displaylist+3, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( xglec_SpringGreen2 );
    glBegin ( GL_LINES );
    switch ( obj->me.spdimen ) {
  case 3:
      if ( dim == 3 ) {      /* non-rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC3d ( deg, (point3d*)obj->cpoints, t, &p3, f3, curv );
          AddVector3Md ( &p3, &f3[2], curv[1]*scf, &q3 );
          glVertex3dv ( &p3.x );
          glVertex3dv ( &q3.x );
        }
      }
      else if ( dim == 3 ) { /* rational */
        for ( j = 0; j <= dens; j++ ) {
          t = (double)j/(double)dens;
          mbs_BCFrenetC3Rd ( deg, (point4d*)obj->cpoints, t, &p3, f3, curv );
          AddVector3Md ( &p3, &f3[2], curv[1]*scf, &q3 );
          glVertex3dv ( &p3.x );
          glVertex3dv ( &q3.x );
        }
      }
      break;
  default:
      break;
    }
    glEnd ();
    glEndList ();
    obj->me.dlistmask |= BEZC_DLM_TORSION;
  }
} /*GeomObjectDrawBezierTorsion*/

void GeomObjectDisplayBezierCurve ( GO_BezierCurve *obj )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->view_curve )
    GeomObjectDrawBezierCurve ( obj );
  if ( obj->view_cpoly )
    GeomObjectDrawBezierCPoly ( obj );
  if ( obj->view_curvature )
    GeomObjectDrawBezierCurvature ( obj );
  if ( obj->view_torsion )
    GeomObjectDrawBezierTorsion ( obj );
} /*GeomObjectDisplayBezierCurve*/

boolean GeomObjectBezierCurveSetDegree ( GO_BezierCurve *obj, int deg )
{
  int    ncp;
  double *cp, *wcp;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return false;
  if ( deg < 1 || deg > MAX_DEGREE )
    return false;
  if ( obj->rational ) {
    wcp = malloc ( deg*obj->me.cpdimen*sizeof(double) );
    if ( !wcp )
      return false;
  }
  else
    wcp = NULL;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  ncp = deg+1;
  cp = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !cp || !mkcp ) {
    goto failure;
  }
  if ( deg > obj->degree ) {
    if ( !mbs_multiBCDegElevd ( 1, obj->me.cpdimen, 0, obj->degree, obj->cpoints,
                                deg-obj->degree, 0, &deg, cp ) )
      goto failure;
  }
  else if ( deg < obj->degree )
    mbs_multiBCDegRedd ( 1, obj->me.cpdimen, 0, obj->degree, obj->cpoints,
                         obj->degree-deg, 0, &deg, cp );
  else
    memcpy ( cp, obj->cpoints, ncp*obj->me.cpdimen*sizeof(double) );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  if ( obj->rational ) {
    free ( obj->weightpoints );
    obj->weightpoints = wcp;
    GeomObjectSetupWeightPoints ( obj->me.cpdimen, deg+1, cp, wcp );
  }
  obj->degree = deg;
  obj->me.dlistmask = 0;
  GeomObjectBezierCurveDisplayInfoText ( obj );
  return true;

failure:
  if ( cp )   free ( cp );
  if ( mkcp ) free ( mkcp );
  if ( wcp )  free ( wcp );
  return false;
} /*GeomObjectBezierCurveSetDegree*/

void GeomObjectBezierCurveFindBBox ( GO_BezierCurve *obj, Box3d *box )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  GeomObjectFindBBox ( obj->me.cpdimen, obj->rational,
                       obj->degree+1, obj->cpoints, box );
} /*GeomObjectBezierCurveFindBBox*/

boolean GeomObjectBezierCurveFindCPoint ( GO_BezierCurve *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist )
{
  boolean ok, ok1;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return false;
  ok = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
                                    obj->degree+1, obj->cpoints,
                                    obj->mkcp, MASK_CP_MOVEABLE, CPos, x, y, dist );
  if ( obj->rational ) {
    ok1 = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
                                       obj->degree, obj->weightpoints,
                                       NULL, 0, CPos, x, y, dist );
    if ( ok1 )
      current_point_ind += obj->degree+1;
  }
  else
    ok1 = false;
  return ok || ok1;
} /*GeomObjectBezierCurveFindCPoint*/

void GeomObjectBezierCurveSetCPoint ( GO_BezierCurve *obj,
                                      CameraRecd *CPos, short x, short y )
{
  int k, c, d;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  d = obj->me.cpdimen;
  if ( current_point_ind >= 0 && current_point_ind <= obj->degree ) {
    k = d*current_point_ind;
    GeomObjectSetPoint ( CPos, x, y, d, obj->me.spdimen,
                         &obj->cpoints[k] );
    if ( obj->rational ) {
      if ( current_point_ind == 0 )
        GeomObjectSetupWeightPoints ( d, 2, obj->cpoints, obj->weightpoints );
      else if ( current_point_ind < obj->degree )
        GeomObjectSetupWeightPoints ( d, 3, &obj->cpoints[k-d],
                                      &obj->weightpoints[k-d] );
      else
        GeomObjectSetupWeightPoints ( d, 2, &obj->cpoints[k-d],
                                      &obj->weightpoints[k-d] );
    }
  }
  else {  /* a Farin point has been moved */
    c = current_point_ind-(obj->degree+1);
    k = d*c;
    GeomObjectSetWeightPoint ( CPos, x, y, d, &obj->cpoints[k],
                               &obj->weightpoints[k] );
    if ( c < obj->degree-1 )
      GeomObjectSetupWeightPoints ( d, 2, &obj->cpoints[k+d],
                                    &obj->weightpoints[k+d] );
  }
  obj->me.dlistmask = 0;
} /*GeomObjectBezierCurveSetCPoint*/

void GeomObjectBezierCurveMarkCPoints ( GO_BezierCurve *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        char mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen, obj->degree+1,
                         obj->mkcp, obj->cpoints, CPos, box, mask, clear );
  obj->me.dlistmask &= ~BEZC_DLM_CPOLY;
} /*GeomObjectBezierCurveMarkCPoints*/

void GeomObjectBezierCurveMarkCPoint ( GO_BezierCurve *obj,
                                       char mask, boolean clear )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  GeomObjectMarkPoint ( obj->degree+1, obj->mkcp, mask, clear );
  obj->me.dlistmask &= ~BEZC_DLM_CPOLY;
} /*GeomObjectBezierCurveMarkCPoint*/

boolean GeomObjectBezierCurveSaveCPoints ( GO_BezierCurve *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return false;
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  size = (obj->degree+1)*(obj->me.cpdimen)*sizeof(double);
  obj->savedcpoints = malloc ( size );
  if ( !obj->savedcpoints )
    return false;
  memcpy ( obj->savedcpoints, obj->cpoints, size );
  obj->savedsize = size;
  return true;
} /*GeomObjectBezierCurveSaveCPoints*/

void GeomObjectBezierCurveUndoLastTransformation ( GO_BezierCurve *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  size = (obj->degree+1)*(obj->me.cpdimen)*sizeof(double);
  if ( obj->savedcpoints && obj->savedsize == size ) {
    memcpy ( obj->cpoints, obj->savedcpoints, size );
    if ( obj->rational )
      GeomObjectSetupWeightPoints ( obj->me.cpdimen, obj->degree+1,
                                    obj->cpoints, obj->weightpoints );
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBezierCurveUndoLastTransformation*/

void GeomObjectBezierCurveTransformCPoints ( GO_BezierCurve *obj,
                                             trans3d *tr, char mask )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->savedcpoints ) {
    GeomObjectTransformPoints ( obj->me.cpdimen, obj->me.spdimen, obj->degree+1,
                                obj->mkcp, mask,
                                obj->savedcpoints, obj->cpoints, tr );
    obj->me.dlistmask = 0;
    if ( obj->rational )
      GeomObjectSetupWeightPoints ( obj->me.cpdimen, obj->degree+1,
                                    obj->cpoints, obj->weightpoints );
  }
} /*GeomObjectBezierCurveTransformCPoints*/

void GeomObjectBezierCurveSetCurvatureGraph ( GO_BezierCurve *obj,
                                boolean view_curvature, double curvature_scale,
                                boolean view_torsion, double torsion_scale,
                                int graph_dens )
{
  if ( (view_curvature && !obj->view_curvature) ||
       curvature_scale != obj->curvature_scale ||
       graph_dens != obj->graph_dens )
    obj->me.dlistmask &= ~BEZC_DLM_CURVATURE;
  obj->view_curvature = view_curvature;
  if ( (view_torsion && !obj->view_torsion) ||
       torsion_scale != obj->torsion_scale ||
       graph_dens != obj->graph_dens )
    obj->me.dlistmask &= ~BEZC_DLM_TORSION;
  obj->curvature_scale = curvature_scale;
  obj->view_torsion = view_torsion && (obj->me.spdimen == 3);
  obj->torsion_scale = torsion_scale;
  obj->graph_dens = graph_dens;
} /*GeomObjectBezierCurveSetCurvatureGraph*/

boolean GeomObjectBezierCurveGetPointCoord ( GO_BezierCurve *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return false;
  if ( p < 0 || p > obj->degree )
    return false;
  *spdimen = obj->me.spdimen;
  *cpdimen = obj->me.cpdimen;
  *pc = &obj->cpoints[obj->me.cpdimen*p];
  return true;
} /*GeomObjectBezierCurveGetPointCoord*/

boolean GeomObjectWriteBCAttributes ( GO_BezierCurve *obj )
{
  return true;
} /*GeomObjectWriteBCAttributes*/

boolean GeomObjectWriteBezierCurve ( GO_BezierCurve *obj )
{
  return bsf_WriteBezierCurved ( obj->me.spdimen, obj->me.cpdimen,
                                 obj->rational, obj->degree,
                                 obj->cpoints, obj->me.name, obj->me.ident,
                                 GeomObjectWriteAttributes,
                                 (void*)&obj->me );
} /*GeomObjectWriteBezierCurve*/

boolean GeomObjectBCResolveDependencies ( GO_BezierCurve *obj )
{
  return true;
} /*GeomObjectBCResolveDependencies*/

void GeomObjectReadBezierCurve ( void *usrdata,
                                 const char *name, int ident, int degree,
                                 const point4d *cpoints, int spdimen,
                                 boolean rational )
{
  GO_BezierCurve *obj;
  double         *cp, *wcp;
  byte           *mkcp;
  int            ncp, cpdimen;

  obj = (GO_BezierCurve*)GeomObjectAddBezierCurve ( name, spdimen, rational );
  if ( obj ) {
    obj->me.ident = ident;
    ncp = degree+1;
    cpdimen = obj->me.cpdimen;
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( rational )
      wcp = malloc ( degree*cpdimen*sizeof(double) );
    else
      wcp = NULL;
    if ( !cp || !mkcp || (rational && !wcp ) ) {
      if ( mkcp ) free ( mkcp );
      if ( cp )   free ( cp );
      if ( wcp )  free ( wcp );
      GeomObjectDeleteBezierCurve ( obj );
      return;
    }
    if ( obj->weightpoints ) {
      free ( obj->weightpoints );
      obj->weightpoints = NULL;
    }
    free ( obj->cpoints );
    free ( obj->mkcp );
    if ( obj->savedcpoints ) {
      free ( obj->savedcpoints );
      obj->savedcpoints = NULL;
    }
    obj->cpoints = cp;
    obj->weightpoints = wcp;
    obj->mkcp = mkcp;
    obj->degree = degree;
    GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                               ncp, (double*)cpoints, cp );
    if ( obj->rational )
      GeomObjectSetupWeightPoints ( obj->me.cpdimen, obj->degree+1,
                                    obj->cpoints, obj->weightpoints );
  }
} /*GeomObjectReadBezierCurve*/

void GeomObjectBezierCurveOutputToRenderer ( GO_BezierCurve *obj )
{
  if ( obj->me.obj_type != GO_BEZIER_CURVE )
    return;
  if ( obj->me.spdimen != 3 )
    return;
  if ( obj->rational )
    RendEnterBezCurve3Rd ( obj->degree, (point4d*)obj->cpoints,
                           0.5*obj->pipe_diameter, obj->me.colour );
  else
    RendEnterBezCurve3d ( obj->degree, (point3d*)obj->cpoints,
                          0.5*obj->pipe_diameter, obj->me.colour );
} /*GeomObjectBezierCurveOutputToRenderer*/

void GeomObjectBezierCurveDisplayInfoText ( GO_BezierCurve *obj )
{
  char s[80];

  sprintf ( s, "deg = %d", obj->degree );
  SetStatusText ( s, true );
} /*GeomObjectBezierCurveDisplayInfoText*/

boolean GeomObjectBezierCurveProcessDep ( GO_BezierCurve *obj, geom_object *go )
{
  return false;
} /*GeomObjectBezierCurveProcessDep*/

void GeomObjectBezierCurveProcessDeletedDep ( GO_BezierCurve *obj, geom_object *go )
{
} /*GeomObjectBezierCurveProcessDeletedDep*/

