
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
#include "pkrender.h"
#include "xgedit.h"
#include "xgledit.h"

#include "editor.h"
#include "edcolours.h"
#include "editor_bezp.h"

extern pkRenderer rend;

boolean GeomObjectInitBezierPatch ( GO_BezierPatch *obj,
                                    char spdimen, boolean rational )
{
  static const double inicp[16] =
    {-1.0,-1.0,0.0,1.0, -1.0,1.0,0.0,1.0,
      1.0,-1.0,0.0,1.0,  1.0,1.0,0.0,1.0};

  obj->me.obj_type = GO_BEZIER_PATCH;
  obj->me.ident = -1;
  obj->me.spdimen = spdimen;
  obj->rational = rational;
  if ( rational )
    obj->me.cpdimen = spdimen+1;
  else
    obj->me.cpdimen = spdimen;
  obj->me.active = false;
  obj->me.name[0] = 0;
  obj->savedsize = 0;
  obj->cpoints = malloc ( 4*obj->me.cpdimen*sizeof(double) );
  obj->mkcp = malloc ( 4 );
  if ( !obj->cpoints || !obj->mkcp ) {
    if ( obj->cpoints ) free ( obj->cpoints );
    if ( obj->mkcp ) free ( obj->mkcp );
    return false;
  }
  GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                             4, inicp, obj->cpoints );
  memset ( obj->mkcp, MASK_CP_MOVEABLE, 4 );
  obj->degree_u = obj->degree_v = 1;
  obj->dens_u = obj->dens_v = 8;
  obj->view_surf = obj->view_cnet = true;
  obj->me.displaylist = glGenLists ( BEZP_NDL );
  obj->me.dlistmask = 0;
  obj->me.colour[0] = obj->me.colour[1] = 1.0;
  obj->me.colour[2] = 0.0;
  obj->me.display_pretrans = false;
  IdentTrans3d ( &obj->me.pretrans );
  return true;
} /*GeomObjectInitBezierPatch*/

geom_object *GeomObjectAddBezierPatch ( const char *name,
                                        char spdimen, boolean rational )
{
  GO_BezierPatch *obj;

  obj = malloc ( sizeof(GO_BezierPatch) );
  if ( obj ) {
    memset ( obj, 0, sizeof(GO_BezierPatch) );
    if ( !GeomObjectInitBezierPatch ( obj, spdimen, rational ) ) {
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
} /*GeomObjectAddBezierPatch*/

geom_object *GeomObjectCopyBezierPatch ( GO_BezierPatch *obj )
{
  GO_BezierPatch *copy;
  int            ncp;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return NULL;
  copy = malloc ( sizeof(GO_BezierPatch) );
  if ( copy ) {
    memset ( copy, 0, sizeof(GO_BezierPatch) );
    ncp = (obj->degree_u+1)*(obj->degree_v+1);
    copy->me.obj_type = GO_BEZIER_PATCH;
    copy->me.ident = -1;
    copy->me.spdimen = obj->me.spdimen;
    copy->me.cpdimen = obj->me.cpdimen;
    copy->me.active = false;
    strcpy ( copy->me.name, obj->me.name );
    copy->me.display_pretrans = false;
    copy->me.pretrans = obj->me.pretrans;
    copy->cpoints = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
    copy->mkcp = malloc ( ncp );
    if ( !copy->cpoints || !copy->mkcp ) {
      if ( copy->cpoints ) free ( copy->cpoints );
      if ( copy->mkcp ) free ( copy->mkcp );
      free ( copy );
      return NULL;
    }
    copy->degree_u = obj->degree_u;
    copy->degree_v = obj->degree_v;
    memcpy ( copy->cpoints, obj->cpoints, ncp*obj->me.cpdimen*sizeof(double) );
    memset ( copy->mkcp, MASK_CP_MOVEABLE, ncp );
    copy->rational = obj->rational;
    copy->dens_u = obj->dens_u;
    copy->dens_v = obj->dens_v;
    copy->view_surf = copy->view_cnet = true;
    copy->me.displaylist = glGenLists ( BEZP_NDL );
    copy->me.dlistmask = 0;
    return &copy->me;
  }
  else
    return NULL;
} /*GeomObjectCopyBezierPatch*/

void GeomObjectDeleteBezierPatch ( GO_BezierPatch *obj )
{
  glDeleteLists ( obj->me.displaylist, BEZP_NDL );
  if ( obj->cpoints ) free ( obj->cpoints );
  if ( obj->mkcp ) free ( obj->mkcp );
  if ( obj->savedcpoints) free ( obj->savedcpoints );
  free ( obj );
} /*GeomObjectDeleteBezierPatch*/

void GeomObjectDrawBezierPatch ( GO_BezierPatch *obj )
{
  int pitch;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  if ( obj->me.dlistmask & BEZP_DLM_PATCH )
    glCallList ( obj->me.displaylist );
  else {
    glNewList ( obj->me.displaylist, GL_COMPILE_AND_EXECUTE );
    pitch = obj->me.cpdimen*(obj->degree_v+1);
    glColor3fv ( OBJC_BEZP_PATCH_BOUNDARY );
    DrawBezierPatchWF ( obj->me.cpdimen, obj->me.spdimen,
                        obj->degree_u, obj->degree_v,
                        pitch, obj->cpoints, 1, 1,
                        8*obj->dens_u, 8*obj->dens_v,
                        true, true, true, true );
    glColor3fv ( OBJC_BEZP_PATCH_CPLINES );
    DrawBezierPatchWF ( obj->me.cpdimen, obj->me.spdimen,
                        obj->degree_u, obj->degree_v,
                        pitch, obj->cpoints, obj->dens_u, obj->dens_v,
                        8*obj->dens_u, 8*obj->dens_v,
                        false, false, false, false );
    glEndList ();
    obj->me.dlistmask |= BEZP_DLM_PATCH;
  }
} /*GeomObjectDrawBezierPatch*/

void GeomObjectDrawBezierPNet ( GO_BezierPatch *obj )
{
  int ncp;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  if ( obj->me.dlistmask & BEZP_DLM_CNET )
    glCallList ( obj->me.displaylist+1 );
  else {
    glNewList ( obj->me.displaylist+1, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( OBJC_BEZP_CNET );
    DrawARectNet ( obj->me.cpdimen, obj->me.spdimen,
                   obj->degree_u+1, obj->degree_v+1,
                   obj->me.cpdimen*(obj->degree_v+1), obj->cpoints );
    ncp = (obj->degree_u+1)*(obj->degree_v+1);
    glColor3fv ( OBJC_BEZP_CP_UNMARKED );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, 0, MASK_CP_MOVEABLE );
    glColor3fv ( OBJC_BEZP_CP_MARKED );
    DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                  obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
    glEndList ();
    obj->me.dlistmask |= BEZP_DLM_CNET;
  }
} /*GeomObjectDrawBezierPNet*/

void GeomObjectDisplayBezierPatch ( GO_BezierPatch *obj )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  if ( obj->view_surf )
    GeomObjectDrawBezierPatch ( obj );
  if ( obj->view_cnet )
    GeomObjectDrawBezierPNet ( obj );
} /*GeomObjectDisplayBezierPatch*/

boolean GeomObjectBezierPatchSetDensityU ( GO_BezierPatch *obj, int densu )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( densu > 0 && densu <= MAX_PNET_DENSITY ) {
    obj->dens_u = densu;
    obj->me.dlistmask &= ~BEZP_DLM_PATCH;
    return true;
  }
  else
    return false;
} /*GeomObjectBezierPatchSetDensityU*/

boolean GeomObjectBezierPatchSetDensityV ( GO_BezierPatch *obj, int densv )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( densv > 0 && densv <= MAX_PNET_DENSITY ) {
    obj->dens_v = densv;
    obj->me.dlistmask &= ~BEZP_DLM_PATCH;
    return true;
  }
  else
    return false;
} /*GeomObjectBezierPatchSetDensityV*/

boolean GeomObjectBezierPatchSetDegreeU ( GO_BezierPatch *obj, int degu )
{
  int    ncp, degv;
  double *cp;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( degu < 1 || degu > MAX_DEGREE )
    return false;
  ncp = (degu+1)*(obj->degree_v+1);
  cp = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !cp || !mkcp )
    goto failure;
  if ( degu > obj->degree_u ) {
    if ( !mbs_BCDegElevPd ( obj->me.cpdimen, obj->degree_u, obj->degree_v,
                            obj->cpoints, degu-obj->degree_u, 0,
                            &degu, &degv, cp  ) )
      goto failure;
  }
  else if ( degu < obj->degree_u ) {
    if ( !mbs_BCDegRedPd ( obj->me.cpdimen, obj->degree_u, obj->degree_v,
                           obj->cpoints, obj->degree_u-degu, 0,
                           &degu, &degv, cp  ) )
      goto failure;
  }
  else
    memcpy ( cp, obj->cpoints, ncp*obj->me.cpdimen*sizeof(double) );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  obj->degree_u = degu;
  obj->degree_v = degv;
  obj->me.dlistmask = 0;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  GeomObjectBezierPatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBezierPatchSetDegreeU*/

boolean GeomObjectBezierPatchSetDegreeV ( GO_BezierPatch *obj, int degv )
{
  int    ncp, degu;
  double *cp;
  byte   *mkcp;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( degv < 1 || degv > MAX_DEGREE )
    return false;
  ncp = (obj->degree_u+1)*(degv+1);
  cp = malloc ( ncp*obj->me.cpdimen*sizeof(double) );
  mkcp = malloc ( ncp );
  if ( !cp || !mkcp )
    goto failure;
  if ( degv > obj->degree_v ) {
    if ( !mbs_BCDegElevPd ( obj->me.cpdimen, obj->degree_u, obj->degree_v,
                            obj->cpoints, 0, degv-obj->degree_v,
                            &degu, &degv, cp ) )
      goto failure;
  }
  else if ( degv < obj->degree_v ) {
    if ( !mbs_BCDegRedPd ( obj->me.cpdimen, obj->degree_u, obj->degree_v,
                           obj->cpoints, 0, obj->degree_v-degv,
                           &degu, &degv, cp ) )
      goto failure;
  }
  else
    memcpy ( cp, obj->cpoints, ncp*obj->me.cpdimen*sizeof(double) );
  free ( obj->cpoints );
  free ( obj->mkcp );
  obj->cpoints = cp;
  obj->mkcp = mkcp;
  memset ( mkcp, MASK_CP_MOVEABLE, ncp );
  obj->degree_u = degu;
  obj->degree_v = degv;
  obj->me.dlistmask = 0;
  if ( obj->savedcpoints ) {
    free ( obj->savedcpoints );
    obj->savedcpoints = NULL;
  }
  GeomObjectBezierPatchDisplayInfoText ( obj );
  return true;

failure:
  if ( cp ) free ( cp );
  if ( mkcp ) free ( mkcp );
  return false;
} /*GeomObjectBezierPatchSetDegreeV*/

boolean GeomObjectBezierPatchFlipUV ( GO_BezierPatch *obj )
{
  void    *sp, *mcp;
  int     ncp, n1, n2, pitch1, pitch2, els, i;

  sp = pkv_GetScratchMemTop ();

  n1 = obj->degree_v+1;
  n2 = obj->degree_u+1;
  ncp = n1*n2;
  els = obj->me.cpdimen*sizeof(double);
  mcp = pkv_GetScratchMem ( ncp*els );
  if ( !mcp )
    goto failure;
  pkv_TransposeMatrixc ( n2, n1, 1, n1, (char*)obj->mkcp, n2, mcp );
  memcpy ( obj->mkcp, mcp, ncp );
  pitch1 = n1*els;
  pitch2 = n2*els;
  pkv_TransposeMatrixc ( n2, n1, els, pitch1, (char*)obj->cpoints, pitch2, mcp );
  memcpy ( obj->cpoints, mcp, ncp*els );

  i = obj->degree_u;    obj->degree_u = obj->degree_v;      obj->degree_v = i;
  i = obj->dens_u;      obj->dens_u = obj->dens_v;          obj->dens_v = i;

  pkv_SetScratchMemTop ( sp );
  GeomObjectBezierPatchDisplayInfoText ( obj );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GeomObjectBezierPatchFlipUV*/

void GeomObjectBezierPatchFindBBox ( GO_BezierPatch *obj, Box3d *box )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  GeomObjectFindBBox ( obj->me.cpdimen, obj->rational,
                       (obj->degree_u+1)*(obj->degree_v+1), obj->cpoints, box );
} /*GeomObjectBezierPatchFindBBox*/

boolean GeomObjectBezierPatchFindCPoint ( GO_BezierPatch *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist )
{
  boolean ok;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  ok = GeomObjectFindNearestPoint ( obj->me.cpdimen, obj->me.spdimen,
             (obj->degree_u+1)*(obj->degree_v+1), obj->cpoints,
             obj->mkcp, MASK_CP_MOVEABLE, CPos, x, y, dist );
  return ok;
} /*GeomObjectBezierPatchFindCPoint*/

void GeomObjectBezierPatchSetCPoint ( GO_BezierPatch *obj,
                                      CameraRecd *CPos, short x, short y )
{
  int ncp, k;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  ncp = (obj->degree_u+1)*(obj->degree_v+1);
  if ( current_point_ind >= 0 && current_point_ind < ncp ) {
    k = obj->me.cpdimen*current_point_ind;
    GeomObjectSetPoint ( CPos, x, y, obj->me.cpdimen, obj->me.spdimen,
                         &obj->cpoints[k] );
  }
  obj->me.dlistmask = 0;
} /*GeomObjectBezierPatchSetCPoint*/

void GeomObjectBezierPatchMarkCPoints ( GO_BezierPatch *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        char mask, int action )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  GeomObjectMarkPoints ( obj->me.cpdimen, obj->me.spdimen,
                         (obj->degree_u+1)*(obj->degree_v+1),
                         obj->mkcp, obj->cpoints, CPos, box, mask, action );
  obj->me.dlistmask &= ~BEZP_DLM_CNET;
} /*GeomObjectBezierPatchMarkCPoints*/

void GeomObjectBezierPatchMarkCPoint ( GO_BezierPatch *obj,
                                       char mask, int action )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  GeomObjectMarkPoint ( (obj->degree_u+1)*(obj->degree_v+1),
                        obj->mkcp, mask, action );
  obj->me.dlistmask &= ~BEZP_DLM_CNET;
} /*GeomObjectBezierPatchMarkCPoint*/

boolean GeomObjectBezierPatchSaveCPoints ( GO_BezierPatch *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( obj->savedcpoints ) free ( obj->savedcpoints );
  size = (obj->degree_u+1)*(obj->degree_v+1)*obj->me.cpdimen*sizeof(double);
  obj->savedcpoints = malloc ( size );
  if ( !obj->savedcpoints )
    return false;
  memcpy ( obj->savedcpoints, obj->cpoints, size );
  obj->savedsize = size;
  return true;
} /*GeomObjectBezierPatchSaveCPoints*/

void GeomObjectBezierPatchUndoLastTransformation ( GO_BezierPatch *obj )
{
  int size;

  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  size = (obj->degree_u+1)*(obj->degree_v+1)*obj->me.cpdimen*sizeof(double);
  if ( obj->savedcpoints && obj->savedsize == size ) {
    memcpy ( obj->cpoints, obj->savedcpoints, size );
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBezierPatchUndoLastTransformation*/

void GeomObjectBezierPatchTransformCPoints ( GO_BezierPatch *obj,
                                             trans3d *tr, char mask )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  if ( obj->savedcpoints ) {
    GeomObjectTransformPoints ( obj->me.cpdimen, obj->me.spdimen,
        (obj->degree_u+1)*(obj->degree_v+1), obj->mkcp, mask,
        obj->savedcpoints, obj->cpoints, tr );
    obj->me.dlistmask = 0;
  }
} /*GeomObjectBezierPatchTransformCPoints*/

boolean GeomObjectBezierPatchGetPointCoord ( GO_BezierPatch *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return false;
  if ( p < 0 || p >= (obj->degree_u+1)*(obj->degree_v+1) )
    return false;
  *spdimen = obj->me.spdimen;
  *cpdimen = obj->me.cpdimen;
  *pc = &obj->cpoints[obj->me.cpdimen*p];
  return true;
} /*GeomObjectBezierPatchGetPointCoord*/

boolean GeomObjectWriteBPAttributes ( GO_BezierPatch *obj )
{
  return true;
} /*GeomObjectWriteBPAttributes*/

boolean GeomObjectWriteBezierPatch ( GO_BezierPatch *obj )
{
  return bsf_WriteBezierPatchd ( obj->me.spdimen, obj->me.cpdimen,
                                 obj->rational, obj->degree_u, obj->degree_v,
                                 obj->me.cpdimen*(obj->degree_v+1),
                                 obj->cpoints, obj->me.name, obj->me.ident,
                                 GeomObjectWriteAttributes,
                                 (void*)&obj->me );
} /*GeomObjectWriteBezierPatch*/

boolean GeomObjectBPResolveDependencies ( GO_BezierPatch *obj )
{
  return true;
} /*GeomObjectBPResolveDependencies*/

void GeomObjectReadBezierPatch ( void *usrdata, const char *name, int ident,
                                 int degreeu, int degreev,
                                 int pitch, const point4d *cpoints, int spdimen,
                                 boolean rational )
{
  GO_BezierPatch       *obj;
  double               *cp;
  byte                 *mkcp;
  int                  ncp, cpdimen, i;
  rw_object_attributes *attrib;

  attrib = (rw_object_attributes*)usrdata;
  obj = (GO_BezierPatch*)attrib->go_being_read;
  if ( obj ) {
    strncpy ( obj->me.name, name, MAX_NAME_LENGTH+1 );
    obj->rational = rational;
    obj->me.spdimen = spdimen;
    obj->me.cpdimen = rational ? spdimen+1 : spdimen;
    obj->me.ident = ident;
    ncp = (degreeu+1)*(degreev+1);
    cpdimen = obj->me.cpdimen;
    cp = malloc ( ncp*cpdimen*sizeof(double) );
    mkcp = malloc ( ncp );
    if ( !cp || !mkcp ) {
      if ( mkcp ) free ( mkcp );
      if ( cp ) free ( cp );
      GeomObjectDeleteBezierPatch ( obj );
      return;
    }
    free ( obj->cpoints );
    free ( obj->mkcp );
    if ( obj->savedcpoints ) {
      free ( obj->savedcpoints );
      obj->savedcpoints = NULL;
    }
    obj->cpoints = cp;
    obj->mkcp = mkcp;
    obj->degree_u = degreeu;
    obj->degree_v = degreev;
    for ( i = 0; i <= degreeu; i++ )
      GeomObjectSetupIniPoints ( spdimen, rational, &obj->me.cpdimen,
                                 degreev+1, &((double*)cpoints)[i*pitch],
                                 &cp[i*cpdimen*(degreev+1)] );
  }
} /*GeomObjectReadBezierPatch*/

void GeomObjectBezierPatchOutputToRenderer3D ( GO_BezierPatch *obj )
{
  if ( obj->me.obj_type != GO_BEZIER_PATCH )
    return;
  if ( obj->me.spdimen != 3 )
    return;
  if ( obj->rational )
    RendEnterBezPatch3Rd ( &rend, obj->degree_u, obj->degree_v,
                           (point4d*)obj->cpoints, obj->me.colour );
  else
    RendEnterBezPatch3d ( &rend, obj->degree_u, obj->degree_v,
                          (point3d*)obj->cpoints, obj->me.colour );
} /*GeomObjectBezierPatchOutputToRenderer3D*/

void GeomObjectBezierPatchDisplayInfoText ( GO_BezierPatch *obj )
{
  char s[80];

  sprintf ( s, "deg = (%d,%d)", obj->degree_u, obj->degree_v );
  SetStatusText ( s, true );
} /*GeomObjectBezierPatchDisplayInfoText*/

boolean GeomObjectBezierPatchProcessDep ( GO_BezierPatch *obj, geom_object *go )
{
  return false;
} /*GeomObjectBezierPatchProcessDep*/

void GeomObjectBezierPatchProcessDeletedDep ( GO_BezierPatch *obj, geom_object *go )
{
} /*GeomObjectBezierPatchProcessDeletedDep*/

