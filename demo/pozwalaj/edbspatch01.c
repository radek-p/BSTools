
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
#include "edcolours.h"
#include "editor_bsp.h"


void GeomObjectDrawBSplinePatch ( GO_BSplinePatch *obj )
{
  int               pitch;
  BSP_TrimmedDomain *trpd;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( obj->me.dlistmask & BSP_DLM_PATCH )
    glCallList ( obj->me.displaylist );
  else {
    glNewList ( obj->me.displaylist, GL_COMPILE_AND_EXECUTE );
    pitch = obj->me.cpdimen*(obj->lastknot_v-obj->degree_v);
    if ( (trpd = obj->trpd) ) { /* a trimmed patch */
      glColor3fv ( OBJC_BSP_TRIMMED_BOUNDARY );
      DrawTrimmedBSplinePatchBoundary ( obj->me.cpdimen, obj->me.spdimen,
                                  obj->degree_u, obj->lastknot_u, obj->knots_u,
                                  obj->degree_v, obj->lastknot_v, obj->knots_v,
                                  pitch, obj->cpoints, obj->dens_u, obj->dens_v,
                                  true, !obj->closed_u, true, !obj->closed_v,
                                  trpd->nelem, trpd->bound );
      glColor3fv ( OBJC_BSP_PATCH_CPLINES );
      DrawTrimmedBSplinePatchWF ( obj->me.cpdimen, obj->me.spdimen,
                                  obj->degree_u, obj->lastknot_u, obj->knots_u,
                                  obj->degree_v, obj->lastknot_v, obj->knots_v,
                                  pitch, obj->cpoints, obj->dens_u, obj->dens_v,
                                  true, !obj->closed_u, true, !obj->closed_v,
                                  trpd->nelem, trpd->bound );
    }
    else {  /* a non-trimmed patch */
      glColor3fv ( OBJC_BSP_PATCH_CPLINES );
      DrawBSplinePatchWF ( obj->me.cpdimen, obj->me.spdimen,
                           obj->degree_u, obj->lastknot_u, obj->knots_u,
                           obj->degree_v, obj->lastknot_v, obj->knots_v,
                           pitch, obj->cpoints, obj->dens_u, obj->dens_v,
                           true, !obj->closed_u, true, !obj->closed_v );
    }
    glEndList ();
    obj->me.dlistmask |= BSP_DLM_PATCH;
  }
} /*GeomObjectDrawBSplinePatch*/

void GeomObjectDrawBSplinePNet ( GO_BSplinePatch *obj )
{
  int ncp;

  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( obj->me.dlistmask & BSP_DLM_CNET )
    glCallList ( obj->me.displaylist+1 );
  else {
    glNewList ( obj->me.displaylist+1, GL_COMPILE_AND_EXECUTE );
    glColor3fv ( OBJC_BSP_CNET );
    DrawARectNet ( obj->me.cpdimen, obj->me.spdimen,
                   obj->lastknot_u-obj->degree_u,
                   obj->lastknot_v-obj->degree_v,
                   obj->me.cpdimen*(obj->lastknot_v-obj->degree_v),
                   obj->cpoints );
    ncp = (obj->lastknot_u-obj->degree_u)*(obj->lastknot_v-obj->degree_v);
    switch ( obj->bsp_type ) {
  case BSP_TYPE_BLENDING_G1:
  case BSP_TYPE_BLENDING_G2:
      glColor3fv ( OBJC_BSP_CP_UNMARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                    obj->mkcp, 0, MASK_CP_MOVEABLE | MASK_CP_BOUNDARY );
      glColor3fv ( OBJC_BSP_CP_MARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                    obj->mkcp, marking_mask,
                    MASK_MARKED | MASK_CP_MOVEABLE | MASK_CP_BOUNDARY );
      break;
  default:
      glColor3fv ( OBJC_BSP_CP_UNMARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                    obj->mkcp, 0, MASK_CP_MOVEABLE );
      glColor3fv ( OBJC_BSP_CP_MARKED );
      DrawCPoints ( obj->me.cpdimen, obj->me.spdimen, ncp, obj->cpoints,
                    obj->mkcp, marking_mask, MASK_MARKED | MASK_CP_MOVEABLE );
      break;
    }
    glEndList ();
    obj->me.dlistmask |= BSP_DLM_CNET;
  }
} /*GeomObjectDrawBSplinePNet*/

void GeomObjectDisplayBSplinePatch ( GO_BSplinePatch *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  if ( obj->view_surf )
    GeomObjectDrawBSplinePatch ( obj );
  if ( obj->view_cnet )
    GeomObjectDrawBSplinePNet ( obj );
} /*GeomObjectDisplayBSplinePatch*/

boolean GeomObjectBSplinePatchSetDensityU ( GO_BSplinePatch *obj, int densu )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( densu > 0 && densu <= MAX_PNET_DENSITY ) {
    obj->dens_u = densu;
    obj->me.dlistmask &= ~BSP_DLM_PATCH;
    return true;
  }
  else
    return false;
} /*GeomObjectBSplinePatchSetDensityU*/

boolean GeomObjectBSplinePatchSetDensityV ( GO_BSplinePatch *obj, int densv )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  if ( densv > 0 && densv <= MAX_PNET_DENSITY ) {
    obj->dens_v = densv;
    obj->me.dlistmask &= ~BSP_DLM_PATCH;
    return true;
  }
  else
    return false;
} /*GeomObjectBSplinePatchSetDensityV*/

