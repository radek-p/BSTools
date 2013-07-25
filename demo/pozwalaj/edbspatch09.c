
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


boolean GeomObjectBSplinePatchProcessDep ( GO_BSplinePatch *obj, geom_object *go )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return false;
  switch ( obj->bsp_type ) {
case BSP_TYPE_GENERAL:
    break;
case BSP_TYPE_SWEPT:
    break;
case BSP_TYPE_SPHERICAL:
    return GeomObjectBSplinePatchGenSphericalProduct ( obj );
case BSP_TYPE_LOFTED:
    break;
case BSP_TYPE_BLENDING_G1:
    break;
case BSP_TYPE_BLENDING_G2:
    break;
default:
    break;
  }
  return false;
} /*GeomObjectBSplinePatchProcessDep*/

void GeomObjectBSplinePatchProcessDeletedDep ( GO_BSplinePatch *obj, geom_object *go )
{
  if ( obj->me.obj_type != GO_BSPLINE_PATCH )
    return;
  switch ( obj->bsp_type ) {
case BSP_TYPE_GENERAL:
    break;
case BSP_TYPE_SWEPT:
    break;
case BSP_TYPE_SPHERICAL:
        /* either the equator or the meridian is deleted */
    GeomObjectBSplinePatchAdjustGeneral ( obj );
    break;
case BSP_TYPE_LOFTED:
    break;
case BSP_TYPE_BLENDING_G1:
    break;
case BSP_TYPE_BLENDING_G2:
    break;
default:
    break;
  }
} /*GeomObjectBSplinePatchProcessDeletedDep*/

