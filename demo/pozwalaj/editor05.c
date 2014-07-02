
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

/*#include "widgets.h"*/
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

void GeomObjectOutputToRenderer3D ( boolean all )
{
  geom_object *go;

  for ( go = first_go; go; go = go->next )
    if ( go->spdimen == 3 && ( all || go->active || go == current_go) ) {
      switch ( go->obj_type ) {
    case GO_BEZIER_CURVE:
        GeomObjectBezierCurveOutputToRenderer ( (GO_BezierCurve*)go );
        break;
    case GO_BSPLINE_CURVE:
        GeomObjectBSplineCurveOutputToRenderer ( (GO_BSplineCurve*)go );
        break;
    case GO_BEZIER_PATCH:
        GeomObjectBezierPatchOutputToRenderer3D ( (GO_BezierPatch*)go );
        break;
    case GO_BSPLINE_PATCH:
        GeomObjectBSplinePatchOutputToRenderer3D ( (GO_BSplinePatch*)go );
        break;
    case GO_BSPLINE_MESH:
        GeomObjectBSplineMeshOutputToRenderer3D ( (GO_BSplineMesh*)go );
        break;
    default:
        break;
      }
    }
} /*GeomObjectOutputToRenderer3D*/

