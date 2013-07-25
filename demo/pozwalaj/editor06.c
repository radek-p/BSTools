
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

/*#include "widgets.h"*/
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"

boolean GeomObjectSelectPoint ( char spdimen, CameraRecd *CPos, short x, short y )
{
  if ( current_go->obj_type == GO_BSPLINE_MESH &&
       spdimen == current_go->spdimen )
    return GeomObjectBSplineMeshSelectPoint ( (GO_BSplineMesh*)current_go,
                                              CPos, x, y );
  else
    return false;
} /*GeomObjectSelectPoint*/

boolean GeomObjectSelectEdge ( char spdimen, CameraRecd *CPos, short x, short y )
{
  if ( current_go->obj_type == GO_BSPLINE_MESH &&
       spdimen == current_go->spdimen )
    return GeomObjectBSplineMeshSelectEdge ( (GO_BSplineMesh*)current_go,
                                             CPos, x, y );
  else
    return false;
} /*GeomObjectSelectEdge*/

void GeomObjectUnselectPoint ( char spdimen )
{
  if ( current_go->obj_type == GO_BSPLINE_MESH &&
       spdimen == current_go->spdimen )
    GeomObjectBSplineMeshUnselectPoint ( (GO_BSplineMesh*)current_go );
} /*GeomObjectUnselectPoint*/

void GeomObjectUnselectEdge ( char spdimen )
{
  if ( current_go->obj_type == GO_BSPLINE_MESH &&
       spdimen == current_go->spdimen )
    GeomObjectBSplineMeshUnselectEdge ( (GO_BSplineMesh*)current_go );
} /*GeomObjectUnselectEdge*/

