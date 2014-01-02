
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
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

#include "widgets.h"
#include "editor.h"  
#include "editor_bezc.h"
#include "editor_bsc.h"
#include "editor_bezp.h"
#include "editor_bsp.h"
#include "editor_bsm.h"
#include "editor_bsh.h"
#include "pozwalaj.h"


void SetupObjectSpecificMenus ( geom_object *obj )
{
  if ( !obj ) {
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    SetStatusText ( "", false );
    SetGeomWin10Empty ();
  }
  else {
        /* show the appropriate editing window in win0: 2D or 3D */
    if ( (side00menu->data1 != side00ewidgets &&
          side00menu->data1 != side00bwidgets &&
          side00menu->data1 != side00cwidgets) ||
          obj->spdimen != 3 ) {
      if ( obj->obj_type == GO_BSPLINE_MESH )
        side00widgets = side00dwidgets;
      else
        side00widgets = side00awidgets;
      xge_SetMenuWidgets ( side00menu, side00widgets, false );
    }
    switch ( obj->spdimen ) {
  case 2:
      SetWin002D ();
      break;
  case 3:
      SetWin003D ();
      break;
  default:
      return;
    }
        /* choose the side menu */
    xge_SetWindow ( win1 );
    GeomObjectDisplayInfoText ( obj );
    switch ( whichside10menu ) {
  case SIDE10MENU_EDIT:
      SetObjectEditMenu ( obj );
      return;
  case SIDE10MENU_VIEW:
      SetObjectViewMenu ( obj );
      return;
  case SIDE10MENU_DATA:
      SetObjectDataMenu ( obj );
      return;
  case SIDE10MENU_OPTIONS:
      SetObjectOptionsMenu ( obj );
      return;
  default:
      xge_SetMenuWidgets ( side10menu, NULL, false );
      return;
    }
  }
} /*SetupObjectSpecificMenus*/

void SetObjectEditMenu ( geom_object *obj )
{
  if ( !obj ) {
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    return;
  }
        /* now show the menu specific for the object */
  whichside10menu = SIDE10MENU_EDIT;
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezc, false );
    SetupBezierCurveWidgets ( (GO_BezierCurve*)obj );
    break;
case GO_BSPLINE_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsc, false );
    SetupBSplineCurveWidgets ( (GO_BSplineCurve*)obj );
    break;
case GO_BEZIER_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezp, false );
    SetupBezierPatchWidgets ( (GO_BezierPatch*)obj );
    break;
case GO_BSPLINE_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsp, false );
    SetupBSplinePatchWidgets ( (GO_BSplinePatch*)obj );
    break;
case GO_BSPLINE_MESH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsm, false );
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    Side10MenuBsmCallBack ( side10menu, xgemsg_RESIZE, 0,
                            side10menu->w, side10menu->h );
    break;
case GO_BSPLINE_HOLE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsh, false );
    SetupBSplineHoleWidgets ( (GO_BSplineHole*)obj );
    break;
default:
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    break;
  }
} /*SetObjectEditMenu*/

void SetObjectViewMenu ( geom_object *obj )
{
  if ( !obj ) {
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    return;
  }
  whichside10menu = SIDE10MENU_VIEW;
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezc_view, false );
    SetupBezierCurveWidgets ( (GO_BezierCurve*)obj );
    break;
case GO_BSPLINE_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsc_view, false );
    SetupBSplineCurveWidgets ( (GO_BSplineCurve*)obj );
    break;
case GO_BEZIER_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezp_view, false );
    SetupBezierPatchWidgets ( (GO_BezierPatch*)obj );
    break;
case GO_BSPLINE_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsp_view, false );
    SetupBSplinePatchWidgets ( (GO_BSplinePatch*)obj );
    break;
case GO_BSPLINE_MESH:
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsm_view, false );
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    break;
case GO_BSPLINE_HOLE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsh_view, false );
    SetupBSplineHoleWidgets ( (GO_BSplineHole*)obj );
    break;
default:
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    break;
  }
} /*SetObjectViewMenu*/

void SetObjectDataMenu ( geom_object *obj )
{
  if ( !obj ) {
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    return;
  }
  whichside10menu = SIDE10MENU_DATA;
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezc_data, false );
    SetupBezierCurveWidgets ( (GO_BezierCurve*)obj );
    break;
case GO_BSPLINE_CURVE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsc_data, false );
    SetupBSplineCurveWidgets ( (GO_BSplineCurve*)obj );
    break;
case GO_BEZIER_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bezp_data, false );
    SetupBezierPatchWidgets ( (GO_BezierPatch*)obj );
    break;
case GO_BSPLINE_PATCH:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsp_data, false );
    SetupBSplinePatchWidgets ( (GO_BSplinePatch*)obj );
    break;
case GO_BSPLINE_MESH:
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsm_data, false );
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    break;
case GO_BSPLINE_HOLE:
    xge_SetMenuWidgets ( side10menu, side10wdg_bsh_data, false );
    SetupBSplineHoleWidgets ( (GO_BSplineHole*)obj );
    break;
default:
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    break;
  }
} /*SetObjectDataMenu*/

void SetObjectOptionsMenu ( geom_object *obj )
{
  if ( !obj ) {
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    return;
  }
        /* now show the menu specific for the object */
  whichside10menu = SIDE10MENU_OPTIONS;
  switch ( obj->obj_type ) {
case GO_BEZIER_CURVE:
    SetupBezierCurveWidgets ( (GO_BezierCurve*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bezc_opt, false );
    break;
case GO_BSPLINE_CURVE:
    SetupBSplineCurveWidgets ( (GO_BSplineCurve*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsc_opt, false );
    break;
case GO_BEZIER_PATCH:
    SetupBezierPatchWidgets ( (GO_BezierPatch*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bezp_opt, false );
    break;
case GO_BSPLINE_PATCH:
    SetupBSplinePatchWidgets ( (GO_BSplinePatch*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsp_opt, false );
    break;
case GO_BSPLINE_MESH:
    SetupBSplineMeshWidgets ( (GO_BSplineMesh*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsm_opt, false );
    Side10MenuBsmCallBack ( side10menu, xgemsg_RESIZE, 0,
                            side10menu->w, side10menu->h );
    break;
case GO_BSPLINE_HOLE:
    SetupBSplineHoleWidgets ( (GO_BSplineHole*)obj );
    xge_SetMenuWidgets ( side10menu, side10wdg_bsh_opt, false );
    break;
default:
    xge_SetMenuWidgets ( side10menu, side10wdg_none, false );
    break;
  }
} /*SetObjectOptionsMenu*/

