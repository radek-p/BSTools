
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
#include "render.h"

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE


static boolean BlBSplineEntire ( int deg, int lknu, int lknv, int range[4] )
{
  if ( range[0] == deg && range[1] == lknu-2*deg-1 &&
       range[2] == deg && range[3] == lknv-2*deg-1 )
    return true;
  else
    return false;
} /*BlBSplineEntire*/

void SetupBSplinePatchWidgets ( GO_BSplinePatch *obj )
{
  GO_BSplineCurve *equator, *meridian;

  InitNameEditor ( &bsp_name_ed, obj->me.name );
  degreeu = obj->degree_u;
  degreev = obj->degree_v;
  intwdensityu.title = &txtDensityU[0];
  intwdensityv.title = &txtDensityV[0];
  density_u = obj->dens_u;
  density_v = obj->dens_v;
  sw_view_surf = obj->view_surf;
  sw_view_cnet = obj->view_cnet;
  bsp_sw_closed_u = obj->closed_u;
  bsp_sw_closed_v = obj->closed_v;
  Geom10winT2SetKnots ( &g10t2knotwin,
                        obj->degree_u, obj->lastknot_u, obj->knots_u,
                        obj->closed_u,
                        obj->degree_v, obj->lastknot_v, obj->knots_v,
                        obj->closed_v );
  switch ( obj->bsp_type ) {
case BSP_TYPE_GENERAL:
general:
    side10wdg_bsp_opt = side10wdg_bsp_opt_general;
    g10t2knotwin.locked_u = g10t2knotwin.locked_v = false;
    bsp_type_button->data0 = txtGeneral;
    SetGeomWin10T2Knotw ( obj );
    break;

case BSP_TYPE_SPHERICAL:
    if ( obj->me.maxdn < 2 || !obj->me.dependencies )
      goto failure_spr;
    equator = (GO_BSplineCurve*)obj->me.dependencies[0];
    meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
    if ( !equator || !meridian )
      goto failure_spr;
    if ( equator->me.obj_type != GO_BSPLINE_CURVE ||
         meridian->me.obj_type != GO_BSPLINE_CURVE )
      goto failure_spr;
    sw_bsp_sproduct_meridian = !sw_bsp_sproduct_equator;
    sw_bsp_sproduct_eqrational = equator->rational;
    sw_bsp_sproduct_eqclosed = equator->closed;
    sw_bsp_sproduct_equniform = equator->uniform;
    bsp_eqdegree = equator->degree;
    InitNameEditor ( &bsp_sproduct_eqname_ed, obj->eqname );
    sw_bsp_sproduct_merrational = meridian->rational;
    sw_bsp_sproduct_merclosed = meridian->closed;
    sw_bsp_sproduct_meruniform = meridian->uniform;
    bsp_merdegree = meridian->degree;
    InitNameEditor ( &bsp_sproduct_mername_ed, obj->mername );
    side10wdg_bsp_opt = side10wdg_bsp_opt_spherical;
    if ( sw_bsp_sproduct_equator ) {
      Geom10winKNSetKnots ( &g10knotwineqmer, equator->degree,
                            equator->lastknot, equator->knots, equator->closed );
      SetGeomWin10SPrEqMer ( equator );
    }
    else {
      Geom10winKNSetKnots ( &g10knotwineqmer, meridian->degree,
                            meridian->lastknot, meridian->knots, meridian->closed );
      SetGeomWin10SPrEqMer ( meridian );
    }
    bsp_type_button->data0 = txtSpherical;
    break;
failure_spr:
    GeomObjectBSplinePatchAdjustGeneral ( obj );
    goto general;

case BSP_TYPE_BLENDING_G1:
    side10wdg_bsp_opt = side10wdg_bsp_opt_blendingG1;
    memcpy ( bsp_bl_opt_range, obj->blp_range, 4*sizeof(int) );
    bsp_blending_entire = BlBSplineEntire ( 2, obj->lastknot_u, obj->lastknot_v,
                                            bsp_bl_opt_range );
    bsp_type_button->data0 = txtBlendingG1;
    goto ifchild;

case BSP_TYPE_BLENDING_G2:
    side10wdg_bsp_opt = side10wdg_bsp_opt_blendingG2;
    memcpy ( bsp_bl_opt_range, obj->blp_range, 4*sizeof(int) );
    bsp_blending_entire = BlBSplineEntire ( 3, obj->lastknot_u, obj->lastknot_v,
                                            bsp_bl_opt_range );
    bsp_type_button->data0 = txtBlendingG2;
ifchild:
    if ( obj->me.bound_with_a_child )
      bsp_bl_optimizeG1->data0 = bsp_bl_optimizeG2->data0 = txtInterrupt;
    else
      bsp_bl_optimizeG1->data0 = bsp_bl_optimizeG2->data0 = txtOptimize;
    sw_bsp_clamped = obj->clamped;
    g10t2knotwin.locked_u = g10t2knotwin.locked_v = true;
    SetGeomWin10T2Knotw ( obj );
    break;

default:  /* to be implemented */
    g10t2knotwin.locked_u = g10t2knotwin.locked_v = false;
    SetGeomWin10T2Knotw ( obj );
    break;
  }
  bsp_bl_nkn1 = obj->nkn1;
  bsp_bl_nkn2 = obj->nkn2;
  bsp_bl_maxit = obj->maxit;
  sl_bsp_bl_param = log ( obj->blp_C/BSP_BL_MIN_CPARAM )/
                    log ( BSP_BL_MAX_CPARAM/BSP_BL_MIN_CPARAM );
} /*SetupBSplinePatchWidgets*/

void BSPatchSetupGeneralOptions ( GO_BSplinePatch *obj )
{
  if ( obj->bsp_type != BSP_TYPE_GENERAL ) {
    GeomObjectBSplinePatchAdjustGeneral ( obj );
    GeomObjectBSPatchMarkCPGeneral ( obj );
  }
  SetupBSplinePatchWidgets ( obj );
  xge_RedrawAll ();
} /*BSPatchSetupGeneralOptions*/

void BSPatchSetupSphericalOptions ( GO_BSplinePatch *obj )
{
  if ( obj->bsp_type != BSP_TYPE_SPHERICAL ) {
    GeomObjectBSplinePatchAdjustSProduct ( obj );
    GeomObjectBSPatchMarkCPGeneral ( obj ); /* ?? */
  }
  SetupBSplinePatchWidgets ( obj );
  xge_RedrawAll ();
} /*BSPatchSetupSphericalOptions*/

void BSPatchSetupBlG1Options ( GO_BSplinePatch *obj )
{
  if ( obj->bsp_type != BSP_TYPE_BLENDING_G1 ) {
    GeomObjectBSplinePatchAdjustBlG1 ( obj );
    xge_T2KnotWindFindMapping ( &g10t2knotwin );
  }
  SetupBSplinePatchWidgets ( obj );
  xge_RedrawAll ();
} /*BSPatchSetupBlG1Options*/

void BSPatchSetupBlG2Options ( GO_BSplinePatch *obj )
{
  if ( obj->bsp_type != BSP_TYPE_BLENDING_G2 ) {
    GeomObjectBSplinePatchAdjustBlG2 ( obj );
    xge_T2KnotWindFindMapping ( &g10t2knotwin );
  }
  SetupBSplinePatchWidgets ( obj );
  xge_RedrawAll ();
} /*BSPatchSetupBlG2Options*/

static boolean SetMinRange ( int minr, int maxr, int newv, int *minvar, int *maxvar )
{
  if ( newv < minr || newv > maxr ) return false;
  *minvar = newv;
  if ( *maxvar < newv ) *maxvar = newv;
  return true;
} /*SetMinRange*/

static boolean SetMaxRange ( int minr, int maxr, int newv, int *minvar, int *maxvar )
{
  if ( newv < minr || newv > maxr ) return false;
  *maxvar = newv;
  if ( *minvar > newv ) *minvar = newv;
  return true;
} /*SetMaxRange*/

int Side10MenuBspCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplinePatch *obj;
  GO_BSplineCurve *equator, *meridian, *saveeqmer;

  if ( current_go->obj_type != GO_BSPLINE_PATCH )
    return 0;
  obj = (GO_BSplinePatch*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BSP_TYPE:
      OpenPopup ( popup15, false );
      return 1;
  case btnM1BSP_UNIFORM_U:
      if ( GeomObjectBSplinePatchSetUniformKnotsU ( obj ) ) {
        xge_T2KnotWindFindMapping ( &g10t2knotwin );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSP_UNIFORM_V:
      if ( GeomObjectBSplinePatchSetUniformKnotsV ( obj ) ) {
        xge_T2KnotWindFindMapping ( &g10t2knotwin );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSP_FLIP:
      if ( GeomObjectBSplinePatchFlipUV ( obj ) ) {
        SetupBSplinePatchWidgets ( obj );
        xge_T2KnotWindFindMapping ( &g10t2knotwin );
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotFlip, 0 );
      return 1;
  case btnM1BSP_PRETRANSFORMATION:
      GetPreTransformation ();
      OpenPopup ( popup14, false );
      obj->me.display_pretrans = true;
      editing_pretrans = true;
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_SetWindow ( win0 );
      xge_Redraw ();
      return 1;
  case btnM1BSP_BLENDING_OPTIMIZE:
      if ( ipc_state == ipcstate_CHILD_BUSY ) {
        if ( obj->me.bound_with_a_child ) {
          IPCInterruptTheChild ();
          bsp_bl_optimizeG1->data0 = bsp_bl_optimizeG2->data0 = txtOptimize;
        }
        else {
          xge_DisplayErrorMessage ( ErrorMsgChildProcessBusy, 0 );
          return 1;
        }
      }
      else {
        if ( BlendingPatchOptimizationPrepareData ( obj ) ) {
          bsp_sw_closed_u = obj->closed_u;
          bsp_bl_nharmonic = false;
          GeomObjectBSplinePatchSetNHarmonic ( obj, bsp_bl_nharmonic );
          InitBlendingPatchOptimization ();
          bsp_bl_optimizeG1->data0 = bsp_bl_optimizeG2->data0 = txtInterrupt;
        }
      }
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case btnM1BSP_BLENDING_REFINE:
      if ( GeomObjectBSplinePatchRefine ( obj ) ) {
        SetupBSplinePatchWidgets ( obj );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSP_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM1BSP_VIEW_SURF:
      obj->view_surf = sw_view_surf;
      RedrawGeom00Win ();
      return 1;
  case swM1BSP_VIEW_CNET:
      obj->view_cnet = sw_view_cnet;
      RedrawGeom00Win ();
      return 1;
  case swM1BSP_CLOSED_U:
      if ( GeomObjectBSplinePatchSetClosedU ( obj, bsp_sw_closed_u ) ) {
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( bsp_sw_closed_u )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      bsp_sw_closed_u = obj->closed_u;
      return 1;
  case swM1BSP_CLOSED_V:
      if ( GeomObjectBSplinePatchSetClosedV ( obj, bsp_sw_closed_v ) ) {
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( bsp_sw_closed_v )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      bsp_sw_closed_v = obj->closed_v;
      return 1;
  case swM1BSP_DOMAIN_COORD:
      g10t2knotwin.display_coord = g10win2D.display_coord =
      g10win2Deqmer.display_coord = g10knotwin.display_coord =
      g10knotwineqmer.display_coord = sw_bsp_dom_coord;
      return 1;
  case swM1BSP_DOMAIN_PANZOOM:
      g10t2knotwin.panning = g10win2D.panning = g10win2Deqmer.panning =
      g10knotwineqmer.panning = sw_bsp_dom_panzoom;
      if ( sw_bsp_dom_panzoom )
        g10t2knotwin.selecting_mode = g10win2D.selecting_mode =
        g10win2Deqmer.selecting_mode = false;
      xge_Redraw ();
      return 1;
  case swM1BSP_BLENDING_CLAMPED:
      if ( sw_bsp_clamped )
        GeomObjectBSplinePatchFreeToClamped ( obj );
      else
        GeomObjectBSplinePatchClampedToFree ( obj );
      sw_bsp_clamped = obj->clamped;
      xge_RedrawAll ();
      return 1;
  case swM1BSP_NHARMONIC:
      if ( GeomObjectBSplinePatchSetNHarmonic ( obj, bsp_bl_nharmonic ) ) {
        SetupBSplinePatchWidgets ( obj );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case swM1BSP_BLENDING_ENTIRE:
      if ( bsp_blending_entire ) {
        GeomObjectBSplinePatchSetEntireBlendingRange ( obj );
        memcpy ( bsp_bl_opt_range, obj->blp_range, 4*sizeof(int) );
        xge_RedrawAll ();
      }
      return 1;
  case swM1BSP_SPRODUCT_EQUATOR:
      sw_bsp_sproduct_meridian = !sw_bsp_sproduct_equator;
      goto switch_eqmer;
  case swM1BSP_SPRODUCT_MERIDIAN:
      sw_bsp_sproduct_equator = !sw_bsp_sproduct_meridian;
switch_eqmer:
      equator = (GO_BSplineCurve*)obj->me.dependencies[0];
      meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
      if ( sw_bsp_sproduct_equator ) {
        Geom10winKNSetKnots ( &g10knotwineqmer, equator->degree,
                              equator->lastknot, equator->knots, equator->closed );
        SetGeomWin10SPrEqMer ( equator );
      }
      else {
        Geom10winKNSetKnots ( &g10knotwineqmer, meridian->degree,
                              meridian->lastknot, meridian->knots, meridian->closed );
        SetGeomWin10SPrEqMer ( meridian );
      }
      xge_Redraw ();
      return 1;
  case swM1BSP_SPRODUCT_EQRATIONAL:
      equator = (GO_BSplineCurve*)obj->me.dependencies[0];
      if ( sw_bsp_sproduct_eqrational )
        GeomObjectBSplineCurveSetRational ( equator );
      else
        GeomObjectBSplineCurveSetNonRational ( equator );
      sw_bsp_sproduct_eqrational = equator->rational;
      goto redraw_eq;
  case swM1BSP_SPRODUCT_EQCLOSED:
      equator = (GO_BSplineCurve*)obj->me.dependencies[0];
      if ( GeomObjectBSplineCurveSetClosed ( equator, sw_bsp_sproduct_eqclosed ) ) {
        bsp_sw_closed_u = sw_bsp_sproduct_eqclosed = equator->closed;
        goto redraw_eq;
      }
      else if ( sw_bsp_sproduct_eqclosed )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      return 1;
  case swM1BSP_SPRODUCT_EQUNIFORM:
      equator = (GO_BSplineCurve*)obj->me.dependencies[0];
      GeomObjectBSplineCurveSetUniformKnots ( equator,
                                      sw_bsp_sproduct_equniform );
      if ( sw_bsp_sproduct_equniform ) {
redraw_eq:
        if ( sw_bsp_sproduct_equator )
          xge_Redraw ();
        xge_SetWindow ( win0 );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_Redraw ();
      }
      return 1;
  case swM1BSP_SPRODUCT_MERRATIONAL:
      meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
      if ( sw_bsp_sproduct_merrational )
        GeomObjectBSplineCurveSetRational ( meridian );
      else
        GeomObjectBSplineCurveSetNonRational ( meridian );
      sw_bsp_sproduct_merrational = meridian->rational;
      goto redraw_mer;
  case swM1BSP_SPRODUCT_MERCLOSED:
      meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
      if ( GeomObjectBSplineCurveSetClosed ( meridian, sw_bsp_sproduct_merclosed ) ) {
        bsp_sw_closed_v = sw_bsp_sproduct_merclosed = meridian->closed;
        goto redraw_mer;
      }
      else if ( sw_bsp_sproduct_merclosed )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      bsp_sw_closed_v = sw_bsp_sproduct_merclosed = meridian->closed;
      return 1;
  case swM1BSP_SPRODUCT_MERUNIFORM:
      meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
      GeomObjectBSplineCurveSetUniformKnots ( meridian,
                                      sw_bsp_sproduct_meruniform );
      if ( sw_bsp_sproduct_meruniform ) {
redraw_mer:
        if ( sw_bsp_sproduct_meridian )
          xge_Redraw ();
        xge_SetWindow ( win0 );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_Redraw ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BSP_DEGU:
      switch ( obj->bsp_type ) {
    case BSP_TYPE_SPHERICAL:
        equator = (GO_BSplineCurve*)obj->me.dependencies[0];
        if ( GeomObjectBSplineCurveSetDegree ( equator, key ) ) {
          bsp_eqdegree = degreeu = equator->degree;
          Geom10winKNSetKnots ( &g10knotwineqmer, equator->degree,
                              equator->lastknot, equator->knots, equator->closed );
          goto proceed_deguv;
        }
        else
          goto error_deguv;
    default:
        if ( GeomObjectBSplinePatchSetDegreeU ( obj, key ) ) {
          degreeu = obj->degree_u;
          goto proceed_deguv;
        }
        else
          goto error_deguv;
      }
  case intwM1BSP_DEGV:
      switch ( obj->bsp_type ) {
    case BSP_TYPE_SPHERICAL:
        meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
        if ( GeomObjectBSplineCurveSetDegree ( meridian, key ) ) {
          bsp_merdegree = degreev = meridian->degree;
          Geom10winKNSetKnots ( &g10knotwineqmer, meridian->degree,
                                meridian->lastknot, meridian->knots, meridian->closed );
          goto proceed_deguv;
        }
        else
          goto error_deguv;
    default:
        if ( GeomObjectBSplinePatchSetDegreeV ( obj, key ) ) {
          degreev = obj->degree_v;
proceed_deguv:
          Geom10winT2SetKnots ( &g10t2knotwin,
                                obj->degree_u, obj->lastknot_u, obj->knots_u,
                                obj->closed_u,
                                obj->degree_v, obj->lastknot_v, obj->knots_v,
                                obj->closed_v );
          if ( RenderingIsOn )
            StopRendering ();
          rendered_picture = false;
          xge_RedrawAll ();
        }
        else {
error_deguv:
          xge_DisplayErrorMessage ( ErrorMessageCannotCangeDegree, 0 );
        }
      }
      return 1;
  case intwM1BSP_DENSITY_U:
      if ( GeomObjectBSplinePatchSetDensityU ( obj, key ) ) {
        density_u = obj->dens_u;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_DENSITY_V:
      if ( GeomObjectBSplinePatchSetDensityV ( obj, key ) ) {
        density_v = obj->dens_v;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G1BLENDING_UMIN:
      if ( SetMinRange ( 2, obj->lastknot_u-5, key,
                         &bsp_bl_opt_range[0], &bsp_bl_opt_range[1] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 2, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G1BLENDING_UMAX:
      if ( SetMaxRange ( 2, obj->lastknot_u-5, key,
                         &bsp_bl_opt_range[0], &bsp_bl_opt_range[1] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 2, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G1BLENDING_VMIN:
      if ( SetMinRange ( 2, obj->lastknot_v-5, key,
                         &bsp_bl_opt_range[2], &bsp_bl_opt_range[3] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 2, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G1BLENDING_VMAX:
      if ( SetMaxRange ( 2, obj->lastknot_v-5, key,
                         &bsp_bl_opt_range[2], &bsp_bl_opt_range[3] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 2, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G2BLENDING_UMIN:
      if ( SetMinRange ( 3, obj->lastknot_u-7, key,
                         &bsp_bl_opt_range[0], &bsp_bl_opt_range[1] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 3, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G2BLENDING_UMAX:
      if ( SetMaxRange ( 3, obj->lastknot_u-7, key,
                         &bsp_bl_opt_range[0], &bsp_bl_opt_range[1] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 3, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G2BLENDING_VMIN:
      if ( SetMinRange ( 3, obj->lastknot_v-7, key,
                         &bsp_bl_opt_range[2], &bsp_bl_opt_range[3] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 3, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_G2BLENDING_VMAX:
      if ( SetMaxRange ( 3, obj->lastknot_v-7, key,
                         &bsp_bl_opt_range[2], &bsp_bl_opt_range[3] ) ) {
        bsp_blending_entire = BlBSplineEntire ( 3, obj->lastknot_u, obj->lastknot_v,
                                                bsp_bl_opt_range );
        memcpy ( obj->blp_range, bsp_bl_opt_range, 4*sizeof(int) );
        GeomObjectBSplinePatchMarkBLRange ( obj );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSP_NKN1:
      obj->nkn1 = bsm_bl_nkn1 = key;
      if ( bsm_bl_nkn2 < key )
        obj->nkn2 = bsm_bl_nkn2 = key;
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case intwM1BSP_NKN2:
      obj->nkn2 = bsm_bl_nkn2 = key;
      if ( bsm_bl_nkn1 > key )
        obj->nkn1 = bsm_bl_nkn1 = key;
        xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case intwM1BSP_MAXIT:
      obj->maxit = bsm_bl_maxit = key;
      return 1;
  case intWM1BSP_SPRODUCT_EQDEG:
      equator = (GO_BSplineCurve*)current_go->dependencies[0];
      if ( GeomObjectBSplineCurveSetDegree ( equator, key ) ) {
        bsp_eqdegree = key;
        Geom10winKNSetKnots ( &g10knotwineqmer, equator->degree, equator->lastknot,
                              equator->knots, equator->closed );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case intWM1BSP_SPRODUCT_MERDEG:
      meridian = (GO_BSplineCurve*)current_go->dependencies[1];
      if ( GeomObjectBSplineCurveSetDegree ( meridian, key ) ) {
        bsp_merdegree = key;
        Geom10winKNSetKnots ( &g10knotwineqmer, meridian->degree, meridian->lastknot,
                              meridian->knots, meridian->closed );
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    switch ( er->id ) {
  case slM1BSP_BLENDING_CPARAM:
      obj->blp_C = xge_LogSlidebarValued ( BSP_BL_MIN_CPARAM, BSP_BL_MAX_CPARAM,
                                           sl_bsm_bl_param );
      NotifyParam2 ( obj->blp_C );
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BSP_NAME:
      return 1;
  case textedM1BSP_SPRODUCT_EQNAME:
      return 1;
  case textedM1BSP_SPRODUCT_MERNAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BSP_NAME:
      return 1;
  case textedM1BSP_SPRODUCT_EQNAME:
      equator = GeomObjectFindBSCurve2DByName ( obj->eqname );
      if ( equator ) {
        saveeqmer = (GO_BSplineCurve*)obj->me.dependencies[0];
        obj->me.dependencies[0] = (geom_object*)equator;
        GeomObjectSortDependencies ();
        if ( GeomObjectBSplinePatchGenSphericalProduct ( obj ) ) {
          SetupBSplinePatchWidgets ( obj );
          if ( RenderingIsOn )
            StopRendering ();
          rendered_picture = false;
          xge_RedrawAll ();
        }
        else {
          obj->me.dependencies[0] = (geom_object*)saveeqmer;
          GeomObjectSortDependencies ();
          goto restore_eq;
        }
      }
      else {
restore_eq:
        equator = (GO_BSplineCurve*)obj->me.dependencies[0];
        strcpy ( obj->eqname, equator->me.name );
        InitNameEditor ( &bsp_sproduct_eqname_ed, obj->eqname );
        xge_DisplayErrorMessage ( ErrorMsgCannotFindObject, 0 );
      }
      return 1;
  case textedM1BSP_SPRODUCT_MERNAME:
      meridian = GeomObjectFindBSCurve2DByName ( obj->mername );
      if ( meridian ) {
        saveeqmer = (GO_BSplineCurve*)obj->me.dependencies[1];
        obj->me.dependencies[1] = (geom_object*)meridian;
        GeomObjectSortDependencies ();
        if ( GeomObjectBSplinePatchGenSphericalProduct ( obj ) ) {
          SetupBSplinePatchWidgets ( obj );
          if ( RenderingIsOn )
            StopRendering ();
          rendered_picture = false;
          xge_RedrawAll ();
        }
        else {
          obj->me.dependencies[1] = (geom_object*)saveeqmer;
          GeomObjectSortDependencies ();
          goto restore_mer;
        }
      }
      else {
restore_mer:
        meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
        strcpy ( obj->mername, meridian->me.name );
        InitNameEditor ( &bsp_sproduct_mername_ed, obj->mername );
        xge_DisplayErrorMessage ( ErrorMsgCannotFindObject, 0 );
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case textedM1BSP_NAME:
      return 1;
  case textedM1BSP_SPRODUCT_EQNAME:
      equator = (GO_BSplineCurve*)obj->me.dependencies[0];
      strcpy ( obj->eqname, equator->me.name );
      InitNameEditor ( &bsp_sproduct_eqname_ed, obj->eqname );
      return 1;
  case textedM1BSP_SPRODUCT_MERNAME:
      meridian = (GO_BSplineCurve*)obj->me.dependencies[1];
      strcpy ( obj->mername, meridian->me.name );
      InitNameEditor ( &bsp_sproduct_mername_ed, obj->mername );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side10MenuBspCallBack*/

