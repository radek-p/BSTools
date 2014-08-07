
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
#include "mengerc.h"
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

#define PARENT_SIDE
#include "pozwalajipc.h"

#define BSC_MC_MAXITER 999


void InitSide10Menu_BSc ( void )
{
  xge_widget *w;
  int        i;

        /* widgets specific for B-spline curves */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBSplineCurve );
  w = xge_NewStringEd ( win1, w, textedM1BSC_NAME, 109, 19, 0, 40,
                        MAX_NAME_LENGTH, objectname, &bsc_name_ed );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_DEG, 109, 19, 0, 60,
                         0, MAX_DEGREE+1, &intw_bscdeg, txtDegree, &degree );
  w = xge_NewButton ( win1, w, btnM1BSC_UNIFORM, 76, 19, 0, 80, txtUniform );
  w = xge_NewButton ( win1, w, btnM1BSC_REFINE, 76, 19, 0, 100, txtRefine );
  w = xge_NewSwitch ( win1, w, swM1BSC_CLOSED, 109, 16, 0, 120,
                      txtClosed, &bsc_sw_closed );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsc_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsc_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_edit = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, swM1BSC_VIEW_CURVE, 109, 16, 0, 22,
                      txtCurve, &sw_view_curve );
  w = xge_NewSwitch ( win1, w, swM1BSC_VIEW_CPOLY, 109, 16, 0, 42,
                      txtControlPolygon, &sw_view_cpoly );
  w = xge_NewSwitch ( win1, w, swM1BSC_VIEW_BPOLY, 109, 16, 0, 62,
                      txtBezierPolygons, &bsc_sw_view_bpoly );
  w = xge_NewSwitch ( win1, w, swM1BSC_VIEW_CURVATURE, 109, 16, 0, 82,
                      txtCurvatureGraph, &sw_view_curvature );
  w = xge_NewSlidebard ( win1, w, slM1BSC_SCALE_CURVATURE, 109, 10, 0, 102,
                         &curvature_scale );
  w = xge_NewSwitch ( win1, w, swM1BSC_VIEW_TORSION, 109, 16, 0, 116,
                      txtTorsionGraph, &sw_view_torsion );
  w = xge_NewSlidebard ( win1, w, slM1BSC_SCALE_TORSION, 109, 10, 0, 136,
                         &torsion_scale );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_GRAPH_DENSITY, 109, 19, 0, 150,
                         1, 32, &bsc_graphdens, txtGraphDensity, &curv_graph_dens );
  w = xge_NewButton ( win1, w, btnM1BSC_COLOUR, 60, 18, 0, 172, txtColour );
  w = xge_NewTextWidget ( win1, w, 0, 109, 19, 0, 192, txtPipeDiameter );
  w = xge_NewSlidebard ( win1, w, slM1BSC_PIPE_DIAMETER, 109, 10, 0, 212,
                         &sl_pipe_diameter );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsc_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsc_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_view = w;
          /* data */
  w = xge_NewSwitch ( win1, NULL, swM1BSC_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsc_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsc_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_data = w;
          /* options */
            /* general */
  w = xge_NewSwitch ( win1, NULL, swM1MSC_MENGERC, 109, 16, 0, 20,
                      txtMengerCurv, &bsc_sw_mengerc );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsc_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsc_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_opt = side10wdg_bsc_opt1 = w;
            /* Menger curvature */
  w = xge_NewSwitch ( win1, NULL, swM1MSC_MENGERC, 109, 16, 0, 20,
                      txtMengerCurv, &bsc_sw_mengerc );
  w = xge_NewTextWidget ( win1, w, 0, 109, 19, 0, 40, txtExponent );
  w = xge_NewSlidebard ( win1, w, slM1MSC_MENGERC_EXP, 109, 10, 0, 60,
                         &bsc_sl_mcexp );
  w = xge_NewTextWidget ( win1, w, 0, 109, 19, 0, 74, txtPenaltyParam );
  for ( i = 0; i < 5; i++ )
    w = xge_NewSlidebard ( win1, w, slM1MSC_MENGERC_P1+i, 109, 10, 0, 94+14*i,
                           &bsc_sl_mcppar[i] );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_MENGERC_QKN, 91, 19, 0, 164,
                         MENGERC_MIN_NQKN, MENGERC_MAX_NQKN, &bsc_mc_qkn,
                         txtQKnots, &bsc_mc_qknots );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_MENGERC_POPT, 91, 19, 0, 184,
                         0, 3, &bsc_mc_popt, txtPenaltyOpt, &bsc_mc_ppopt );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_MENGERC_MAXIT, 91, 19, 0, 204,
                         1, BSC_MC_MAXITER, &bsc_mc_maxit, txtMaxIter,
                         &bsc_mc_maxiter );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_MENGERC_NPTHREADS, 91, 19, 0, 224,
                         1, MAX_PTHREADS, &bscnpthreads,
                         txtNPThreads, &bsc_npthreads );
  w = xge_NewSwitch ( win1, w, swM1BSC_MENGERC_LOG, 109, 16, 0, 245,
                      txtLogIt, &bsc_sw_mc_logit );
  w = bsc_mc_optimize =
      xge_NewButton ( win1, w, btnM1BSC_MENGERC_OPTIMIZE, 58, 19, 0, 263,
                      txtOptimize );
  bsc_npthreads = ncpu;
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsc_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSC_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsc_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_opt2 = w;
} /*InitSide10Menu_BSc*/

boolean ChangeSide10MenuWidth_BSc ( short h )
{
  boolean result;

  result = side10menu_wide;
  side10menu_wide = false;
  return result;
} /*ChangeSide10MenuWidth_BSc*/

void SetupBSplineCurveWidgets ( GO_BSplineCurve *obj )
{
  int i;

  InitNameEditor ( &bsc_name_ed, obj->me.name );
  degree = obj->degree;
  bsc_sw_closed = obj->closed;
  sw_view_curve = obj->view_curve;
  sw_view_cpoly = obj->view_cpoly;
  bsc_sw_view_bpoly = obj->view_bpoly;
  sw_view_curvature = obj->view_curvature;
  sw_view_torsion = obj->view_torsion;
  curvature_scale = obj->curvature_scale;
  torsion_scale = obj->torsion_scale;
  curv_graph_dens = obj->graph_dens;
  sl_pipe_diameter = xge_LogSlidebarPosd ( BSC_MIN_PIPE_DIAMETER,
                                           BSC_MAX_PIPE_DIAMETER,
                                           obj->pipe_diameter );
  Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                        obj->knots, obj->closed );
  SetGeomWin10Knotw ( obj );
  sw_bsc_dom_coord = g10knotwin.display_coord;
  sw_bsc_dom_panzoom = g10knotwin.panning;
        /* Menger curvature related stuff */
  bsc_sw_mengerc = obj->mengerc;
  bsc_sl_mcexp = xge_LogSlidebarPosd ( BSC_MC_MIN_EXP, BSC_MC_MAX_EXP,
                                       obj->mc_exponent );
  for ( i = 0; i < 5; i++ )
    bsc_sl_mcppar[i] = xge_LogSlidebarPosd ( BSC_MC_MIN_PPAR, BSC_MC_MAX_PPAR,
                                             obj->mc_pparam[i] );
  bsc_mc_qknots = obj->mc_nqkn;
  bsc_mc_ppopt = obj->mc_ppopt;
  if ( bsc_sw_mengerc ) {
    if ( obj->me.bound_with_a_child )
      bsc_mc_optimize->data0 = txtInterrupt;
    else
      bsc_mc_optimize->data0 = txtOptimize;
    side10wdg_bsc_opt = side10wdg_bsc_opt2;
    xge_KnotWindFindMapping ( &g10knotwin );
  }
  else
    side10wdg_bsc_opt = side10wdg_bsc_opt1;
  if ( whichside10menu == SIDE10MENU_OPTIONS )
    xge_SetMenuWidgets ( side10menu, side10wdg_bsc_opt, true );
} /*SetupBSplineCurveWidgets*/

int Side10MenuBscCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplineCurve *obj;
  int             i;
  boolean         success;

  if ( current_go->obj_type != GO_BSPLINE_CURVE )
    return 0;
  obj = (GO_BSplineCurve*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BSC_UNIFORM:
      if ( GeomObjectBSplineCurveSetUniformKnots ( obj, true ) ) {
        xge_KnotWindFindMapping ( &g10knotwin );
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSC_REFINE:
      if ( GeomObjectBSplineCurveRefine ( obj ) ) {
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_KnotWindFindMapping ( &g10knotwin );
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSC_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  case btnM1BSC_MENGERC_OPTIMIZE:
      if ( ipc_state == ipcstate_CHILD_BUSY ) {
        if ( obj->me.bound_with_a_child ) {
          IPCInterruptTheChild ();
          bsc_mc_optimize->data0 = txtOptimize;
        }
        else {
          xge_DisplayErrorMessage ( ErrorMsgChildProcessBusy, 0 );
          return 1;
        }
      }
      else {
        if ( MengerCurvatureOptimizationPrepareData ( obj ) ) {
          InitMengerCurvatureOptimization ();
          bsc_mc_optimize->data0 = txtInterrupt;
        }
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM1BSC_VIEW_CURVE:
      obj->view_curve = sw_view_curve;
      RedrawGeom00Win ();
      return 1;
  case swM1BSC_VIEW_CPOLY:
      obj->view_cpoly = sw_view_cpoly;
      RedrawGeom00Win ();
      return 1;
  case swM1BSC_CLOSED:
      if ( GeomObjectBSplineCurveSetClosed ( obj, bsc_sw_closed ) ) {
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_RedrawAll ();
      }
      else if ( bsc_sw_closed )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      bsc_sw_closed = obj->closed;
      return 1;
  case swM1BSC_DOMAIN_COORD:
      g10knotwin.display_coord = sw_bsc_dom_coord;
      return 1;
  case swM1BSC_DOMAIN_PANZOOM:
      g10knotwin.panning = sw_bsc_dom_panzoom;
      return 1;
  case swM1BSC_VIEW_BPOLY:
      obj->view_bpoly = bsc_sw_view_bpoly;
      RedrawGeom00Win ();
      return 1;
  case swM1BSC_VIEW_CURVATURE:
      obj->view_curvature = sw_view_curvature;
      if ( sw_view_curvature )
        GeomObjectBSplineCurveSetCurvatureGraph ( obj,
                        sw_view_curvature, curvature_scale,
                        sw_view_torsion, torsion_scale, curv_graph_dens );
      RedrawGeom00Win ();
      return 1;
  case swM1BSC_VIEW_TORSION:
      if ( sw_view_torsion ) {
        if ( obj->me.spdimen == 3 ) {
          GeomObjectBSplineCurveSetCurvatureGraph ( obj,
                          sw_view_curvature, curvature_scale,
                          sw_view_torsion, torsion_scale, curv_graph_dens );
          RedrawGeom00Win ();
        }
        else
          sw_view_torsion = false;
      }
      else {
        obj->view_torsion = false;
        RedrawGeom00Win ();
      }
      return 1;
  case swM1MSC_MENGERC:
      success = GeomObjectBSplineCurveSetMengerc ( obj, bsc_sw_mengerc );
      if ( !success )
        bsc_sw_mengerc = obj->mengerc;
      SetupBSplineCurveWidgets ( obj );
      xge_RedrawAll ();
      if ( !success && !bsc_sw_mengerc ) {
        if ( obj->degree < 3 || !obj->closed )
          xge_DisplayErrorMessage ( ErrorMsgCurveMustBeCubicAndClosed, 0 );
        else if ( obj->rational || obj->me.spdimen != 3 )
          xge_DisplayErrorMessage ( ErrorMsgCurveMustBeNonRationalAnd3D, 0 );
      }
      return 1;
  case swM1BSC_MENGERC_LOG:
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    switch ( er->id ) {
  case slM1BSC_SCALE_CURVATURE:
      if ( sw_view_curvature ) {
        GeomObjectBSplineCurveSetCurvatureGraph ( obj,
                        sw_view_curvature, curvature_scale,
                        sw_view_torsion, torsion_scale, curv_graph_dens );
        RedrawGeom00Win ();
      }
      return 1;
  case slM1BSC_SCALE_TORSION:
      if ( sw_view_torsion ) {
        GeomObjectBSplineCurveSetCurvatureGraph ( obj,
                        sw_view_curvature, curvature_scale,
                        sw_view_torsion, torsion_scale, curv_graph_dens );
        RedrawGeom00Win ();
      }
      return 1;
  case slM1BSC_PIPE_DIAMETER:
      obj->pipe_diameter = xge_LogSlidebarValued ( BSC_MIN_PIPE_DIAMETER,
                                  BSC_MAX_PIPE_DIAMETER, sl_pipe_diameter );
      NotifyParam2 ( obj->pipe_diameter );
      return 1;
  case slM1MSC_MENGERC_EXP:
      obj->mc_exponent = xge_LogSlidebarValued ( BSC_MC_MIN_EXP,
                                 BSC_MC_MAX_EXP, bsc_sl_mcexp );
      NotifyParam2 ( obj->mc_exponent );
      return 1;
  case slM1MSC_MENGERC_P1:
  case slM1MSC_MENGERC_P2:
  case slM1MSC_MENGERC_P3:
  case slM1MSC_MENGERC_P4:
  case slM1MSC_MENGERC_P5:
      i = er->id-slM1MSC_MENGERC_P1;
      obj->mc_pparam[i] = xge_LogSlidebarValued ( BSC_MC_MIN_PPAR,
                                 BSC_MC_MAX_PPAR, bsc_sl_mcppar[i] );
      NotifyParam2 ( obj->mc_pparam[i] );
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BSC_DEG:
      if ( GeomObjectBSplineCurveSetDegree ( obj, key ) ) {
        degree = obj->degree;
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSC_GRAPH_DENSITY:
      curv_graph_dens = key;
      if ( sw_view_curvature || (sw_view_torsion && obj->me.spdimen == 3) ) {
        GeomObjectBSplineCurveSetCurvatureGraph ( obj,
                        sw_view_curvature, curvature_scale,
                        sw_view_torsion, torsion_scale, curv_graph_dens );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSC_MENGERC_QKN:
      obj->mc_nqkn = bsc_mc_qknots = key;
      return 1;
  case intwM1BSC_MENGERC_POPT:
      obj->mc_ppopt = bsc_mc_ppopt = key;
      return 1;
  case intwM1BSC_MENGERC_MAXIT:
      bsc_mc_maxiter = key;
      return 1;
  case intwM1BSC_MENGERC_NPTHREADS:
      bsc_npthreads = key;
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BSC_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BSC_NAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side10MenuBscCallBack*/

/* ////////////////////////////////////////////////////////////////////////// */
/* data buffers for Menger curvature optimization */
static ipc_bscmc_size    bsc_size;
static ipc_bscmc_options bsc_options;

boolean MengerCurvatureOptimizationPrepareData ( GO_BSplineCurve *obj )
{
  int cpdimen, lkn, deg;

  BindChildToGeomObject ( (geom_object*)obj );
  ResetIPCBuffer ();
  bsc_size.spdimen = obj->me.spdimen;
  bsc_size.cpdimen = cpdimen = obj->me.cpdimen;
  bsc_size.degree = deg = obj->degree;
  bsc_size.lkn = lkn = obj->lastknot;
  bsc_size.closed = obj->closed;
  IPCAppendDataItem ( ipcd_BSC_SIZE, sizeof(ipc_bscmc_size), &bsc_size );
  IPCAppendDataItem ( ipcd_BSC_KNOTS, (lkn+1)*sizeof(double), obj->knots );
  IPCAppendDataItem ( ipcd_BSC_CPOINTS, (lkn-deg)*cpdimen*sizeof(double),
                      obj->cpoints );
  IPCAppendDataItem ( ipcd_BSC_MKCP, (lkn-deg)*sizeof(char), obj->mkcp );
  bsc_options.nqkn = obj->mc_nqkn;
  bsc_options.maxiter = bsc_mc_maxiter;
  bsc_options.ppopt = obj->mc_ppopt;
  bsc_options.exponent = obj->mc_exponent;
  memcpy ( bsc_options.pparam, obj->mc_pparam, 5*sizeof(double) );
  bsc_options.npthreads = bsc_npthreads;
  IPCAppendDataItem ( ipcd_BSC_OPTIMIZE, sizeof(ipc_bscmc_options),
                      &bsc_options );
  return true;
} /*MengerCurvatureOptimizationPrepareData*/

void InitMengerCurvatureOptimization ( void )
{
  switch ( ipc_state ) {
case ipcstate_NO_CHILD: 
    if ( LaunchAChildProcess () )
      xge_PostIdleCommand ( IDLE_COMMAND_BSC_MENGERC_OPT_INIT, 0, 0 );
    else
      xge_DisplayErrorMessage ( ErrorMsgCannotLaunchAChild, 0 );
    break;
case ipcstate_CHILD_LAUNCHED:
    xge_PostIdleCommand ( IDLE_COMMAND_BSC_MENGERC_OPT_INIT, 0, 0 );
    break;
case ipcstate_CHILD_READY:
    IPCWakeUpChild ();
    break;
case ipcstate_CHILD_BUSY:
    xge_DisplayErrorMessage ( ErrorMsgChildProcessBusy, 0 );
    break;
default:
    break;
  }
} /*InitMengerCurvatureOptimization*/

