
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
#include <pthread.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/gl.h>  
#include <GL/glu.h> 
#include <GL/glx.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"
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

#define BSM_BL_MAXITER 999

void InitSide10Menu_BSm ( void )
{
  xge_widget *w, *ww, *sw;

        /* widgets specific for B-spline meshes */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 0, txtBSplineMesh );
  w = xge_NewStringEd ( win1, w, textedM1BSM_NAME, 109, 19, 0, 20,
                        MAX_NAME_LENGTH, objectname, &bsm_name_ed );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_DEG, 109, 19, 0, 40,
                         0, MAX_DEGREE+1, &intw_bsmdeg, txtDegree, &degree );
  w = xge_NewButton ( win1, w, btnM1BSM_REFINEMENT, 79, 19, 0, 60, txtRefine );
  w = xge_NewButton ( win1, w, btnM1BSM_DOUBLING, 79, 19, 0, 80, txtDouble );
  w = xge_NewButton ( win1, w, btnM1BSM_AVERAGING, 79, 19, 0, 100, txtAverage );
  w = xge_NewButton ( win1, w, btnM1BSM_EXTRACTSUBMESH, 79, 19, 0, 120,
                      txtGetSubmesh );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_VERTEX0, 75, 19, 0, 144,
                         -2, 10, &intw_bsm_vert0, txtVertex, &bsm_vertex_num0 );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_VERTEX1, 35, 19, 74, 144,
                         -2, 10, &intw_bsm_vert1, NULL, &bsm_vertex_num1 );
  w = xge_NewButton ( win1, w, btnM1BSM_MARK_VERT, 27, 19, 0, 164,
                      txtMark );
  w = xge_NewButton ( win1, w, btnM1BSM_UNMARK_VERT, 39, 19, 28, 164,
                      txtUnmark );
  w = xge_NewTextWidget ( win1, w, 0, 45, 19, 70, 164, txtVert );
  w = xge_NewButton ( win1, w, btnM1BSM_MARK_HEDGE, 27, 19, 0, 184,
                      txtMark );
  w = xge_NewButton ( win1, w, btnM1BSM_UNMARK_HEDGE, 39, 19, 28, 184,
                      txtUnmark );
  w = xge_NewTextWidget ( win1, w, 0, 45, 19, 70, 184, txtHEdg );
  w = xge_NewButton ( win1, w, btnM1BSM_FILTER, 39, 19, 0, 204, txtFilter );
  w = xge_NewButton ( win1, w, btnM1BSM_ENTER_LINE, 39, 19, 40, 204, txtLine );
  w = xge_NewButton ( win1, w, btnM1BSM_REMOVE_VERTEX, 79, 19, 0, 224,
                      txtRemove );
  w = xge_NewButton ( win1, w, btnM1BSM_DIVIDE_FACET, 79, 19, 0, 244,
                      txtDivideFacet );
  w = xge_NewButton ( win1, w, btnM1BSM_DOUBLE_LOOP, 79, 19, 0, 264,
                      txtDoubleLoop );
  w->state = xgestate_BUTTON_INACTIVE;
  w = xge_NewIntWidget ( win1, w, intwM1BSM_EDGE0, 75, 19, 0, 288,
                         -2, 10, &intw_bsm_edge0, txtEdge, &bsm_edge_num0 );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_EDGE1, 35, 19, 74, 288,
                         -2, 10, &intw_bsm_edge1, NULL, &bsm_edge_num1 );
  w = xge_NewButton ( win1, w, btnM1BSM_SHRINK_EDGE, 79, 19, 0, 308,
                      txtShrink );
  w = xge_NewButton ( win1, w, btnM1BSM_CONTRACT_EDGE, 79, 19, 0, 328,
                      txtContract );
  w = xge_NewButton ( win1, w, btnM1BSM_GLUE_EDGES, 79, 19, 0, 348,
                      txtGlueEdges );
  w = xge_NewButton ( win1, w, btnM1BSM_GLUE_EDGE_LOOPS, 79, 19, 0, 368,
                      txtGlueLoops );
  w = xge_NewButton ( win1, w, btnM1BSM_SEAL_HOLE, 79, 19, 0, 388,
                      txtSealHole );
  w = xge_NewButton ( win1, w, btnM1BSM_SPLIT_BOUNDARY_EDGE, 79, 19, 0, 408,
                      txtSplitEdge );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_FACET0, 75, 19, 0, 432,
                         -2, 10, &intw_bsm_fac0, txtFacet, &bsm_facet_num0 );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_FACET1, 35, 19, 74, 432,
                         -2, 10, &intw_bsm_fac1, NULL, &bsm_facet_num1 );
  w = xge_NewButton ( win1, w, btnM1BSM_REMOVE_FACET, 79, 19, 0, 452,
                      txtRemove );
  w = xge_NewButton ( win1, w, btnM1BSM_DOUBLE_FAC_EDGES, 79, 19, 0, 472,
                      txtDoubleEdges );
  for ( ww = w; ww; ww = ww->prev )
    xge_SetWidgetPositioning ( ww, 0, ww->x, ww->y );
  side10wdg_bsm_editcontents = xge_NewMenu ( win1, NULL, scwM1BSM_ECONTENTS,
                      SIDEMENUWIDTH0, 495, 0, 20, w );
  side10wdg_bsm_editscroll = xge_NewScrollWidget ( win1, NULL, scwM1BSM_ESCROLL,
                      SIDEMENUWIDTH0, xge_HEIGHT-TOPMENUHEIGHT-20, 0, 20,
                      &side10_bsm_editsw, side10wdg_bsm_editcontents );
  w = xge_NewSwitch ( win1, side10wdg_bsm_editscroll, swM11STATUS,
                      16, 16, 0, xge_HEIGHT-16, txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_edit = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, swM1BSM_VIEW_CNET, 109, 16, 0, 22,
                      txtControlNet, &sw_view_cnet );
  w = xge_NewSwitch ( win1, w, swM1BSM_VIEW_SURF, 109, 16, 0, 42,
                      txtSurface, &sw_view_surf );
  w = xge_NewSwitch ( win1, w, swM1BSM_HOLE_FILLING, 109, 16, 0, 62,
                      txtHoleFilling, &sw_view_hole_filling );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_DENSITY, 109, 19, 0, 82,
                         1, MAX_PNET_DENSITY, &intwdensityu, txtDensity,
                         &density_u );
  w = xge_NewSwitch ( win1, w, swM1BSM_VIEW_SPECIAL, 109, 16, 0, 105,
                      txtSpecialNets, &bsm_sw_view_special );
  w = xge_NewButton ( win1, w, btnM1BSM_COLOUR, 60, 18, 0, 125, txtColour );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_view = w;
          /* data */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBSplineMesh );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_DATA_KGON, 16, 19, 0, 40,
                         3, MAX_BSM_DEGREE, &intw_bsmdeg1, txtNull, &bsm_degree1 );
  w = xge_NewButton ( win1, w, btnM1BSM_DATA_KGON, 75, 19, 15, 40, txt_gon );
  w = xge_NewButton ( win1, w, btnM1BSM_DATA_TETRAHEDRON, 90, 19, 0, 60,
                      txtTetrahedron );
  w = xge_NewButton ( win1, w, btnM1BSM_DATA_CUBE, 90, 19, 0, 80,
                      txtCube );
  w = xge_NewButton ( win1, w, btnM1BSM_DATA_DODECAHEDRON, 90, 19, 0, 100,
                      txtDodecahedron );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_DATA_KPRISM, 16, 19, 0, 120,
                         3, MAX_BSM_DEGREE, &intw_bsmdeg2, txtNull, &bsm_degree2 );
  w = xge_NewButton ( win1, w, btnM1BSM_DATA_KPRISM, 75, 19, 15, 120,
                      txt_gonalPrism );
  w = xge_NewSwitch ( win1, w, swM1BSM_DATA_REPLACE, 109, 16, 0, 142,
                      txtReplace, &bsm_sw_data_replace );
  w = xge_NewSwitch ( win1, w, swM1BSM_DATA_ADD, 109, 16, 0, 162,
                      txtAdd, &bsm_sw_data_add );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_data = w;
          /* options */
            /* general */
  w = xge_NewSwitch ( win1, NULL, swM1BSM_SUBDIVISION, 109, 16, 0, 20,
                      txtSubdivision, &bsm_sw_subdivision );
  w = xge_NewSwitch ( win1, w, swM1BSM_BLENDING, 109, 16, 0, 40,
                      txtBlending, &bsm_sw_blending );
  w = xge_NewSwitch ( win1, w, swM1BSM_G1, 109, 16, 0, 64,
                      txtG1, &sw_hfill_g1 );
  w = xge_NewSwitch ( win1, w, swM1BSM_G2, 109, 16, 0, 84,
                      txtG2, &sw_hfill_g2 );
  w = xge_NewSwitch ( win1, w, swM1BSM_G1_QUASI_G2, 109, 16, 0, 104,
                      txtG1quasiG2, &sw_hfill_g1q2 );
  w = xge_NewSlidebard ( win1, w, slM1BSM_G1Q2_PARAM, 109, 10, 0, 124,
                         &sl_g1q2_param );
  w = xge_NewSwitch ( win1, w, swM1BSM_HOLEFILL_COONS, 109, 16, 0, 138,
                      txtCoons, &sw_hfill_coons );
  w = xge_NewSwitch ( win1, w, swM1BSM_HOLEFILL_BEZIER, 109, 16, 0, 159,
                      txtBezier, &sw_hfill_bezier );
/*  w = xge_NewButton ( win1, w, btnM1BSM_OPTIMIZE_SPECIALS, 76, 19, 0, 179,
                      txtSpecials ); */
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_opt = side10wdg_bsm_opt1 = w;
            /* subdivision */
  w = xge_NewSwitch ( win1, NULL, swM1BSM_SUBDIVISION, 109, 16, 0, 20,
                      txtSubdivision, &bsm_sw_subdivision );
  w = xge_NewSwitch ( win1, w, swM1BSM_BLENDING, 109, 16, 0, 40,
                      txtBlending, &bsm_sw_blending );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_opt2 = w;
            /* blending */
  w = xge_NewSwitch ( win1, NULL, swM1BSM_SUBDIVISION, 109, 16, 0, 20,
                      txtSubdivision, &bsm_sw_subdivision );
  sw = xge_NewSwitch ( win1, w, swM1BSM_BLENDING, 109, 16, 0, 40,
                      txtBlending, &bsm_sw_blending );
  w = xge_NewButton ( win1, NULL, btnM1BSM_PRETRANSFORMATION, 76, 19, 0, 1,
                      txtPreTransf );
  w = xge_NewSlidebard ( win1, w, slM1BSM_BLENDING_CPARAM, 109, 10, 0, 23,
                         &sl_bsm_bl_param );
  w = xge_NewSwitch ( win1, w, swM1BSM_BLENDING_CONSTRAINTS, 109, 16, 0, 36,
                      txtConstraints, &bsm_sw_constr );
  w = xge_NewSwitch ( win1, w, swM1BSM_SHAPE_ONLY, 109, 16, 0, 56,
                      txtShapeOnly, &bsm_sw_shape_only );
  w = xge_NewSwitch ( win1, w, swM1BSM_ALT_MULTILEVEL, 109, 16, 0, 76,
                      txtAltML, &bsm_sw_alt_multilevel );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_NKN1, 91, 19, 0, 96,
                         3, 10, &bsm_bl_qknots1, txtQKnots, &bsm_bl_nkn1 );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_NKN2, 17, 19, 92, 96,
                         3, 10, &bsm_bl_qknots2, txtNull, &bsm_bl_nkn2 );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_NLEVELS, 91, 19, 0, 116,
                         0, G2MBL_MAX_LEVELS, &bsm_bl_nlevels, txtLevels,
                         &bsm_bl_nlev );
  w = xge_NewButton ( win1, w, btnM1SUGGESTLEVELS, 17, 19, 92, 116,
                      txtQuestionmark );
  w = xge_NewButton ( win1, w, btnM1SUGGESTBLOCKS, 17, 19, 92, 136,
                      txtQuestionmark );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_NBLOCKS, 91, 19, 0, 136,
                         1, G2MBL_MAX_BLOCKS, &bsm_bl_nblocks, txtBlocks,
                         &bsm_bl_nbl );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_MAXIT, 91, 19, 0, 156,
                         1, BSM_BL_MAXITER, &bsm_bl_maxiter, txtMaxIter,
                         &bsm_bl_maxit );
  w = xge_NewSwitch ( win1, w, swM1BSM_COARSE_PRECOND, 109, 16, 0, 177,
                      txtUseCoarseMesh, &bsm_sw_use_coarse );
  w = xge_NewStringEd ( win1, w, textedM1BSM_COARSE_NAME, 109, 19, 0, 195,
                        MAX_NAME_LENGTH, objectname, &bsm_coarse_editor );
  w = xge_NewIntWidget ( win1, w, intwM1BSM_NPTHREADS, 91, 19, 0, 215,
                         1, MAX_PTHREADS, &bsmnpthreads,
                         txtNPThreads, &bsm_npthreads );
  bsm_npthreads = ncpu;
  w = xge_NewSwitch ( win1, w, swM1BSM_LOG_IT, 109, 16, 0, 237,
                      txtLogIt, &bsm_sw_log_it );
  w = bsm_bl_optimize =
      xge_NewButton ( win1, w, btnM1BSM_BLENDING_OPTIMIZE, 76, 19, 0, 257,
                      txtOptimize );
  w = xge_NewButton ( win1, w, btnM1BSM_OPTIMIZE_SPECIALS, 76, 19, 0, 277,
                      txtSpecials );
  for ( ww = w; ww; ww = ww->prev )
    xge_SetWidgetPositioning ( ww, 0, ww->x, ww->y );
  side10wdg_bsm_optcontents = xge_NewMenu ( win1, NULL, scwM1BSM_ECONTENTS,
                      SIDEMENUWIDTH0, 300, 0, 20, w );
  side10wdg_bsm_optscroll = xge_NewScrollWidget ( win1, sw, scwM1BSM_ESCROLL,
                      SIDEMENUWIDTH0, xge_HEIGHT-TOPMENUHEIGHT-60, 0, 60,
                      &side10_bsm_optsw, side10wdg_bsm_optcontents );
  w = xge_NewSwitch ( win1, side10wdg_bsm_optscroll, swM11STATUS,
                      16, 16, 0, xge_HEIGHT-16, txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsm_opt3 = w;
} /*InitSide10Menu_BSm*/

boolean ChangeSide10MenuWidth_BSm ( short h )
{
  boolean wide, result;

  wide = false;
  if ( side10menu->data1 == side10wdg_bsm_edit )
    wide = side10_bsm_editsw.contents->h > h-20;
  else if ( side10menu->data1 == side10wdg_bsm_opt3 )
    wide = side10_bsm_optsw.contents->h > h-60;
  else
    wide = false;
  result = wide != side10menu_wide;
  side10menu_wide = wide;
  return result;
} /*ChangeSide10MenuWidth_BSm*/

void SetupBSplineMeshVEFnum ( GO_BSplineMesh *obj )
{
  bsm_vertex_num0 = obj->current_vertex[0];
  bsm_vertex_num1 = obj->current_vertex[1];
  intw_bsm_vert0.maxvalue = intw_bsm_vert1.maxvalue = max ( 10, obj->nv );
  bsm_edge_num0 = obj->current_edge[0];
  bsm_edge_num1 = obj->current_edge[1];
  intw_bsm_edge0.maxvalue = intw_bsm_edge1.maxvalue = max ( 10, obj->nhe );
  bsm_facet_num0 = obj->current_facet[0];
  bsm_facet_num1 = obj->current_facet[1];
  intw_bsm_fac0.maxvalue = intw_bsm_fac1.maxvalue = max ( 10, obj->nfac );
} /*SetupBSplineMeshVEFnum*/

static void SetupMaxIter ( GO_BSplineMesh *obj )
{
  bsm_bl_maxiter.title = txtMaxIter;
  bsm_bl_maxiter.minvalue = 1;
  bsm_bl_maxiter.maxvalue = BSM_BL_MAXITER;
  bsm_bl_maxit = obj->maxit;
} /*SetupMaxIter*/

static void SetupStartFrom ( GO_BSplineMesh *obj )
{
  bsm_bl_maxiter.title = txtStartFrom;
  bsm_bl_maxiter.minvalue = 0;
  bsm_bl_maxiter.maxvalue = obj->nblocks;
  bsm_bl_maxit = obj->startfrom;
} /*SetupStartFrom*/

void SetupBSplineMeshWidgets ( GO_BSplineMesh *obj )
{
  if ( obj->me.obj_type != GO_BSPLINE_MESH )
    return;
  InitNameEditor ( &bsm_name_ed, obj->me.name );
  degree = obj->degree;
  intwdensityu.title = &txtDensity[0];
  bsm_sw_subdivision = obj->subdivision;
  bsm_sw_blending = obj->blending;
  if ( bsm_sw_subdivision ) {
    intwdensityu.title = &txtLevel[0];
    density_u = obj->subdivl;
    side10wdg_bsm_opt = side10wdg_bsm_opt2;
  }
  else {
    intwdensityu.title = &txtDensity[0];
    density_u = obj->density;
    if ( bsm_sw_blending )
      side10wdg_bsm_opt = side10wdg_bsm_opt3;
    else
      side10wdg_bsm_opt = side10wdg_bsm_opt1;
  }
  sw_view_surf = obj->view_surf;
  sw_view_cnet = obj->view_cnet;
  bsm_sw_view_special = obj->view_special;
  sw_view_hole_filling = obj->view_holefill;
  sw_hfill_g1 = obj->fill_G1;
  sw_hfill_g2 = obj->fill_G2;
  sw_hfill_g1q2 = obj->fill_G1Q2;
  sw_hfill_coons = obj->fill_Coons;
  sw_hfill_bezier = obj->fill_Bezier;
  sl_g1q2_param = obj->sl_g1q2param;
  bsm_bl_nkn1 = obj->nkn1;
  bsm_bl_nkn2 = obj->nkn2;
  bsm_bl_nlev = obj->nlevels;
  bsm_bl_nbl = obj->nblocks;
  sl_bsm_bl_param = log ( obj->bsm_bl_C/BSM_BL_MIN_CPARAM )/
                    log ( BSM_BL_MAX_CPARAM/BSM_BL_MIN_CPARAM );
  bsm_sw_constr = obj->bl_constr;
  bsm_sw_shape_only = obj->bl_shape_only;
  bsm_sw_use_coarse = obj->bl_use_coarse;
  InitNameEditor ( &bsm_coarse_editor, obj->coarse_name );
  if ( bsm_sw_blending ) {
    if ( obj->me.bound_with_a_child )
      bsm_bl_optimize->data0 = txtInterrupt;
    else
      bsm_bl_optimize->data0 = txtOptimize;
    if ( obj->nlevels )
      SetupStartFrom ( obj );
    else
      SetupMaxIter ( obj );
  }
  SetupBSplineMeshVEFnum ( obj );
  SetGeomWin10Empty ();
  bsm_bl_nlevels.maxvalue = G2MBL_MAX_LEVELS;
} /*SetupBSplineMeshWidgets*/

void SuggestOptLevels ( GO_BSplineMesh *obj, boolean multilevel )
{
  void    *sp;
  int     i, nv, nvcp, minl, maxl;
  byte    *mkcp;
  boolean res;

  sp = pkv_GetScratchMemTop ();
  nv = obj->nv;
  if ( bsm_sw_constr ) {
    mkcp = pkv_GetScratchMem ( nv );
    if ( !mkcp )
      goto failure;
    memcpy ( mkcp, obj->mkcp, nv );
    for ( i = 0; i < nv; i++ )
      mkcp[i] &= marking_mask;
  }
  else
    mkcp = NULL;
  if ( multilevel ) {
    if ( obj->bl_shape_only ) {
      if ( obj->bl_use_coarse )
        res = g2mbl_MLCPSSuggestNLevels ( obj->nv, obj->meshv, obj->meshvhei,
                                  obj->nhe, obj->meshhe,
                                  obj->nfac, obj->meshfac, obj->meshfhei,
                                  mkcp, &minl, &maxl );
      else
        res = g2mbl_MLSSuggestNLevels ( obj->nv, obj->meshv, obj->meshvhei,
                                  obj->nhe, obj->meshhe,
                                  obj->nfac, obj->meshfac, obj->meshfhei,
                                  mkcp, &minl, &maxl );
    }
    else {
      if ( obj->bl_use_coarse )
        res = g2mbl_MLCPSuggestNLevels ( obj->nv, obj->meshv, obj->meshvhei,
                                  obj->nhe, obj->meshhe,
                                  obj->nfac, obj->meshfac, obj->meshfhei,
                                  mkcp, &minl, &maxl );
      else
        res = g2mbl_MLSuggestNLevels ( obj->nv, obj->meshv, obj->meshvhei,
                                  obj->nhe, obj->meshhe,
                                  obj->nfac, obj->meshfac, obj->meshfhei,
                                  mkcp, &minl, &maxl );
    }
    if ( res ) {
      bsm_bl_nlevels.maxvalue = maxl;
      bsm_bl_nbl = obj->nblocks = (0x00001 << minl) - 1;
      bsm_bl_nblocks.maxvalue = bsm_bl_nbl+1;
      bsm_bl_nlev = obj->nlevels = minl;
      obj->startfrom = bsm_bl_nbl-1;
      SetupStartFrom ( obj );
    }
  }
  else {
        /* a primitive method, it does not take into account the */
        /* number of boundary vertices */
    if ( mkcp ) {
      for ( i = nvcp = 0;  i < nv;  i++ )
        if ( mkcp[i] & marking_mask )
          nvcp ++;
    }
    else
      nvcp = nv;
    bsm_bl_nbl = obj->nblocks = nvcp / 2500 + 1;
    bsm_bl_nblocks.maxvalue = bsm_bl_nbl+1;
    bsm_bl_nlev = obj->nlevels = 0;
    SetupMaxIter ( obj );
  }
  xge_SetClipping ( side10menu );
  side10menu->redraw ( side10menu, true );
failure:
  pkv_SetScratchMemTop ( sp );
} /*SuggestOptLevels*/

void BlendingMeshOptSpecialPatches ( void )
{
  GO_BSplineMesh *obj;

  if ( current_go && current_go->obj_type == GO_BSPLINE_MESH ) {
    obj = (GO_BSplineMesh*)current_go;
    GeomObjectBSplineMeshOptSpecialPatches ( obj );
    rendered_picture = false;
    xge_SetWindow ( win0 );
    xge_Redraw ();
  }
  xge_ReleaseFocus ( xge_null_widget );
} /*BlendingMeshOptSpecialPatches*/

void Side10MenuBsmResize ( short x, short y )
{
  if ( side10menu->data1 == side10wdg_bsm_edit )
    side10wdg_bsm_editscroll->msgproc ( side10wdg_bsm_editscroll,
                                        xgemsg_RESIZE, 0, x, y-20 );
  else if ( side10menu->data1 == side10wdg_bsm_opt3 )
    side10wdg_bsm_optscroll->msgproc ( side10wdg_bsm_optscroll,
                                       xgemsg_RESIZE, 0, x, y-60 );
} /*Side10MenuBsmResize*/

int Side10MenuBsmCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplineMesh *obj;
  point3d        p[BSM_NCURRENT_VERTICES];

  if ( current_go->obj_type != GO_BSPLINE_MESH )
    return 0;
  obj = (GO_BSplineMesh*)current_go;
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM1BSM_VIEW_SURF:
      obj->view_surf = sw_view_surf;
      RedrawGeom00Win ();
      return 1;
  case swM1BSM_VIEW_CNET:
      obj->view_cnet = sw_view_cnet;
      RedrawGeom00Win ();
      return 1;
  case swM1BSM_HOLE_FILLING:
      obj->view_holefill = sw_view_hole_filling;
      if ( obj->degree == 3 ) {
        rendered_picture = false;
        RedrawGeom00Win ();
      }
      return 1;
  case swM1BSM_VIEW_SPECIAL:
      obj->view_special = bsm_sw_view_special;
      RedrawGeom00Win ();
      return 1;
  case swM1BSM_DATA_REPLACE:
      bsm_sw_data_add = !bsm_sw_data_replace;
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case swM1BSM_DATA_ADD:
      bsm_sw_data_replace = !bsm_sw_data_add;
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case swM1BSM_SUBDIVISION:
      if ( GeomObjectBSplineMeshSetSubdiv ( obj, bsm_sw_subdivision ) ) {
        rendered_picture = false;
        RedrawGeom00Win ();
      }
      else
        bsm_sw_subdivision = obj->subdivision;
      if ( bsm_sw_subdivision ) {
        if ( bsm_sw_blending ) {
          bsm_sw_blending = false;
        }
        side10wdg_bsm_opt = side10wdg_bsm_opt2;
      }
      else
        side10wdg_bsm_opt = side10wdg_bsm_opt1;
      xge_SetMenuWidgets ( side10menu, side10wdg_bsm_opt, false );
      bsm_sw_blending = obj->blending;
      if ( ChangeSide10MenuWidth ( side00menu->h ) )
        ResizeWindow1 ( xge_current_width, xge_current_height );
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_BLENDING:
      if ( GeomObjectBSplineMeshSetBlending ( obj, bsm_sw_blending ) ) {
        rendered_picture = false;
        RedrawGeom00Win ();
      }
      else
        bsm_sw_blending = obj->blending;
      SetupBSplineMeshWidgets ( obj );
      if ( bsm_sw_blending )
        side10wdg_bsm_opt = side10wdg_bsm_opt3;
      else
        side10wdg_bsm_opt = side10wdg_bsm_opt1;
      xge_SetMenuWidgets ( side10menu, side10wdg_bsm_opt, false );
      Side10MenuBsmResize ( SIDEMENUWIDTH0, xge_current_height-TOPMENUHEIGHT );
      if ( ChangeSide10MenuWidth ( side00menu->h ) )
        ResizeWindow1 ( xge_current_width, xge_current_height );
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_G1:
      if ( sw_hfill_g1 ) sw_hfill_g2 = sw_hfill_g1q2 = false;
      else sw_hfill_g2 = true, sw_hfill_g1q2 = false;
      if ( obj->degree == 3 ) {
        GeomObjectBSplineMeshSetHoleFillingOpt ( obj, sw_hfill_coons, sw_hfill_bezier,
                sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2 );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_G2:
      if ( sw_hfill_g2 ) sw_hfill_g1 = sw_hfill_g1q2 = false;
      else sw_hfill_g1q2 = true, sw_hfill_g1 = false;
      if ( obj->degree == 3 ) {
        GeomObjectBSplineMeshSetHoleFillingOpt ( obj, sw_hfill_coons, sw_hfill_bezier,
                sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2 );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_G1_QUASI_G2:
      if ( sw_hfill_g1q2 ) sw_hfill_g1 = sw_hfill_g2 = false;
      else sw_hfill_g1 = true, sw_hfill_g2 = false;
      if ( obj->degree == 3 ) {
        GeomObjectBSplineMeshSetHoleFillingOpt ( obj, sw_hfill_coons, sw_hfill_bezier,
                sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2 );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_HOLEFILL_COONS:
      sw_hfill_bezier = !sw_hfill_coons;
      if ( obj->degree == 3 ) {
        GeomObjectBSplineMeshSetHoleFillingOpt ( obj, sw_hfill_coons, sw_hfill_bezier,
                sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2 );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_HOLEFILL_BEZIER:
      sw_hfill_coons = !sw_hfill_bezier;
      if ( obj->degree == 3 ) {
        GeomObjectBSplineMeshSetHoleFillingOpt ( obj, sw_hfill_coons, sw_hfill_bezier,
                sw_hfill_g1, sw_hfill_g2, sw_hfill_g1q2 );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else {
        xge_SetClipping ( side10menu );
        side10menu->redraw ( side10menu, true );
      }
      return 1;
  case swM1BSM_SHAPE_ONLY:
      obj->bl_shape_only = bsm_sw_shape_only;
      SuggestOptLevels ( obj, obj->nlevels > 0 );
      return 1;
  case swM1BSM_BLENDING_CONSTRAINTS:
      obj->bl_constr = bsm_sw_constr;
      return 1;
  case swM1BSM_ALT_MULTILEVEL:
      return 1;
  case swM1BSM_COARSE_PRECOND:
      obj->bl_use_coarse = bsm_sw_use_coarse;
      SuggestOptLevels ( obj, obj->nlevels > 0 );
      return 1;
  default:
      return 0;
    }

case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BSM_REFINEMENT:
      if ( GeomObjectBSplineMeshRefinement ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        xge_SetWindow ( win0 );
        rendered_picture = false;
        xge_Redraw ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      return 1;
  case btnM1BSM_DOUBLING:
      if ( GeomObjectBSplineMeshDoubling ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        xge_SetWindow ( win0 );
        rendered_picture = false;
        xge_Redraw ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      return 1;
  case btnM1BSM_AVERAGING:
      if ( GeomObjectBSplineMeshAveraging ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        xge_SetWindow ( win0 );
        rendered_picture = false;
        xge_Redraw ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      return 1;
  case btnM1BSM_EXTRACTSUBMESH:
      if ( GeomObjectBSplineMeshExtractSubmesh ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj ); 
        xge_SetWindow ( win0 );
        rendered_picture = false;
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCouldNotExtractMesh, 0 );
      return 1;
  case btnM1BSM_MARK_VERT:
      if ( GeomObjectBSplineMeshMarkBetweenVertices ( obj, 1 ) ) {
        xge_SetWindow ( win0 );
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_UNMARK_VERT:
      if ( GeomObjectBSplineMeshMarkBetweenVertices ( obj, 0 ) ) {
        xge_SetWindow ( win0 );
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_MARK_HEDGE:
      if ( GeomObjectBSplineMeshMarkHalfedgesBetweenVertices ( obj, 1 ) ) {
        xge_SetWindow ( win0 );
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_UNMARK_HEDGE:
      if ( GeomObjectBSplineMeshMarkHalfedgesBetweenVertices ( obj, 0 ) ) {
        xge_SetWindow ( win0 );
        xge_Redraw ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_FILTER:
      if ( GeomObjectBSplineMeshFilterPolyline ( obj ) ) {
        xge_SetClipping ( geom00win );
        geom00win->redraw ( geom00win, true );
      }
      return 1;
  case btnM1BSM_ENTER_LINE:
      if ( GeomObjectBSplineMeshGetCurrentVertices ( obj, 2, p ) ) {
        EnterRefLine ( p );
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_REMOVE_VERTEX:
      if ( GeomObjectBSplineMeshRemoveCurrentVertex ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveVertex, 0 );
      return 1;
  case btnM1BSM_DIVIDE_FACET:
      if ( GeomObjectBSplineMeshDivideFacet ( obj ) ) {
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgNoValidPairOfVertices, 0 );
      return 1;
  case btnM1BSM_DOUBLE_LOOP:
      if ( GeomObjectBSplineMeshDoubleEdgeLoop ( obj ) ) {
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoubleLoop, 0 );
      return 1;
  case btnM1BSM_SHRINK_EDGE:
      if ( GeomObjectBSplineMeshShrinkCurrentEdge ( obj ) ) {
        rendered_picture = false;
        xge_SetClipping ( geom00win );
        geom00win->redraw ( geom00win, true );
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoIt, 0 );
      return 1;
  case btnM1BSM_CONTRACT_EDGE:
      if ( GeomObjectBSplineMeshContractCurrentEdge ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotContractEdge, 0 );
      return 1;
  case btnM1BSM_GLUE_EDGES:
      if ( GeomObjectBSplineMeshGlueEdges ( obj ) ) {
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotGlueEdges, 0 );
      return 1;
  case btnM1BSM_GLUE_EDGE_LOOPS:
      if ( GeomObjectBSplineMeshGlueEdgeLoops ( obj ) ) {
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotGlueLoops, 0 );
      return 1;
  case btnM1BSM_SEAL_HOLE:
      if ( GeomObjectBSplineMeshSealHole ( obj ) ) {
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotSealHole, 0 );
      return 1;
  case btnM1BSM_SPLIT_BOUNDARY_EDGE:
      if ( GeomObjectBSplineMeshSplitBoundaryEdge ( obj ) ) {
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMsgMustBeBoundaryEdge, 0 );
      return 1;
  case btnM1BSM_REMOVE_FACET:
      if ( GeomObjectBSplineMeshRemoveCurrentFacet ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveFacet, 0 );
      return 1;
  case btnM1BSM_DOUBLE_FAC_EDGES:
      if ( GeomObjectBSplineMeshDoubleCurrentFacEdges ( obj ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      else if ( !obj->integrity_ok )
        xge_DisplayErrorMessage ( ErrorMsgMeshIntegrity, 0 );
      else
        xge_DisplayErrorMessage ( ErrorMsgCannotDoubleFacetEdges, 0 );
      return 1;
  case btnM1BSM_DATA_KGON:
      if ( GeomObjectBSplineMeshInitKGon ( obj, bsm_degree1, bsm_sw_data_add ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSM_DATA_TETRAHEDRON:
      if ( GeomObjectBSplineMeshInitTetrahedron ( obj, bsm_sw_data_add ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSM_DATA_CUBE:
      if ( GeomObjectBSplineMeshInitCube ( obj, bsm_sw_data_add ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSM_DATA_DODECAHEDRON:
      if ( GeomObjectBSplineMeshInitDodecahedron ( obj, bsm_sw_data_add ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSM_DATA_KPRISM:
      if ( GeomObjectBSplineMeshInitKPrism ( obj, bsm_degree2, bsm_sw_data_add ) ) {
        bsm_sw_blending = obj->blending;
        SetupBSplineMeshVEFnum ( obj );
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSM_PRETRANSFORMATION:
      GetPreTransformation ();
      OpenPopup ( popup14, false );
      obj->me.display_pretrans = true;
      editing_pretrans = true;
      rendered_picture = false;
      xge_SetWindow ( win0 );
      xge_Redraw ();
      return 1;
  case btnM1BSM_BLENDING_OPTIMIZE:
      if ( ipc_state == ipcstate_CHILD_BUSY ) {
        if ( obj->me.bound_with_a_child ) {
          IPCInterruptTheChild ();
          bsm_bl_optimize->data0 = txtOptimize;
        }
        else {
          xge_DisplayErrorMessage ( ErrorMsgChildProcessBusy, 0 );
          return 1;
        }
      }
      else {
        if ( BlendingMeshOptimizationPrepareData ( obj ) ) {
          InitBlendingMeshOptimization ();
          bsm_bl_optimize->data0 = txtInterrupt;
        }
      }
      xge_SetClipping ( er );
      er->redraw ( er, true );
      return 1;
  case btnM1SUGGESTLEVELS:
      SuggestOptLevels ( obj, true );
      return 1;
  case btnM1SUGGESTBLOCKS:
      obj->bl_shape_only = bsm_sw_shape_only = false;
      SuggestOptLevels ( obj, false );
      return 1;
  case btnM1BSM_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  case btnM1BSM_OPTIMIZE_SPECIALS:
      xge_GrabFocus ( xge_null_widget, true );
      xge_SetWindowCursor ( win0, xgeCURSOR_WATCH );
      xge_SetWindowCursor ( win1, xgeCURSOR_WATCH );
      xge_PostIdleCommand ( IDLE_COMMAND_BSMESH_OPT_SPECIALS, 0, 0 );
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    switch ( er->id ) {
  case slM1BSM_G1Q2_PARAM:
      GeomObjectBSplineMeshSetG1Q2param ( obj, sl_g1q2_param );
      NotifyParam1 ( obj->g1q2param );
      if ( obj->degree == 3 && obj->fill_G1Q2 ) {
        rendered_picture = false;
        RedrawGeom00Win ();
      }
      return 1;
  case slM1BSM_BLENDING_CPARAM:
      obj->bsm_bl_C = xge_LogSlidebarValued ( BSM_BL_MIN_CPARAM, BSM_BL_MAX_CPARAM,
                                              sl_bsm_bl_param );
      NotifyParam2 ( obj->bsm_bl_C );
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BSM_DEG:
      if ( GeomObjectBSplineMeshSetDegree ( obj, key ) ) {
        bsm_sw_blending = obj->blending;
        degree = obj->degree;
        rendered_picture = false;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_VERTEX0:
      if ( GeomObjectBSplineMeshSetCurrentVertex ( obj, key, 0 ) ) {
        bsm_vertex_num0 = obj->current_vertex[0];
        BottomDisplayPoint ( win1, obj->me.spdimen, obj->me.cpdimen, bsm_vertex_num0,
                             &obj->meshvpc[obj->me.cpdimen*bsm_vertex_num0], false );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_VERTEX1:
      if ( GeomObjectBSplineMeshSetCurrentVertex ( obj, key, 1 ) ) {
        bsm_vertex_num1 = obj->current_vertex[1];
        BottomDisplayPoint ( win1, obj->me.spdimen, obj->me.cpdimen, bsm_vertex_num1,
                             &obj->meshvpc[obj->me.cpdimen*bsm_vertex_num1], false );
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_EDGE0:
      if ( GeomObjectBSplineMeshSetCurrentEdge ( obj, key, 0 ) ) {
        bsm_edge_num0 = obj->current_edge[0];
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_EDGE1:
      if ( GeomObjectBSplineMeshSetCurrentEdge ( obj, key, 1 ) ) {
        bsm_edge_num1 = obj->current_edge[1];
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_FACET0:
      if ( GeomObjectBSplineMeshSetCurrentFacet ( obj, key, 0 ) ) {
        bsm_facet_num0 = obj->current_facet[0];
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_FACET1:
      if ( GeomObjectBSplineMeshSetCurrentFacet ( obj, key, 1 ) ) {
        bsm_facet_num1 = obj->current_facet[1];
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_DATA_KGON:
      bsm_degree1 = key;
      xge_SetClipping ( er );
      er->redraw ( er, true );
      return 1;
  case intwM1BSM_DATA_KPRISM:
      bsm_degree2 = key;
      xge_SetClipping ( er );
      er->redraw ( er, true );
      return 1;
  case intwM1BSM_DENSITY:
      if ( GeomObjectBSplineMeshSetDensLevel ( obj, key ) ) {
        if ( obj->subdivision )
          density_u = obj->subdivl;
        else
          density_u = obj->density;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BSM_NKN1:
      obj->nkn1 = bsm_bl_nkn1 = key;
      if ( bsm_bl_nkn2 < key ) {
        obj->nkn2 = bsm_bl_nkn2 = key;
        xge_SetClipping ( bsm_bl_qknots2.er );
        bsm_bl_qknots2.er->redraw ( bsm_bl_qknots2.er, true );
      }
      return 1;
  case intwM1BSM_NKN2:
      obj->nkn2 = bsm_bl_nkn2 = key;
      if ( bsm_bl_nkn1 > key ) {
        obj->nkn1 = bsm_bl_nkn1 = key;
        xge_SetClipping ( bsm_bl_qknots1.er );
        bsm_bl_qknots1.er->redraw ( bsm_bl_qknots1.er, true );
      }
      return 1;
  case intwM1BSM_NLEVELS:
      obj->nlevels = bsm_bl_nlev = key;
      bsm_bl_nbl = (0x0001 << bsm_bl_nlev) - 1;
      if ( bsm_bl_nbl < 1 )
        bsm_bl_nbl = 1;
      else
        bsm_bl_nblocks.maxvalue = max ( bsm_bl_nblocks.maxvalue, bsm_bl_nbl );
      obj->nblocks = bsm_bl_nbl;
      if ( obj->nlevels ) {
        obj->startfrom = obj->nblocks-1;
        SetupStartFrom ( obj );
      }
      else
        SetupMaxIter ( obj );
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case intwM1BSM_NBLOCKS:
      obj->nlevels = bsm_bl_nlev = 0;
      obj->bl_shape_only = false;
      bsm_bl_nlevels.maxvalue = G2MBL_MAX_LEVELS;
      bsm_bl_nblocks.maxvalue = G2MBL_MAX_BLOCKS;
      key = min ( key, G2MBL_MAX_BLOCKS );
      obj->nblocks = bsm_bl_nbl = key;
      SetupMaxIter ( obj );
      xge_SetClipping ( side10menu );
      side10menu->redraw ( side10menu, true );
      return 1;
  case intwM1BSM_MAXIT:
      if ( obj->nlevels ) {
        key = max ( 0, key );
        key = min ( obj->nblocks-1, key );
        obj->startfrom = bsm_bl_maxit = key;
      }
      else {
        key = max ( 1, key );
        obj->maxit = bsm_bl_maxit = key;
      }
      return 1;
  case intwM1BSM_NPTHREADS:
      bsm_npthreads = key;
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BSM_NAME:
      return 1;
  case textedM1BSM_COARSE_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BSM_NAME:
      return 1;
  case textedM1BSM_COARSE_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_RESIZE:
    Side10MenuBsmResize ( x, y );
    return 1;

default:
    return 0;
  }
} /*Side10MenuBsmCallBack*/

