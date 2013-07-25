
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

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE


void InitSide10Menu_BSp ( void )
{
  xge_widget *w;

        /* widgets specific for B-spline patches */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBSplinePatch );
  w = xge_NewStringEd ( win1, w, textedM1BSP_NAME, 109, 19, 0, 40,
                        MAX_NAME_LENGTH, objectname, &bsp_name_ed );
  w = bsp_type_button =
           xge_NewButton ( win1, w, btnM1BSP_TYPE, 76, 19, 0, 60, txtGeneral );
  if ( w ) w->state = xgestate_BUTTON_COMBO_0;
  w = xge_NewIntWidget ( win1, w, intwM1BSP_DEGU, 76, 19, 0, 80,
                         0, MAX_DEGREE+1, &intw_bspdegu, txtDegreeU, &degreeu );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_DEGV, 76, 19, 0, 100,
                         0, MAX_DEGREE+1, &intw_bspdegv, txtDegreeV, &degreev );
  w = xge_NewButton ( win1, w, btnM1BSP_FLIP, 76, 19, 0, 120, txtFlip );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, swM1BSP_VIEW_SURF, 109, 16, 0, 22,
                      txtSurface, &sw_view_surf );
  w = xge_NewSwitch ( win1, w, swM1BSP_VIEW_CNET, 109, 16, 0, 42,
                      txtControlNet, &sw_view_cnet );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_DENSITY_U, 109, 19, 0, 62,
                         1, MAX_PNET_DENSITY, &intwdensityu, txtDensityU,
                         &density_u );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_DENSITY_V, 109, 19, 0, 82,
                         1, MAX_PNET_DENSITY, &intwdensityv, txtDensityV,
                         &density_v );
  w = xge_NewButton ( win1, w, btnM1BSP_COLOUR, 60, 18, 0, 102, txtColour );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_view = w;
          /* data */
  w = xge_NewSwitch ( win1, NULL, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_data = w;
          /* options */
            /* general */
  w = xge_NewButton ( win1, NULL, btnM1BSP_UNIFORM_U, 76, 19, 0, 20,
                      txtUniformU );
  w = xge_NewButton ( win1, w, btnM1BSP_UNIFORM_V, 76, 19, 0, 40,
                      txtUniformV );
  w = xge_NewSwitch ( win1, w, swM1BSP_CLOSED_U, 109, 16, 0, 60,
                      txtClosedU, &bsp_sw_closed_u );
  w = xge_NewSwitch ( win1, w, swM1BSP_CLOSED_V, 109, 16, 0, 80,
                      txtClosedV, &bsp_sw_closed_v );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_opt = side10wdg_bsp_opt_general = w;
            /* spherical product */
  w = xge_NewSwitch ( win1, NULL, swM1BSP_SPRODUCT_EQUATOR, 109, 16, 0, 22,
                      txtEquator, &sw_bsp_sproduct_equator );
  w = xge_NewStringEd ( win1, w, textedM1BSP_SPRODUCT_EQNAME, 109, 19, 0, 42,
                        MAX_NAME_LENGTH, objectname, &bsp_sproduct_eqname_ed );
  w = xge_NewIntWidget ( win1, w, intWM1BSP_SPRODUCT_EQDEG, 76, 19, 0, 62,
                         0, MAX_DEGREE+1, &intw_bsp_sproduct_eqdeg, txtDegree,
                         &bsp_eqdegree );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_EQRATIONAL, 109, 16, 0, 82,
                      txtRational, &sw_bsp_sproduct_eqrational );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_EQCLOSED, 109, 16, 0, 102,
                      txtClosed, &sw_bsp_sproduct_eqclosed );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_EQUNIFORM, 109, 16, 0, 122,
                      txtUniform, &sw_bsp_sproduct_equniform );

  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_MERIDIAN, 109, 16, 0, 152,
                      txtMeridian, &sw_bsp_sproduct_meridian );
  w = xge_NewStringEd ( win1, w, textedM1BSP_SPRODUCT_MERNAME, 109, 19, 0, 172,
                        MAX_NAME_LENGTH, objectname, &bsp_sproduct_mername_ed );
  w = xge_NewIntWidget ( win1, w, intWM1BSP_SPRODUCT_MERDEG, 76, 19, 0, 192,
                         0, MAX_DEGREE+1, &intw_bsp_sproduct_merdeg, txtDegree,
                         &bsp_merdegree );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_MERRATIONAL, 109, 16, 0, 212,
                      txtRational, &sw_bsp_sproduct_merrational );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_MERCLOSED, 109, 16, 0, 232,
                      txtClosed, &sw_bsp_sproduct_merclosed );
  w = xge_NewSwitch ( win1, w, swM1BSP_SPRODUCT_MERUNIFORM, 109, 16, 0, 252,
                      txtUniform, &sw_bsp_sproduct_meruniform );

  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_opt_spherical = w;
            /* blending G1 */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtG1BlendingPatch );
  w = xge_NewSwitch ( win1, w, swM1BSP_BLENDING_CLAMPED, 109, 16, 0, 42,
                      txtClamped, &sw_bsp_clamped );
  w = xge_NewSwitch ( win1, w, swM1BSP_CLOSED_U, 109, 16, 0, 62,
                      txtClosedU, &bsp_sw_closed_u );
  w = xge_NewButton ( win1, w, btnM1BSP_BLENDING_REFINE, 58, 19, 0, 82,
                      txtRefine );
  w = xge_NewSwitch ( win1, w, swM1BSP_BLENDING_ENTIRE, 16, 16, 29, 125,
                      NULL, &bsp_blending_entire );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G1BLENDING_UMIN, 23, 19, 0, 124,
                         3, MAX_BSP_KNOTS, &bsp_blG1_minu,
                         NULL, &bsp_bl_opt_range[0] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G1BLENDING_UMAX, 23, 19, 50, 124,
                         3, MAX_BSP_KNOTS, &bsp_blG1_maxu,
                         NULL, &bsp_bl_opt_range[1] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G1BLENDING_VMIN, 23, 19, 25, 144,
                         3, MAX_BSP_KNOTS, &bsp_blG1_minv,
                         NULL, &bsp_bl_opt_range[2] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G1BLENDING_VMAX, 23, 19, 25, 104,
                         3, MAX_BSP_KNOTS, &bsp_blG1_maxv,
                         NULL, &bsp_bl_opt_range[3] );
  w = xge_NewSwitch ( win1, w, swM1BSP_NHARMONIC, 109, 16, 0, 166,
                      txtBiharmonic, &bsp_bl_nharmonic );
  w = xge_NewButton ( win1, w, btnM1BSP_PRETRANSFORMATION, 58, 19, 0, 186,
                      txtPreTransf );
  w = xge_NewSlidebard ( win1, w, slM1BSP_BLENDING_CPARAM, 109, 10, 0, 206,
                         &sl_bsm_bl_param );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_NKN1, 91, 19, 0, 220,
                         3, 10, &bsm_bl_qknots1, txtQKnots, &bsm_bl_nkn1 );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_NKN2, 17, 19, 92, 220,
                         3, 10, &bsm_bl_qknots2, txtNull, &bsm_bl_nkn2 );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_MAXIT, 91, 19, 0, 240,
                         1, 50, &bsm_bl_maxiter, txtMaxIter, &bsm_bl_maxit );
  w = bsp_bl_optimizeG1 =
      xge_NewButton ( win1, w, btnM1BSP_BLENDING_OPTIMIZE, 58, 19, 0, 260,
                      txtOptimize );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_opt_blendingG1 = w;
            /* blending G2 */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtG2BlendingPatch );
  w = xge_NewSwitch ( win1, w, swM1BSP_BLENDING_CLAMPED, 109, 16, 0, 42,
                      txtClamped, &sw_bsp_clamped );
  w = xge_NewSwitch ( win1, w, swM1BSP_CLOSED_U, 109, 16, 0, 62,
                      txtClosedU, &bsp_sw_closed_u );
  w = xge_NewButton ( win1, w, btnM1BSP_BLENDING_REFINE, 58, 19, 0, 82,
                      txtRefine );
  w = xge_NewSwitch ( win1, w, swM1BSP_BLENDING_ENTIRE, 16, 16, 29, 125,
                      NULL, &bsp_blending_entire );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G2BLENDING_UMIN, 23, 19, 0, 124,
                         3, MAX_BSP_KNOTS, &bsp_blG2_minu,
                         NULL, &bsp_bl_opt_range[0] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G2BLENDING_UMAX, 23, 19, 50, 124,
                         3, MAX_BSP_KNOTS, &bsp_blG2_maxu,
                         NULL, &bsp_bl_opt_range[1] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G2BLENDING_VMIN, 23, 19, 25, 144,
                         3, MAX_BSP_KNOTS, &bsp_blG2_minv,
                         NULL, &bsp_bl_opt_range[2] );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_G2BLENDING_VMAX, 23, 19, 25, 104,
                         3, MAX_BSP_KNOTS, &bsp_blG2_maxv,
                         NULL, &bsp_bl_opt_range[3] );
  w = xge_NewSwitch ( win1, w, swM1BSP_NHARMONIC, 109, 16, 0, 166,
                      txtTriharmonic, &bsp_bl_nharmonic );
  w = xge_NewButton ( win1, w, btnM1BSP_PRETRANSFORMATION, 58, 19, 0, 186,
                      txtPreTransf );
  w = xge_NewSlidebard ( win1, w, slM1BSP_BLENDING_CPARAM, 109, 10, 0, 206,
                         &sl_bsm_bl_param );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_NKN1, 91, 19, 0, 220,
                         3, 10, &bsm_bl_qknots1, txtQKnots, &bsm_bl_nkn1 );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_NKN2, 17, 19, 92, 220,
                         3, 10, &bsm_bl_qknots2, txtNull, &bsm_bl_nkn2 );
  w = xge_NewIntWidget ( win1, w, intwM1BSP_MAXIT, 91, 19, 0, 240,
                         1, 50, &bsm_bl_maxiter, txtMaxIter, &bsm_bl_maxit );
  w = bsp_bl_optimizeG2 =
      xge_NewButton ( win1, w, btnM1BSP_BLENDING_OPTIMIZE, 58, 19, 0, 260,
                      txtOptimize );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_COORD, 109, 16, 0, xge_HEIGHT-56,
                      txtCoordinates, &sw_bsp_dom_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -56 );
  w = xge_NewSwitch ( win1, w, swM1BSP_DOMAIN_PANZOOM, 109, 16, 0, xge_HEIGHT-36,
                      txtPanZoom, &sw_bsp_dom_panzoom );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsp_opt_blendingG2 = w;
} /*InitSide10Menu_BSp*/

