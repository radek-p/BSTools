
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "ed3dswidgets.h"
#include "pomnijipc.h"

void UpdateBlendingRangeWidgets ( void )
{
  int i;

  if ( kwind.closed_u ) {
    blending_opt_range[0].minvalue = blending_opt_range[1].minvalue = -1;
    blending_opt_range[0].maxvalue = blending_opt_range[1].maxvalue =
      lastknot_u-2*degree_u;
  }
  else {
    blending_opt_range[0].minvalue = blending_opt_range[1].minvalue = degree_u;
    blending_opt_range[0].maxvalue = blending_opt_range[1].maxvalue =
      lastknot_u-2*degree_u-1;
  }
  blending_opt_part[0] = min ( blending_opt_part[0], blending_opt_range[0].maxvalue );
  blending_opt_part[0] = max ( blending_opt_part[0], blending_opt_range[0].minvalue );
  blending_opt_part[1] = min ( blending_opt_part[1], blending_opt_range[1].maxvalue );
  blending_opt_part[1] = max ( blending_opt_part[1], blending_opt_range[1].minvalue );
  if ( kwind.closed_u ) {
        /* for closed blending surfaces wrapping is possible */
    if ( blending_opt_part[0] < 0 )
      blending_opt_part[0] += lastknot_u-2*degree_u;
    else if ( blending_opt_part[0] >= lastknot_u-2*degree_u )
      blending_opt_part[0] -= lastknot_u-2*degree_u;
    if ( blending_opt_part[1] < 0 )
      blending_opt_part[1] += lastknot_u-2*degree_u;
    else if ( blending_opt_part[1] >= lastknot_u-2*degree_u )
      blending_opt_part[1] -= lastknot_u-2*degree_u;
  }
  else {
    if ( blending_opt_part[0] > blending_opt_part[1] ) {
      i = blending_opt_part[0];
      blending_opt_part[0] = blending_opt_part[1];
      blending_opt_part[1] = i;
    }
  }

  blending_opt_range[2].minvalue = blending_opt_range[3].minvalue = degree_v;
  blending_opt_range[2].maxvalue = blending_opt_range[3].maxvalue =
    lastknot_v-2*degree_v-1;
  blending_opt_part[2] = min ( blending_opt_part[2], blending_opt_range[2].maxvalue );
  blending_opt_part[2] = max ( blending_opt_part[2], blending_opt_range[2].minvalue );
  blending_opt_part[3] = min ( blending_opt_part[3], blending_opt_range[3].maxvalue );
  blending_opt_part[3] = max ( blending_opt_part[3], blending_opt_range[3].minvalue );
  if ( blending_opt_part[2] > blending_opt_part[3] ) {
    i = blending_opt_part[2];
    blending_opt_part[2] = blending_opt_part[3];
    blending_opt_part[3] = i;
  }
  if ( kwind.closed_u ) {
    sw_blending_opt_entire =
        (blending_opt_part[0] == blending_opt_range[0].minvalue+1) &&
        (blending_opt_part[1] == blending_opt_range[1].maxvalue-1) &&
        (blending_opt_part[2] == blending_opt_range[2].minvalue) &&
        (blending_opt_part[3] == blending_opt_range[3].maxvalue);
  }
  else {
    sw_blending_opt_entire =
        (blending_opt_part[0] == blending_opt_range[0].minvalue) &&
        (blending_opt_part[1] == blending_opt_range[1].maxvalue) &&
        (blending_opt_part[2] == blending_opt_range[2].minvalue) &&
        (blending_opt_part[3] == blending_opt_range[3].maxvalue);
  }
  blending_mat_valid = false;
} /*UpdateBlendingRangeWidgets*/

boolean IsValidDoubleString ( char *s )
{
  double x;

  return sscanf ( s, "%lf", &x ) == 1;
} /*IsValidDoubleString*/

int Win1CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  boolean      dd;
  xge_2Dwind   *cwind;
  xge_KnotWind *ckwind;

  cwind = eq_cwind.er->data0;
  ckwind = eq_ckwind.er->data0;
  if ( er ) {
    switch ( msg ) {
  case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
    case btn14aTYPE:
    case btn14bTYPE:
        xge_SetWindow ( win1 );
        xge_AddPopup ( popup10 );
        xge_GrabFocus ( popup10, true );
        break;
    case btn14bEDIT:
        xge_SetMenuWidgets ( menu6, menu60list, true );
        break;
    case btn14bVIEW:
        xge_SetMenuWidgets ( menu6, menu61list, true );
        break;
    case btn14bARCS:
        xge_SetMenuWidgets ( menu6, menu62list, true );
        break;
    case btn16GENERAL:  /* General */
        if ( win1_contents != WIN1_GENERAL ) {
          xge_RemovePopup ( true );
          sw_bind_blending = sw_triharmonic_blending = sw_nonlin_blending = false;
          BreakNLBlending ();
          sw_blending_constraints = false;
          xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
          win1_contents = WIN1_GENERAL;
          domwind->next = menu3;
          menu3->prev = domwind;
          xge_SetWindow ( win1 );
          xge_SetWinEdRect ( menu7 );
          xge_RedrawAll ();
        }
        break;
    case btn16SPHERICAL:  /* Spherical */
        if ( win1_contents != WIN1_SPHERICAL ) {
          xge_RemovePopup ( true );
          sw_bind_blending = sw_triharmonic_blending = sw_nonlin_blending = false;
          BreakNLBlending ();
          sw_blending_constraints = false;
          xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
          win1_contents = WIN1_SPHERICAL;
          if ( equator )
            SelectEquator ();
          else 
            SelectMeridian ();
          ProjectEqMerCurve ();
          xge_SetWindow ( win1 );
          xge_SetWinEdRect ( menu8 );
          xge_RedrawAll ();
        }
        break;
    case btn16SWEPT:     /* Swept */
        if ( win1_contents != WIN1_SWEPT ) {
          xge_RemovePopup ( true );
          sw_bind_blending = sw_triharmonic_blending = sw_nonlin_blending = false;
          BreakNLBlending ();
          sw_blending_constraints = false;
          xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
          win1_contents = WIN1_SWEPT;
          domwind->next = menu9;
          menu9->prev = domwind;
          xge_SetWindow ( win1 );
          xge_SetWinEdRect ( menu13 );
          xge_RedrawAll ();
        }
        break;
    case btn16BLENDING:  /* Blending */
        if ( win1_contents != WIN1_BLENDING ) {
          xge_RemovePopup ( true );
          win1_contents = WIN1_BLENDING;
          domwind->next = menu11;
          menu11->prev = domwind;
          UpdateBlendingRangeWidgets ();
          if ( degree_u == 2 && degree_v == 2 ) {
            sw_blending_g1 = true;
            sw_blending_g2 = false;
            bl_trihsw->data0 = txtBiharmonic;
          }
          else if ( degree_u == 3 && degree_v == 3 ) {
            sw_blending_g1 = false;
            sw_blending_g2 = true;
            bl_trihsw->data0 = txtTriharmonic;
          }
          xge_SetWindow ( win1 );
          xge_SetWinEdRect ( menu14 );
          xge_RedrawAll ();
        }
        break;
    case btn13EQUIDIST_U:
        SetEquidistantU ();
        goto redraw_equidist;
    case btn13EQUIDIST_V:
        SetEquidistantV ();
redraw_equidist:
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        xge_RedrawAll ();
        break;
    case btn13FLIP:
        bind_spr = sw_bind_blending = sw_triharmonic_blending =
          sw_nonlin_blending = false;
        FlipPatch ();
        xge_Redraw ();
        break;
    case btn15aRESET:
        ResetEqMer ();
        if ( bind_spr )
          xge_RedrawAll ();
        else
          xge_Redraw ();
        break;
    case btn15cQUARTER_CIRCLE:
        EqMerSetQuarterCircleArc ();
        goto bind_spr_update2;
    case btn15cHALF_CIRCLE:
        EqMerSetHalfCircleArc ();
        goto bind_spr_update2;
    case btn15cFULL_CIRCLE:
        EqMerSetFullCircleArc ();
        goto bind_spr_update2;
    case btn112INIT_BLENDING:
        sw_blending_constraints = false;
        xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
        if ( sw_blending_g1 )
          InitG1BlendingSurface ();
        else
          InitG2BlendingSurface ();
        SetKWindNKnots ();
        UpdateBlendingRangeWidgets ();
        xge_SetWindow ( win0 );
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        xge_RedrawAll ();
        break;
    case btn112BLENDING_REFINE:
        bind_spr = sw_bind_blending = sw_triharmonic_blending =
          sw_nonlin_blending = false;
        if ( sw_blending_g1 )
          dd = RefineG1BlendingSurface ();
        else
          dd = RefineG2BlendingSurface ();
        if ( dd ) {
          UpdateBlendingRangeWidgets ();
          xge_SetWindow ( win0 );
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          UpdateBlendingRangeWidgets ();
          xge_RedrawAll ();
        }
        else
          xge_DisplayErrorMessage ( ErrorMsgRefineBlending, -1 );
        break;
    case btn112BLENDING_LMT_ITER:
        bind_spr = sw_triharmonic_blending = false;
/*
        if ( !sw_blending_g2 ) {
          xge_DisplayErrorMessage ( ErrorMsgNotImplemented, -1 );
          return 1;
        }
*/
        if ( !xge_ChildIsActive () ) {
          if ( !MakeTheChildProcess () ) {
            xge_DisplayErrorMessage ( ErrorChildProcessNotActive, -1 );
            return 1;
          }
        }
        switch ( ipc_state ) {
      case ipcSTATE_NOTHING:  /* initiate sending data to the child process */
          xge_CallTheChild ( cmdGET_OPTIONS, sizeof(ipc_options) );
          xge_PostIdleCommand ( IDLE_COMMAND_SEND_OPTIONS, 0, 0 );
          ipc_state = ipcSTATE_OPTIONS_SENT;
          break;
      case ipcSTATE_G2BLOPT_LAUNCHED: /* interrupt the blending surface optimization */
          BreakNLBlending ();
          break;
      default:
          break;
        }
        xge_SetClipping ( bl_optbtn );
        bl_optbtn->redraw ( bl_optbtn, true );
        break;
    case btn111INFO:
        xge_DisplayErrorMessage ( ErrorMsgNotImplemented, -1 );
        break;
    case btn111TRANSFORM:
        xge_AddPopup ( popup11 );
        display_trans_net = true;
        xge_SetWindow ( win0 );
        xge_Redraw ();
        xge_SetWindow ( win1 );
        xge_GrabFocus ( popup11, true );
        break;
    case btnP11IDENTITY:
        ResetBlendingOptTrans ();
        xge_SetWindow ( win0 );
        xge_Redraw ();
        xge_SetWindow ( win1 );
        break;
    case btnP11OK:
        xge_RemovePopup ( true );
        display_trans_net = false;
        xge_SetWindow ( win0 );
        xge_Redraw ();
        xge_SetWindow ( win1 );
        break;
    case btn111OPTIONS:
        xge_AddPopup ( popup12 );
        xge_GrabFocus ( popup12, true );
        break;
    case btnP12OK:
        xge_RemovePopup ( true );
        break;
    default:
        break;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
    case sw13CLOSED_U:
        if ( SetClosedUSurface () )
          RedrawGeomWindows ( true, true );
        break;
    case sw13CLOSED_V:
        if ( SetClosedVSurface () )
          RedrawGeomWindows ( true, true );
        break;
    case sw13DOMAIN_NET:  /* display domain net */
        xge_SetClipping ( domwind );
        domwind->redraw ( domwind, true );
        break;
    case sw13MARK_UNMARK:  /* mark/unmark switch */
        if ( kwind.selecting_mode ) {
          dd = display_domain_net;
          if ( kwind.panning || !dd ) {
            display_domain_net = true;
            kwind.panning = false;
            xge_SetClipping ( menu2 );
            menu3->redraw ( menu2, true );
          }
          if ( !dd )
            RedrawGeomWindows ( false, true );
        }
        break;
    case sw13PAN_ZOOM:  /* pan & zoom */
        if ( kwind.selecting_mode ) {
          kwind.selecting_mode = false;
          xge_SetClipping ( menu2 );
          menu3->redraw ( menu2, true );
        }
        break;
    case sw15aBIND:
        if ( bind_spr ) {
          BindSphericalProduct ();
          SetKWindNKnots ();
          xge_SetWindow ( win0 );
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          xge_RedrawAll ();
        }
        break;
    case sw15aEQUATOR:  /* Equator */
        if ( equator )
          SelectEquator ();
        else
          SelectMeridian ();
        goto switch_spr;;
    case sw15aMERIDIAN:  /* Meridian */
        if ( meridian )
          SelectMeridian ();
        else
          SelectEquator ();
switch_spr:


        xge_Redraw ();
        break;
    case sw15aCLOSED:
        SetClosedEqMer ();
        if ( bind_spr ) {
          BindSphericalProduct ();
          SetKWindNKnots ();
          xge_RedrawAll ();
        }
        else
          xge_Redraw ();
        break;
    case sw15aNURBS:
        SetupNURBSEqMer ();
        if ( bind_spr ) {
          BindSphericalProduct ();
          SetKWindNKnots ();
          xge_RedrawAll ();
        }
        else
          xge_Redraw ();
        break;
    case sw15aMARK_UNMARK:
        if ( (mer_cwind.selecting_mode = eq_cwind.selecting_mode) ) {
          eq_ckwind.panning = eq_cwind.panning = eq_cwind.moving_tool =
          eq_cwind.scaling_tool = eq_cwind.rotating_tool =
          eq_cwind.shear_tool = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_SELECTING_TOOL );
          mer_ckwind.panning = mer_cwind.panning = mer_cwind.moving_tool =
          mer_cwind.scaling_tool = mer_cwind.rotating_tool =
          mer_cwind.shear_tool = false;
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_SELECTING_TOOL );
        }
        else {
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_NO_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_NO_TOOL );
        }
        xge_Redraw ();
        break;
    case sw15aMOVE:
        if ( (mer_cwind.moving_tool = eq_cwind.moving_tool) ) {
          eq_ckwind.panning = mer_ckwind.panning = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_MOVING_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_MOVING_TOOL );
        }
        else {
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_NO_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_NO_TOOL );
        }
        xge_Redraw ();
        break;
    case sw15aSCALE:
        if ( (mer_cwind.scaling_tool = eq_cwind.scaling_tool) ) {
          eq_ckwind.panning = mer_ckwind.panning = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_SCALING_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_SCALING_TOOL );
        }
        else {
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_NO_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_NO_TOOL );
        }
        xge_Redraw ();
        break;
    case sw15aROTATE:
        if ( (mer_cwind.rotating_tool = eq_cwind.rotating_tool) ) {
          eq_ckwind.panning = mer_ckwind.panning = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_ROTATING_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_ROTATING_TOOL );
        }
        else {
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_NO_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_NO_TOOL );
        }
        xge_Redraw ();
        break;
    case sw15aSHEAR:
        if ( (mer_cwind.shear_tool = eq_cwind.shear_tool) ) {
          eq_ckwind.panning = mer_ckwind.panning = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_SHEAR_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_SHEAR_TOOL );
        }
        else {
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_NO_TOOL );
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_NO_TOOL );
        }
        xge_Redraw ();
        break;
    case sw15aPAN_ZOOM:
        if ( (mer_ckwind.panning = mer_cwind.panning =
              eq_ckwind.panning = eq_cwind.panning) ) {
          eq_cwind.selecting_mode = eq_cwind.moving_tool =
          eq_cwind.scaling_tool = eq_cwind.rotating_tool =
          eq_cwind.shear_tool = false;
          xge_2DwindEnableGeomWidget ( &eq_cwind, xge_2DWIN_PANNING_TOOL );
          mer_cwind.selecting_mode = mer_cwind.moving_tool =
          mer_cwind.scaling_tool = mer_cwind.rotating_tool =
          mer_cwind.shear_tool = false;
          xge_2DwindEnableGeomWidget ( &mer_cwind, xge_2DWIN_PANNING_TOOL );
          xge_Redraw ();
        }
        break;
    case sw15aCOORDINATES:
        mer_ckwind.display_coord = eq_ckwind.display_coord =
        mer_cwind.display_coord = eq_cwind.display_coord;
        break;
    case sw15aMOVE_MANY_KNOTS:
        mer_ckwind.moving_many = eq_ckwind.moving_many;
        break;
    case sw15bCURVE:
    case sw15bCONTROL_POLYGON:
    case sw15bBEZIER_POLYGONS:
    case sw15bTICKS:
        xge_SetClipping ( eq_cwind.er );
        eq_cwind.er->redraw ( eq_cwind.er, true );
        break;
    case sw16STATUS:
        ResizeWinStatus ( win1 );
        xge_Redraw ();
        break;
    case sw112TRIHARMONIC_BLENDING:
        bind_spr = false;
        if ( sw_triharmonic_blending ) {
          if ( sw_blending_g1 ) {
            if ( degree_u == 2 && degree_v == 2 ) {
              blending_mat_valid = false;
              if ( ConstructG1BlendingSurface () ) {
                sw_bind_blending = true;
                swind_picture = false;
              }
              else
                sw_bind_blending = sw_triharmonic_blending = false;
            }
            else {
              sw_bind_blending = sw_triharmonic_blending = false;
              xge_DisplayErrorMessage ( ErrorMsgBiquadratic, -1 );
            }
          }
          else {
            if ( degree_u == 3 && degree_v == 3 ) {
              blending_mat_valid = false;
              if ( ConstructG2BlendingSurface () ) {
                sw_bind_blending = true;
                swind_picture = false;
              }
              else
                sw_bind_blending = sw_triharmonic_blending = false;
            }
            else {
              sw_bind_blending = sw_triharmonic_blending = false;
              xge_DisplayErrorMessage ( ErrorMsgBicubic, -1 );
            }
          }
        }
        else
          sw_bind_blending = false;
        xge_RedrawAll ();
        break;
    case sw112CLAMPED_BLENDING:
        if ( sw_blending_g1 ) {
          if ( sw_clamped_blending )
            G1FreeBoundaryToClamped ();
          else
            G1ClampedBoundaryToFree ();
        }
        else {
          if ( sw_clamped_blending )
            G2FreeBoundaryToClamped ();
          else
            G2ClampedBoundaryToFree ();
        }
        ResizeObject ();
        xge_RedrawAll ();
        break;
    case sw112CLOSED_BLENDING:
        if ( !sw_blending_g2 ) {
          if ( kwind.closed_u ) {
            kwind.closed_u = false;
            xge_DisplayErrorMessage ( ErrorMsgNotImplemented, -1 );
          }
        }
        else {
          UpdateBlendingRangeWidgets ();
          blending_mat_valid = false;
          if ( SetClosedUSurface () )
            RedrawGeomWindows ( true, true );
        }
        break;
    case sw112CONSTRAINTS:
        if ( !sw_blending_g2 ) {
          if ( sw_blending_constraints ) {
            sw_blending_constraints = false;
            xge_DisplayErrorMessage ( ErrorMsgNotImplemented, -1 );
          }
          break;
        }
        if ( sw_blending_constraints ) {
          xge_T2KnotWindSwitchAltKnots ( &kwind, true, false );
        }
        else {
          xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
        }
        sw_bind_blending = sw_triharmonic_blending = sw_nonlin_blending = false;
        xge_Redraw ();
        break;
    case sw112BLENDING_OPT_ENTIRE:
        if ( sw_blending_opt_entire ) {
          sw_triharmonic_blending = sw_nonlin_blending = false;
          if ( sw_blending_g1 ) {
            if ( degree_u == 2 && degree_v == 2 &&
                 lastknot_u >= 7 && lastknot_v >= 7 ) {
              blending_opt_part[0] = blending_opt_part[2] = 2;
              blending_opt_part[1] = blending_opt_range[0].maxvalue =
                blending_opt_range[1].maxvalue = lastknot_u-5;
              blending_opt_part[3] = blending_opt_range[2].maxvalue =
                blending_opt_range[3].maxvalue = lastknot_v-5;
            }
            else
              sw_blending_opt_entire = false;
          }
          else {
            if ( degree_u == 3 && degree_v == 3 &&
                 lastknot_u >= 10 && lastknot_v >= 10 ) {
              if ( kwind.closed_u ) {
                blending_opt_range[0].minvalue = blending_opt_range[1].minvalue = -1;
                blending_opt_range[0].maxvalue = blending_opt_range[1].maxvalue = lastknot_u-5;
                blending_opt_part[0] = 0;
              }
              else {
                blending_opt_range[0].minvalue = blending_opt_range[1].minvalue = 3;
                blending_opt_range[0].maxvalue = blending_opt_range[1].maxvalue = lastknot_u-7;
                blending_opt_part[0] = 3;
              }
              blending_opt_part[1] = lastknot_u-7;
              blending_opt_range[2].minvalue = blending_opt_range[3].minvalue = 3;
              blending_opt_range[2].maxvalue = blending_opt_range[3].maxvalue = lastknot_v-7;
              blending_opt_part[2] = 3;
              blending_opt_part[3] = lastknot_v-7;
            }
            else
              sw_blending_opt_entire = false;
          }
          blending_mat_valid = false;
          xge_RedrawAll ();
        }
        break;
    case swP12DUMPDATA:
        break;
    case sw112BLENDING_G1:
        sw_blending_g2 = !sw_blending_g1;
        goto switchg1g2;
    case sw112BLENDING_G2:
        sw_blending_g1 = !sw_blending_g2;
switchg1g2:
        if ( sw_blending_g1 )
          bl_trihsw->data0 = txtBiharmonic;
        else
          bl_trihsw->data0 = txtTriharmonic;
        sw_bind_blending = blending_mat_valid =
        sw_triharmonic_blending = sw_nonlin_blending =
        sw_blending_constraints = false;
        xge_SetClipping ( menu12 );
        menu12->redraw ( menu12, true );
        break;
    default:
        break;
      }
      break;

  case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
    case intw13DEGREE_U: /* Degree u elevation or reduction */
        if ( key > degree_u ) {
          if ( key > MAX_DEGREE )
            goto deg_too_high;
          if ( DegreeElevationU () ) {
            SetKWindNKnots ();
            goto redraw_it;
          }
        }
        else if ( key < degree_u ) {
          if ( key < 1 )
            goto deg_too_low;
          if ( DegreeReductionU () ) {
            goto redraw_it;
            SetKWindNKnots ();
          }
        }
        break;
    case intw13DEGREE_V: /* Degree v elevation or reduction */
        if ( key > degree_v ) {
          if ( key > MAX_DEGREE ) {
deg_too_high:
            xge_DisplayErrorMessage ( ErrorMsgCannotRaiseDeg, -1 );
            return 0;
          }
          else if ( DegreeElevationV () ) {
            SetKWindNKnots ();
            goto redraw_it;
          }
        }
        else if ( key < degree_v ) {
          if ( key < 1 ) {
deg_too_low:
            xge_DisplayErrorMessage ( ErrorMsgCannotReduceDeg, -1 );
            return 0;
          }
          else if ( DegreeReductionV () ) {
            SetKWindNKnots ();
redraw_it:
            swind_picture = false;
            if ( RenderingIsOn )
              BreakRendering ( false );
            xge_RedrawAll ();
          }
        }
        break;
    case intw15aDEGREE: /* equator od meridian degree elevation */
        if ( key > ckwind->degree ) {
          if ( key > MAX_DEGREE ) {
            xge_DisplayErrorMessage ( ErrorMsgCannotRaiseDeg, -1 );
            return 0;
          }
          else if ( EqMerDegreeElevation () )
            goto bind_spr_update2;
        }
        else if ( key < ckwind->degree ) {
          if ( key < 1 ) {
            xge_DisplayErrorMessage ( ErrorMsgCannotReduceDeg, -1 );
            return 0;
          }
          else if ( EqMerDegreeReduction () ) {
bind_spr_update2:
            if ( bind_spr ) {
              BindSphericalProduct ();
              SetKWindNKnots ();
              swind_picture = false;
              if ( RenderingIsOn )
                BreakRendering ( false );
              xge_RedrawAll ();
            }
            else
              xge_Redraw ();
          }
        }
        break;
    case intwP12BLENDING_LMT_NITER:
        blending_lmt_iter = key;
        break;
    case intwP12BLENDING_QUAD1:
        blending_quad1 = key;
        if ( blending_quad1 > blending_quad2 ) {
          blending_quad2 = blending_quad1;
          xge_SetClipping ( wdg_blending_quad2.er );
          wdg_blending_quad2.er->redraw ( wdg_blending_quad2.er, true );
        }
        break;
    case intwP12BLENDING_QUAD2:
        blending_quad2 = key;
        if ( blending_quad2 < blending_quad1 ) {
          blending_quad1 = blending_quad2;
          xge_SetClipping ( wdg_blending_quad1.er );
          wdg_blending_quad1.er->redraw ( wdg_blending_quad1.er, true );
        }
        break;
    case intw112BLENDING_UMIN:
        if ( kwind.closed_u ) {
          if ( key < 0 )
            blending_opt_part[0] = key + kwind.clcKu;
          else if ( key >= kwind.clcKu )
            blending_opt_part[0] = key - kwind.clcKu;
          else
            blending_opt_part[0] = key;
        }
        else {
          if ( key < blending_opt_part[0] )
            blending_opt_part[0] = key;
          else if ( key > blending_opt_part[0] && key < lastknot_u-2*degree_u ) {
            blending_opt_part[0] = key;
            blending_opt_part[1] = max ( blending_opt_part[1], key );
          }
        }
        goto redraw_opt_range;
    case intw112BLENDING_UMAX:
        if ( kwind.closed_u ) {
          if ( key < 0 )
            blending_opt_part[1] = key + kwind.clcKu;
          else if ( key >= kwind.clcKu )
            blending_opt_part[1] = key - kwind.clcKu;
          else
            blending_opt_part[1] = key;
        }
        else {
          if ( key > blending_opt_part[1] && key < lastknot_u-2*degree_u )
            blending_opt_part[1] = key;
          else if ( key < blending_opt_part[1] && key >= degree_u ) {
            blending_opt_part[1] = key;
            blending_opt_part[0] = min ( blending_opt_part[0], key );
          }
        }
        goto redraw_opt_range;
    case intw112BLENDING_VMIN:
        if ( key < blending_opt_part[2] )
          blending_opt_part[2] = key;
        else if ( key > blending_opt_part[2] && key < lastknot_v-2*degree_v ) {
          blending_opt_part[2] = key;
          blending_opt_part[3] = max ( blending_opt_part[3], key );
        }
        goto redraw_opt_range;
    case intw112BLENDING_VMAX:
        if ( key > blending_opt_part[3] && key < lastknot_v-2*degree_v )
          blending_opt_part[3] = key;
        else if ( key < blending_opt_part[3] && key >= degree_v ) {
          blending_opt_part[3] = key;
          blending_opt_part[2] = min ( blending_opt_part[2], key );
        }
redraw_opt_range:
        if ( kwind.closed_u ) {
          sw_blending_opt_entire =
              (blending_opt_part[2] == degree_v) &&
              (blending_opt_part[3] == lastknot_v-2*degree_v-1);
          if ( (blending_opt_part[0]-blending_opt_part[1] != 1) &&
               (blending_opt_part[0] != 0 || blending_opt_part[1] != kwind.clcKu-1) )
            sw_blending_opt_entire = false;
        }
        else {
          sw_blending_opt_entire =
              (blending_opt_part[0] == degree_u) &&
              (blending_opt_part[1] == lastknot_u-2*degree_u-1) &&
              (blending_opt_part[2] == degree_v) &&
              (blending_opt_part[3] == lastknot_v-2*degree_v-1);
        }
        sw_triharmonic_blending = sw_nonlin_blending = false;
        blending_mat_valid = false;
        xge_RedrawAll ();
        break;
    default:
        break;
      }
      break;

  case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
    case sl112NONLIN_BLENDING_C:
        NotifyDoubleNumber ( win1, xge_LogSlidebarValued (
                                     NLBLENDING_CMIN, NLBLENDING_CMAX,
                                     blending_factor ), true );
        break;
    default:
        break;
      }
      break;

  case xgemsg_DIAL_COMMAND:
      switch ( er->id ) {
    case dial15cARC_ANGLE:
        EqMerSetCircleArc ( arc_angle );
        goto bind_spr_update2;
    default:
        break;
      }
      break;

  case xgemsg_T2KNOTWIN_CHANGE_KNOT_U:
  case xgemsg_T2KNOTWIN_CHANGE_KNOT_V:
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( false );
      if ( bind_spr ) {
        memcpy ( equator_knots, knots_u, (lastknot_u+1)*sizeof(double) );
        memcpy ( meridian_knots, knots_v, (lastknot_v+1)*sizeof(double) );
      }
      RedrawGeomWindows ( true, false );
      break;

  case xgemsg_T2KNOTWIN_INSERT_KNOT_U:
      if ( !key ) {  /* incorrect place, to be handled by the application */
        xge_DisplayErrorMessage ( ErrorMsgCannotInsertKnot, -1 );
        return 0;
      }
      if ( InsertKnotU ( kwind.newknot ) ) {
        SetKWindNKnots ();
        UpdateBlendingRangeWidgets ();
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        RedrawGeomWindows ( true, false );
        return 1;
      }
      else
        return 0;

  case xgemsg_T2KNOTWIN_INSERT_KNOT_V:
      if ( !key ) {  /* incorrect place, to be handled by the application */
        xge_DisplayErrorMessage ( ErrorMsgCannotInsertKnot, -1 );
        return 0;
      }
      if ( InsertKnotV ( kwind.newknot ) ) {
        SetKWindNKnots ();
        UpdateBlendingRangeWidgets ();
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        RedrawGeomWindows ( true, false );
        return 1;
      }
      else
        return 0;

  case xgemsg_T2KNOTWIN_REMOVE_KNOT_U:
      if ( key ) {  /* key is nonzero if the conditions for knot removal */
                    /* checked by the knot window are satisfied */
        if ( RemoveKnotU ( kwind.current_item ) ) {
          SetKWindNKnots ();
          UpdateBlendingRangeWidgets ();
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          RedrawGeomWindows ( true, false );
          return 1;
        }
        else
          return 0;
      }
      else {
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveKnot, -1 );
        return 0;
      }

  case xgemsg_T2KNOTWIN_REMOVE_KNOT_V:
      if ( key ) {  /* key is nonzero if the conditions for knot removal */
                    /* checked by the knot window are satisfied */
        if ( RemoveKnotV ( kwind.current_item ) ) {
          SetKWindNKnots ();
          UpdateBlendingRangeWidgets ();
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          RedrawGeomWindows ( true, false );
          return 1;
        }
        else
          return 0;
      }
      else {
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveKnot, -1 );
        return 0;
      }

  case xgemsg_T2KNOTWIN_INSERT_ALTKNOT_U:
      if ( InsertG2BlendingConstrKnot ( kwind.newknot ) ) {
        RedrawGeomWindows ( true, true );
        return 1;
      }
      else {
        xge_DisplayErrorMessage ( ErrorMsgCannotAddConstraint, -1 );
        return 0;
      }
      break;

  case xgemsg_T2KNOTWIN_REMOVE_ALTKNOT_U:
      if ( RemoveG2BlendingConstrKnot ( kwind.current_item ) ) {
        RedrawGeomWindows ( true, true );
        return 1;
      }
      break;

  case xgemsg_T2KNOTWIN_CHANGE_ALTKNOT_U:
      ChangeG2BlendingConstrKnot ( key, kwind.current_item );
      RedrawGeomWindows ( true, true );
      return 1;

  case xgemsg_T2KNOTWIN_SELECT_POINTS:
      if ( display_domain_net ) {
        SelectDomPoints ( &kwind.selection_rect, true );
        if ( display_control_net )
          RedrawGeomWindows ( true, false );
        return 1;
      }
      else
        return 0;

  case xgemsg_T2KNOTWIN_UNSELECT_POINTS:
      if ( display_domain_net ) {
        SelectDomPoints ( &kwind.selection_rect, false );
        if ( display_control_net )
          RedrawGeomWindows ( true, false );
        return 1;
      }
      else
        return 0;

  case xgemsg_T2KNOTWIN_CHANGE_MAPPING:
      break;

  case xgemsg_T2KNOTWIN_ERROR:
      switch ( key ) {
    case 0:
    case 2:
        xge_DisplayErrorMessage ( ErrorMsgCannotInsertKnot, -1 );
        break;
    case 1:
    case 3:
        xge_DisplayErrorMessage ( ErrorMsgTooManyKnots, -1 );
        break;
/* not sent any more in this form */
/*
    case 4:
    case 5:
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveKnot, -1 );
        break;
*/
    default:
        break;
      }
      break;

  case xgemsg_2DWIN_RESIZE:
  case xgemsg_2DWIN_PROJCHANGE:
      ProjectEqMerCurve ();
      break;

  case xgemsg_2DWIN_PICK_POINT:
      return FindNearestEqMerCPoint ( x, y, xge_MINDIST );

  case xgemsg_2DWIN_MOVE_POINT:
      EqMerSetCPoint ( x, y );
      goto bind_spr_update2;

  case xgemsg_2DWIN_SELECT_POINTS:
      SelectEqMerPoints ( &cwind->selection_rect, true );
      break;

  case xgemsg_2DWIN_UNSELECT_POINTS:
      SelectEqMerPoints ( &cwind->selection_rect, false );
      break;

  case xgemsg_2DWIN_SAVE_POINTS:
      SaveEqMerControlPoints ();
      break;

  case xgemsg_2DWIN_TRANSFORM_POINTS:
      TransformEqMerMarkedControlPoints ( &cwind->gwtrans );
      goto bind_spr_update2;

  case xgemsg_2DWIN_FIND_REFBBOX:
      FindEqMerRefBox ( &cwind->RefBBox );
      break;

  case xgemsg_2DWIN_ERROR:
      break;

  case xgemsg_KNOTWIN_CHANGE_KNOT:
      goto bind_spr_update2;

  case xgemsg_KNOTWIN_INSERT_KNOT:
      InsertEqMerKnot ();
      goto bind_spr_update2;

  case xgemsg_KNOTWIN_REMOVE_KNOT:
      RemoveEqMerKnot ();
      goto bind_spr_update2;

  case xgemsg_KNOTWIN_CHANGE_MAPPING:
      break;

  case xgemsg_KNOTWIN_ERROR:
      switch ( key ) {
    case 0:
        xge_DisplayErrorMessage ( ErrorMsgCannotInsertKnot, -1 );
        break;
    case 1:
        xge_DisplayErrorMessage ( ErrorMsgTooManyKnots, -1 );
        break;
    case 2:
        xge_DisplayErrorMessage ( ErrorMsgCannotRemoveKnot, -1 );
        break;
    default:
        break;
      }
      break;

  case xgemsg_TEXT_EDIT_VERIFY:
      if ( er->id >= txtedP11A0 && er->id <= txtedP11A8 )
        return IsValidDoubleString ( er->data0 );
      else
        return true;

  case xgemsg_TEXT_EDIT_ENTER:
      if ( er->id >= txtedP11A0 && er->id <= txtedP11A8 ) {
        if ( !EnterBlendingOptTransCoefficient ( er->id ) )
          xge_DisplayErrorMessage ( ErrorMsgIncorrectNumber, -1 );
        else {
          xge_SetWindow ( win0 );
          xge_Redraw ();
          xge_SetWindow ( win1 );
        }
      }
      break;

  default:
      break;
    }
  }
  else
    return ProcessOtherMsg ( msg, key, x, y );
  return 1;
} /*Win1CallBack*/

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int win;

  win = xge_CurrentWindow ();
  if ( win == win0 )
    return Win0CallBack ( er, msg, key, x, y );
  else if ( win == win1 )
    return Win1CallBack ( er, msg, key, x, y );
  return 0;
} /*CallBack*/

