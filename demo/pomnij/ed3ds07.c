
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

/* ///////////////////////////////////////////////////////////////////////// */
int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int i;

  if ( er ) {
    switch ( msg ) {
  case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
    case btn00aRESET: /* Reset */
        swind_picture = false;
/*
        if ( xge_ChildIsActive () )
          xge_CallTheChild ( cmdTERMINATE, 0 );
*/
        if ( RenderingIsOn )
          BreakRendering ( false );
        BreakNLBlending ();
        ResetObject ();
        UpdateBlendingRangeWidgets ();
        sw_blending_constraints = false;
        xge_T2KnotWindSwitchAltKnots ( &kwind, false, false );
        SetKWindNKnots ();
        xge_3DwindResetGeomWidgetPos ( &swind );
        ClearPointMarking ( (lastknot_u-degree_u)*(lastknot_v-degree_v), mkpoints );
        xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        FindBoundingBox ( &swind.RefBBox );
        win1_contents = WIN1_GENERAL;
        domwind->next = menu3;
        menu3->prev = domwind;
        xge_SetWindow ( win1 );
        xge_SetWinEdRect ( menu7 );
        xge_RedrawAll ();
        break;

    case btn01ABOUT: /* About */
        xge_DisplayInfoMessage ( InfoMsg, -1 );
        break;

    case btn01FILEpopup: /* File */
        xge_AddPopup ( popup00 );
        xge_GrabFocus ( popup00, true );
        break;

    case btn02OPENpopup: /* Open - open box */
        xge_RemovePopup ( true );
        getcwd ( current_directory, MAX_PATH_LGT+1 );
        xge_SetupFileList ( &filelist, ".", file_filter );
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup01 );
        xge_GrabFocus ( popup01, true );
        break;

    case btn02SAVEpopup: /* Save */
        xge_RemovePopup ( true );
process_save_command:
        if ( FilenameCorrect ( filename ) ) {
          if ( !SaveBSPatch ( filename ) )
            xge_DisplayErrorMessage ( ErrorMsgCannotSave, -1 );
        }
        else
          goto open_popup02;
        break;

    case btn02SAVEASpopup: /* Save as - open box */
        xge_RemovePopup ( true );
open_popup02:
        getcwd ( current_directory, MAX_PATH_LGT+1 );
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup02 );
        xge_GrabFocus ( popup02, true );
        break;

    case btn02EXPORT: /* Export */
        xge_RemovePopup ( true );
/*
        CExport ();
*/
        if ( PovRayExport () ) {
          if ( degree_u > 3 || degree_v > 3 )
            xge_DisplayWarningMessage ( WarningMsgPovExport, -1 );
        }
        else
          xge_DisplayErrorMessage ( ErrorMsgFileWritingError, -1 );
        break;

    case btn02EXIT: /* Exit - open the Exit popup menu */
        xge_RemovePopup ( true );
        xge_AddPopup ( popup04 );
        xge_GrabFocus ( popup04, true );
/*        xge_done = 1; */
        break;

    case btn03OPEN:   /* Open */
        BreakNLBlending ();
        OpenFile ();
        break;

    case btn03CANCEL:  /* Cancel open */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &filelist );
        xge_ClearListBox ( &dirlist );
        break;

    case btn01EDITpopup:  /* Edit */
        xge_AddPopup ( popup03 );
        xge_GrabFocus ( popup03, true );
        break;

    case btn05SURFACE:
        xge_RemovePopup ( true );
        xge_SetMenuWidgets ( menu0, menu00list, true );
        break;

    case btn05LIGHT:
        xge_RemovePopup ( true );
        xge_SetMenuWidgets ( menu0, menu03list, true );
        break;

    case btn01VIEW:  /* View */
        xge_SetMenuWidgets ( menu0, menu01list, true );
        break;

    case btn01PICTURE:
        xge_SetMenuWidgets ( menu0, menu02list, true );
        break;

    case btn04SAVE:  /* Save */
            /* **************** */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &dirlist );
        if ( !SaveBSPatch ( filename ) )
          xge_DisplayErrorMessage ( ErrorMsgCannotSave, -1 );
        break;

    case btn04CANCEL:  /* Cancel save as */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &dirlist );
        break;

    case btn00cRENDER_STOP:
    case btn00dRENDER_STOP:
        if ( RenderingIsOn )
          BreakRendering ( true );
        else {
          if ( swind.fww.zoomwin != -1 && swind.fww.zoomwin != 3 ) {
            swind.fww.zoomwin = 3;
            xge_CompSizeFourWW ( swind.fww.er, 1 );
            xge_3DwindUpdatePerspProj ( &swind );
          }
          renderbtn0->data0 = renderbtn1->data0 = txtInterrupt;
          if ( menu0->data1 == menu02list ) {
            xge_SetClipping ( renderbtn0 );
            renderbtn0->redraw ( renderbtn0, true );
          }
          else if ( menu0->data1 == menu03list ) {
            xge_SetClipping ( renderbtn1 );
            renderbtn1->redraw ( renderbtn1, true );
          }
          xge_PostIdleCommand ( IDLE_COMMAND_START_RENDERING, 0, 0 );
        }
        break;


    case btn06EXIT:
        xge_done = 1;
        break;

    case btn06SAVE:
        xge_RemovePopup ( true );
        goto process_save_command;

    case btn06CANCEL:
        xge_RemovePopup ( true );
        break;

    default:
        break;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
    case sw00bCONTROLNET:
    case sw00bSURFACE:
    case sw00bBEZIER_NETS:
        xge_SetWindow ( win0 );
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        break;
    case sw00bCONSTRPOLY:
        if ( win1_contents == WIN1_BLENDING ) {
          if ( display_constr_poly ) {
            for ( i = 0; i < 4; i++ )
              ProjectG2BlendingConstrCPoly ( i );;
          }
          xge_SetWindow ( win0 );
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
        }
        else if ( display_constr_poly ) {
          display_constr_poly = false;
          xge_SetWindow ( win0 );
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        break;
    case sw00bTRANSFNET:
        if ( win1_contents == WIN1_BLENDING ) {
          xge_SetWindow ( win0 );
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
        }
        else if (display_trans_net ) {
          display_trans_net = false;
          xge_SetWindow ( win0 );
          xge_SetClipping ( er );
          er->redraw ( er, true );
        }
        break;
    case sw00aMOVE:  /* move */
        if ( swind.moving_tool )
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_MOVING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        goto redraw_it;
    case sw00aSCALE:  /* scale */
        if ( swind.scaling_tool )
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SCALING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        goto redraw_it;
    case sw00aROTATE:  /* rotate */
        if ( swind.rotating_tool )
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_ROTATING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        goto redraw_it;
    case sw00aSHEAR:  /* shear */
        if ( swind.shear_tool )
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SHEAR_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        goto redraw_it;
    case sw00aMARK: /* mark/unmark */
        if ( swind.selecting_mode ) {
          if ( !display_constr_poly )
            display_control_net = true;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SELECTING_TOOL );
        }
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
redraw_it:
        xge_3DwindResetGeomWidgets ( &swind );
        xge_SetWindow ( win0 );
        xge_Redraw ();
        break;
    case sw00aPAN_ZOOM: /* panning */
    case sw00bPAN_ZOOM:
        if ( swind.panning )
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_PANNING_TOOL );
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_SetWindow ( win0 );
        xge_Redraw ();
        break;
    case sw02STATUS:
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        ResizeWinStatus ( win0 );
        xge_Redraw ();
        break;

    case sw00cGAUSSIAN_C:
        if ( swGaussian_c )
          swMean_c = swLambert_c = swReflection_c = swHighlight_c =
          swSections_c = swParam_c = false;
        goto redraw_menu;
    case sw00cMEAN_C:
        if ( swMean_c )
          swGaussian_c = swLambert_c = swReflection_c = swHighlight_c =
          swSections_c = swParam_c = false;
        goto redraw_menu;
    case sw00cLAMBERTISO_C:
        if ( swLambert_c )
          swGaussian_c = swMean_c = swReflection_c = swHighlight_c =
          swSections_c = swParam_c = false;
        goto redraw_menu;
    case sw00cREFLECTION_C:
        if ( swReflection_c )
          swGaussian_c = swMean_c = swLambert_c = swHighlight_c =
          swSections_c = swParam_c = false;
        goto redraw_menu;
    case sw00cHIGHLIGHT_C:
        if ( swHighlight_c )
          swGaussian_c = swMean_c = swLambert_c = swReflection_c =
          swSections_c = swParam_c = false;
        goto redraw_menu;
    case sw00cSECTIONS_C:
        if ( swSections_c )
          swGaussian_c = swMean_c = swLambert_c = swReflection_c =
          swHighlight_c = swParam_c = false;
        goto redraw_menu;
    case sw00cPARAM_C:
        if ( swParam_c )
          swGaussian_c = swMean_c = swLambert_c = swReflection_c =
          swHighlight_c = swSections_c = false;
        goto redraw_menu;
    case sw00cGAUSSIAN_D:
        if ( swGaussian_d )
          swMean_d = swLambert_d = swReflection_d = swHighlight_d =
          swSections_d = swParam_d = false;
        goto redraw_menu;
    case sw00cMEAN_D:
        if ( swMean_d )
          swGaussian_d = swLambert_d = swReflection_d = swHighlight_d =
          swSections_d = swParam_d = false;
        goto redraw_menu;
    case sw00cLAMBERTISO_D:
        if ( swLambert_d )
          swGaussian_d = swMean_d = swReflection_d = swHighlight_d =
          swSections_d = swParam_d = false;
        goto redraw_menu;
    case sw00cREFLECTION_D:
        if ( swReflection_d )
          swGaussian_d = swMean_d = swLambert_d = swHighlight_d =
          swSections_d = swParam_d = false;
        goto redraw_menu;
    case sw00cHIGHLIGHT_D:
        if ( swHighlight_d )
          swGaussian_d = swMean_d = swLambert_d = swReflection_d =
          swSections_d = swParam_d = false;
        goto redraw_menu;
    case sw00cSECTIONS_D:
        if ( swSections_d )
          swGaussian_d = swMean_d = swLambert_d = swReflection_d =
          swHighlight_d = swParam_d = false;
        goto redraw_menu;
    case sw00cPARAM_D:
        if ( swParam_d )
          swGaussian_d = swMean_d = swLambert_d = swReflection_d =
          swHighlight_d = swSections_d = false;
redraw_menu:
        xge_SetClipping ( menu1 );
        menu0->redraw ( menu0, true );
        return 1;
    case sw00cSHADOWS:
    case sw00cANTIALIAS:
        return 1;

    case sw00dLIGHT0DIR:
    case sw00dLIGHT1DIR:
    case sw00dLIGHT2DIR:
    case sw00dLIGHT3DIR:
    case sw00dREFLECTIONFRAME:
    case sw00dHIGHLIGHTFRAME:
    case sw00dSECTIONSFRAME:
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        return 1;

    default:
        break;
      }
      break;

  case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
    case sl00cREND_DFSF:
        break;
    default:
        break;
      }
      break;

  case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
    case intw00bU_DENSITY:
        display_bez_dens_u = key;
        RedrawGeomWindows ( true, false );
        break;
    case intw00bV_DENSITY:
        display_bez_dens_v = key;
        RedrawGeomWindows ( true, false );
        break;
    default:
        break;
      }
      break;

  case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
    case lb03DIRLIST:  /* change directory - popup 0 3 */
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupFileList ( &filelist, ".", file_filter );
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup01 );
          popup01->redraw ( popup01, true );
        }
        break;

    case lb03FILELIST:  /* open file */
        BreakNLBlending ();
        OpenFile ();
        break;

    case lb04DIRLIST:  /* change directory - popup 0 4 */
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup01 );
          popup02->redraw ( popup02, true );
        }
        break;

    default:
        break;
      }
      break;

  case xgemsg_3DWIN_RESIZE:
  case xgemsg_3DWIN_PROJCHANGE:
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( false );
      for ( i = 0; i < 4; i++ )
        ProjectSurface ( i );
      if ( display_constr_poly ) {
        for ( i = 0; i < 4; i++ )
          ProjectG2BlendingConstrCPoly ( i );;
      }
      break;

  case xgemsg_3DWIN_PICK_POINT:
      return FindNearestCPoint ( er->id, x, y );

  case xgemsg_3DWIN_MOVE_POINT:
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( false );
      SetCPoint ( er->id, x, y );
      if ( sw_nonlin_blending ) {
        sw_bind_blending = sw_nonlin_blending = sw_triharmonic_blending = false;
        BreakNLBlending ();
        xge_RedrawAll ();
        xge_SetWindow ( win0 );
      }
      else
        BreakSprBind ();
      break;

  case xgemsg_3DWIN_SELECT_POINTS:
      SelectPoints ( er->id, &swind.selection_rect, true );
      RedrawGeomWindows ( false, true );
      break;

  case xgemsg_3DWIN_UNSELECT_POINTS:
      SelectPoints ( er->id, &swind.selection_rect, false );
      RedrawGeomWindows ( false, true );
      break;

  case xgemsg_3DWIN_CHANGE_TRANS:
      Notify3DTransChange ( win0, er->data1 );
      break;

  case xgemsg_3DWIN_SAVE_POINTS:
      SaveControlPoints ();
      break;

  case xgemsg_3DWIN_TRANSFORM_POINTS:
      if ( sw_nonlin_blending ) {
        sw_bind_blending = sw_nonlin_blending = sw_triharmonic_blending = false;
        BreakNLBlending ();
        xge_RedrawAll ();
        xge_SetWindow ( win0 );
      }
      Notify3DTrans ( win0, er->data1 );
      TransformMarkedControlPoints ( &swind.gwtrans );
      BreakSprBind ();
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( false );
      break;

  case xgemsg_3DWIN_FIND_REFBBOX:
      FindBoundingBox ( &swind.RefBBox );
      break;

  case xgemsg_3DWIN_ERROR:
      return 0;

  default:
      break;
    }
  }
  else
    return ProcessOtherMsg ( msg, key, x, y );
  return 1;
} /*Win0CallBack*/

