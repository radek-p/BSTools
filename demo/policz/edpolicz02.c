
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
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
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"
#include "datagend.h"


int ProcessKey ( int key )
{
  switch ( key ) {
case 'm':
    xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
    return 1;
case 'M':
    xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
    return 1;
/*
case 'q':  case 'Q':
    xge_done = 1;
    return 1;
*/
default:
    return 0;
  }
} /*ProcessKey*/

void ProcessIdleCommand ( int key, short x, short y )
{
  switch ( key ) {
case IDLE_COMMAND_FIND_SURFACE1:
    if ( UpdateFinalSurface1 () )
      goto redraw;
    break;
case IDLE_COMMAND_FIND_SURFACE2:
    if ( !UpdateFinalSurface2 () )
      break;
redraw:
    xge_SetWindow ( win1 );
    xge_SetCurrentWindowCursor ( None );
    xge_SetWindow ( win0 );
    xge_SetCurrentWindowCursor ( None );
    xge_SetClipping ( swind.fww.er );
    swind.fww.er->redraw ( swind.fww.er, true );
    break;
case IDLE_COMMAND_START_RENDERING:
    StartRendering ();
    break;
case IDLE_COMMAND_RENDER:
    ContRendering ();
    break;
default:
    break;
  }
} /*ProcessIdleCommand*/

boolean SetupHoleSides ( int key )
{
  if ( InitGHObject ( key ) ) {
    KnotWindSetHoleK ( &knwind, hole_k, knots );
    ConfigureConstraintWidgets ( true );
    ResizeWinStatus ( win1 );
    view_surf_1 = view_surf_2 = swind_picture = false;
    options1.constr_matrix_valid = options2.constr_matrix_valid = false;
    if ( RenderingIsOn )
      BreakRendering ( false );
    xge_RedrawAll ();
    return 1;
  }
  else
    return 0;
} /*SetupHoleSides*/

void UpdateSurfPicture ( void )
{
  boolean ch;

  InvalFinalSurfaces ();
  CountTheConstraints ();
  ch = false;
  if ( view_surf_1 ) {
    if ( SurfFast ( &options1 ) ) {
      UpdateFinalSurface1 ();
    }
    else {
      view_surf_1 = false;
      ch = true;
    }
  }
  if ( view_surf_2 ) {
    if ( SurfFast ( &options2 ) ) {
      UpdateFinalSurface2 ();
    }
    else {
      view_surf_2 = false;
      ch = true;
    }
  }
  if ( ch && menu1->data1 == menu01dlist ) {
    xge_SetClipping ( menu1 );
    menu1->redraw ( menu1, true );
  }
} /*UpdateSurfPicture*/

static void UpdateConstrOpt ( void )
{
  GHoptions *opt;

  if ( constraints1 )
    opt = &options1;
  else
    opt = &options2;
  if ( constraints1st )
    opt->constr_type = 1;
  else if ( constraints2nd )
    opt->constr_type = 2;
  else
    opt->constr_type = 0;
  opt->constr_matrix_valid = false;
} /*UpdateConstrOpt*/

static void OpenTheFile ( void )
{
  xge_RemovePopup ( true );
  if ( xge_GetCurrentListBoxString ( &filelist, filename ) ) {
    xge_ClearListBox ( &filelist );
    xge_ClearListBox ( &dirlist );
    OpenFile ( filename );
  }
} /*OpenTheFile*/

int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int surfno;

  if ( er ) {
    switch ( msg ) {
  case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
    case btn00FILE:
        xge_AddPopup ( popup00 );
        xge_GrabFocus ( popup00, true );
        return 1;
    case btn00DATA:
        EditSetSurface ();
        xge_SetMenuWidgets ( menu1, menu01alist, false );
        xge_Redraw ();
        return 1;
    case btn00EDIT:
        xge_AddPopup ( popup03 );
        xge_GrabFocus ( popup03, true );
        return 1;
    case btn00VIEW:
        EditSetSelectView ();
        xge_SetMenuWidgets ( menu1, menu01dlist, false );
        xge_Redraw ();
        return 1;
    case btn00PICTURE:
        xge_SetMenuWidgets ( menu1, menu01elist, true );
        return 1;
    case btn00ABOUT:
        xge_DisplayInfoMessage ( InfoMsg, -1 );
        return 1;
    case btnP00NEW:
        xge_RemovePopup ( true );
        xge_DisplayErrorMessage ( ErrMsgNoAction, -1 );
        return 1;
    case btnP00OPEN:
        xge_RemovePopup ( true );
        getcwd ( current_directory, MAX_PATH_LGT+1 );
        xge_SetupFileList ( &filelist, ".", file_filter, false );
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
        xge_AddPopup ( popup01 );
        xge_GrabFocus ( popup01, true );
        return 1;
    case btnP00SAVE:
        xge_RemovePopup ( true );
process_save_command:
        if ( FilenameCorrect ( filename ) ) {
          SaveFile ( filename );
          return 1;
        }
        else goto open_save_as;
    case btnP00SAVEAS:
        xge_RemovePopup ( true );
open_save_as:
        getcwd ( current_directory, MAX_PATH_LGT+1 );
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
        xge_AddPopup ( popup02 );
        xge_GrabFocus ( popup02, true );
        return 1;
    case btnP00EXIT:
        xge_RemovePopup ( true );
        xge_AddPopup ( popup04 );
        xge_GrabFocus ( popup04, true );
        return 1;
    case btnP01OPEN:
        OpenTheFile ();
        return 1;
    case btnP01CANCEL:
        xge_RemovePopup ( true );
        return 1;
    case btnP02SAVE:
        xge_RemovePopup ( true );
        SaveFile ( filename );
        return 1;
    case btnP02CANCEL:
        xge_RemovePopup ( true );
        return 1;
    case btn01eRENDER_STOP:
    case btn01fRENDER_STOP:
        if ( RenderingIsOn )
          BreakRendering ( true );
        else {
          renderbtn0->data0 = renderbtn1->data0 = txtStop;
          if ( menu1->data1 == menu01elist ) {
            xge_SetClipping ( renderbtn0 );
            renderbtn0->redraw ( renderbtn0, true );
          }
          else if ( menu1->data1 == menu01flist ) {
            xge_SetClipping ( renderbtn1 );
            renderbtn1->redraw ( renderbtn1, true );
          }
          xge_PostIdleCommand ( IDLE_COMMAND_START_RENDERING, 0, 0 );
        }
        return 1;
    case btnP03SURFACE:
        EditSetSurface ();
        xge_RemovePopup ( true );
        xge_SetMenuWidgets ( menu1, menu01blist, false );
        xge_Redraw ();
        return 1;
    case btnP03CONSTRAINTS:
        xge_RemovePopup ( true );
        EditSetConstraints ();
        xge_SetMenuWidgets ( menu1, menu01clist, false );
        xge_Redraw ();
        return 1;
    case btnP03LIGHT:
        xge_RemovePopup ( true );
        EditSetLight ();
        xge_SetMenuWidgets ( menu1, menu01flist, false );
        xge_Redraw ();
        return 1;
    case btn01cGET_CURRENT:
        if ( constraints1 )
          surfno = 1;
        else
          surfno = 2;
        if ( GetCurrentConstraintFrame ( surfno ) )
          xge_Redraw ();
        else
          xge_DisplayErrorMessage ( ErrMagBadConstrFrame, -1 );
        return 1;
  case btnP04EXIT:
        xge_done = 1;
        break;
  case btnP04SAVE:
        xge_RemovePopup ( true );
        goto process_save_command;
  case btnP04CANCEL:
        xge_RemovePopup ( true );
        break;
    default:
        return 0;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      if ( er->id >= sw01cCONSTR_SWITCH &&
           er->id < sw01cCONSTR_SWITCH+NUM_CONSTR_SWITCHES ) {
        SwitchTheConstraint ( er->id - sw01cCONSTR_SWITCH );
        UpdateSurfPicture ();
        xge_Redraw ();
        return 1;
      }
      else {
        switch ( er->id ) {
      case sw01bMARK_UNMARK:
      case sw01cMARK_UNMARK:
          if ( swind.selecting_mode )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SELECTING_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01bMOVE:
      case sw01cMOVE:
          if ( swind.moving_tool )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_MOVING_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01bSCALE:
      case sw01cSCALE:
          if ( swind.scaling_tool )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SCALING_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01bROTATE:
      case sw01cROTATE:
          if ( swind.rotating_tool )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_ROTATING_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01bSHEAR:
      case sw01cSHEAR:
          if ( swind.shear_tool )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SHEAR_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01dFIRST:
          if ( view_surf_1 && !final_cp1 ) {
            if ( SurfFast ( &options1 ) ) {
              UpdateFinalSurface1 ();
            }
            else {
              xge_SetWindow ( win1 );
              xge_SetCurrentWindowCursor ( xgecursor[5] );
              xge_SetWindow ( win0 );
              xge_SetCurrentWindowCursor ( xgecursor[5] );
              xge_PostIdleCommand ( IDLE_COMMAND_FIND_SURFACE1, 0, 0 );
              return 1;
            }
          }
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
          return 1;
      case sw01dSECOND:
          if ( view_surf_2 && !final_cp2 ) {
            if ( SurfFast ( &options2 ) ) {
              UpdateFinalSurface2 ();
            }
            else {
              xge_SetWindow ( win1 );
              xge_SetCurrentWindowCursor ( xgecursor[5] );
              xge_SetWindow ( win0 );
              xge_SetCurrentWindowCursor ( xgecursor[5] );
              xge_PostIdleCommand ( IDLE_COMMAND_FIND_SURFACE2, 0, 0 );
              return 1;
            }
          }
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
          return 1;
      case sw01dCPOINTS:
      case sw01dSURFACE:
      case sw01dNUMBERS:
      case sw01dCONSTRAINT_FRAME:
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
          return 1;
      case sw01PAN_ZOOM:
          if ( swind.panning )
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_PANNING_TOOL );
          else
            xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
          xge_3DwindResetGeomWidgets ( &swind );
          xge_Redraw ();
          return 1;
      case sw01COORDINATES:
          return 1;
      case sw02STATUS:
          swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( true );
          StatusLineOnOff ( win0 );
          return 1;

      case sw01eGAUSSIAN_C:
          if ( swGaussian_c )
            swMean_c = swLambert_c = swReflection_c = swHighlight_c =
            swSections_c = false;
          goto redraw_menu;
      case sw01eMEAN_C:
          if ( swMean_c )
            swGaussian_c = swLambert_c = swReflection_c = swHighlight_c =
            swSections_c = false;
          goto redraw_menu;
      case sw01eLAMBERTISO_C:
          if ( swLambert_c )
            swGaussian_c = swMean_c = swReflection_c = swHighlight_c =
            swSections_c = false;
          goto redraw_menu;
      case sw01eREFLECTION_C:
          if ( swReflection_c )
            swGaussian_c = swMean_c = swLambert_c = swHighlight_c =
            swSections_c = false;
          goto redraw_menu;
      case sw01eHIGHLIGHT_C:
          if ( swHighlight_c )
            swGaussian_c = swMean_c = swLambert_c = swReflection_c =
            swSections_c = false;
          goto redraw_menu;
      case sw01eSECTIONS_C:
          if ( swSections_c )
            swGaussian_c = swMean_c = swLambert_c = swReflection_c =
            swHighlight_c = false;
          goto redraw_menu;
      case sw01eGAUSSIAN_D:
          if ( swGaussian_d )
            swMean_d = swLambert_d = swReflection_d = swHighlight_d =
            swSections_d = false;
          goto redraw_menu;
      case sw01eMEAN_D:
          if ( swMean_d )
            swGaussian_d = swLambert_d = swReflection_d = swHighlight_d =
            swSections_d = false;
          goto redraw_menu;
      case sw01eLAMBERTISO_D:
          if ( swLambert_d )
            swGaussian_d = swMean_d = swReflection_d = swHighlight_d =
            swSections_d = false;
          goto redraw_menu;
      case sw01eREFLECTION_D:
          if ( swReflection_d )
            swGaussian_d = swMean_d = swLambert_d = swHighlight_d =
            swSections_d = false;
          goto redraw_menu;
      case sw01eHIGHLIGHT_D:
          if ( swHighlight_d )
            swGaussian_d = swMean_d = swLambert_d = swReflection_d =
            swSections_d = false;
          goto redraw_menu;
      case sw01eSECTIONS_D:
          if ( swSections_d )
            swGaussian_d = swMean_d = swLambert_d = swReflection_d =
            swHighlight_d = false;
  redraw_menu:
          xge_SetClipping ( menu1 );
          menu1->redraw ( menu1, true );
          return 1;
      case sw01eSHADOWS:
      case sw01eANTIALIAS:
          return 1;

      case sw01fLIGHT0DIR:
      case sw01fLIGHT1DIR:
      case sw01fLIGHT2DIR:
      case sw01fLIGHT3DIR:
      case sw01fREFLECTIONFRAME:
      case sw01fHIGHLIGHTFRAME:
      case sw01fSECTIONSFRAME:
          xge_SetClipping ( swind.fww.er );
          swind.fww.er->redraw ( swind.fww.er, true );
          return 1;

      case sw01c1ST_SURF_CONSTR:
          constraints2 = !constraints1;
          goto set_constr_type_sw;
      case sw01c2ND_SURF_CONSTR:
          constraints1 = !constraints2;
set_constr_type_sw:
          if ( constraints1 ) {
            switch ( options1.constr_type ) {
          case 0: constraints1st = constraints2nd = false;         break;
          case 1: constraints1st = true;  constraints2nd = false;  break;
          case 2: constraints1st = false;  constraints2nd = true;  break;
            }
          }
          else {
            switch ( options2.constr_type ) {
          case 0: constraints1st = constraints2nd = false;         break;
          case 1: constraints1st = true;  constraints2nd = false;  break;
          case 2: constraints1st = false;  constraints2nd = true;  break;
            }
          }
          ConfigureConstraintWidgets ( false );
          xge_Redraw ();
          return 1;
      case sw01c1ST_CONSTR:
          if ( constraints1st )
            constraints2nd = false;
          UpdateConstrOpt ();
          goto update_constr_menu;
      case sw01c2ND_CONSTR:
          if ( constraints2nd )
            constraints1st = false;
          UpdateConstrOpt ();
          goto update_constr_menu;
      case sw01cZERO_CONSTR:
update_constr_menu:
          if ( constraints1 )
            options1.constr_matrix_valid = false;
          else
            options2.constr_matrix_valid = false;
          UpdateSurfPicture ();
          xge_Redraw ();
          return 1;

      default:
          return 0;
        }
      }

  case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
    case sl01aPARAM0:
    case sl01aPARAM1:
    case sl01aPARAM2:
    case sl01aPARAM3:
    case sl01aPARAM4:
        InitGHSurfNetd ( hole_k, surfcvect, 5, surfcparam, hole_cp );
        InitConstraintFrame ( 1 );
        InitConstraintFrame ( 2 );
        ProjectSurfaceNet ();
        UpdateSurfPicture ();
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        return 1;
    case sl01fLIGHT0INT:
    case sl01fLIGHT1INT:
    case sl01fLIGHT2INT:
    case sl01fLIGHT3INT:
    case sl01fLIGHTAMB:
        return 1;
    default:
        break;
      }
      return 0;

  case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
    case intw01aHOLE_SIDES:
        if ( SetupHoleSides ( key ) ) {
          return 1;
        }
        else
          return 0;
    default:
        return 0;
      }

  case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
    case lbP01DIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupFileList ( &filelist, ".", file_filter, false );
          xge_SetupDirList ( &dirlist, ".", NULL, false, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup01 );
          popup01->redraw ( popup01, true );
        }
        break;
    case lbP01FILELIST:
        OpenTheFile ();
        break;
    case lbP02DIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupDirList ( &dirlist, ".", NULL, false, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup02 );
          popup02->redraw ( popup02, true );
        }
        break;
    default:
        break;
      }
      return 1;

  case xgemsg_TEXT_EDIT_VERIFY:
      return 1;

  case xgemsg_TEXT_EDIT_ENTER:
      return 1;

  case xgemsg_TEXT_EDIT_ESCAPE:
      return 1;

  case xgemsg_3DWIN_RESIZE:
  case xgemsg_3DWIN_PROJCHANGE:
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
      ProjectSurfaceNet ();
printf ( "w = %d, h = %d\n", swind.CPos[3].width, swind.CPos[3].height );
      return 1;

  case xgemsg_3DWIN_PICK_POINT:
      return FindNearestSWinPoint ( er->id, x, y, xge_MINDIST );

  case xgemsg_3DWIN_MOVE_POINT:
      SetSWinPoint ( er->id, x, y );
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
      UpdateSurfPicture ();
      return 1;

  case xgemsg_3DWIN_SELECT_POINTS:
      SelectSWinPoints ( er->id );
      return 1;

  case xgemsg_3DWIN_UNSELECT_POINTS:
      UnselectSWinPoints ( er->id );
      return 1;

  case xgemsg_3DWIN_CHANGE_TRANS:
      Notify3DTransChange ( win0, er->data1 );
      return 1;

  case xgemsg_3DWIN_SAVE_POINTS:
      SaveSWinPoints ();
      return 1;

  case xgemsg_3DWIN_TRANSFORM_POINTS:
      Notify3DTrans ( win0, er->data1 );
      TransformSWinPoints ();
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
      UpdateSurfPicture ();
      return 1;

  case xgemsg_3DWIN_FIND_REFBBOX:
      FindBoundingBox ( &swind.RefBBox );
      return 1;

  case xgemsg_3DWIN_ERROR:
      return 0;

  default:
      return 0;
    }
  }
  else {
    switch ( msg ) {
  case xgemsg_RESIZE:
      ResizeWinStatus ( win0 );
      swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
      if ( key )
        xge_Redraw ();
      return 1;

  case xgemsg_KEY:
      return ProcessKey ( key );

  case xgemsg_IDLE_COMMAND:
      ProcessIdleCommand ( key, x, y );
      return 1;

  default:
      return 0;
    }
  }
  return 0;
} /*Win0CallBack*/

int Win1CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
    switch ( msg ) {
  case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
    case btn10OPTIONS:
        xge_SetMenuWidgets ( menu4, menu14alist, true );
        return 1;
    case btn10DATA:
        xge_SetMenuWidgets ( menu4, menu14blist, true );
        return 1;
    case btn10EDIT:
        xge_SetMenuWidgets ( menu4, menu14clist, true );
        return 1;
    case btn10VIEW:
        xge_SetMenuWidgets ( menu4, menu14dlist, true );
        return 1;
    case btn10INFO:
        DisplayBasisInfo ();
        return 1;
    default:
        return 0;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
    case sw11cMARK_UNMARK:
        if ( domwind.selecting_mode ) {
          if ( !view_dom_cp )
            view_dom_cp = true;
          knwind.panning = domwind.panning = domwind.moving_tool =
          domwind.scaling_tool = domwind.rotating_tool =
          domwind.shear_tool = false;
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
          xge_Redraw ();
        }
        return 1;
    case sw11cMOVE:
        if ( domwind.moving_tool ) {
          knwind.panning = domwind.panning = false;
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_MOVING_TOOL );
        }
        else
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return 1;
    case sw11cSCALE:
        if ( domwind.scaling_tool ) {
          knwind.panning = domwind.panning = false;
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_SCALING_TOOL );
        }
        else
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return 1;
    case sw11cROTATE:
        if ( domwind.rotating_tool ) {
          knwind.panning = domwind.panning = false;
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_ROTATING_TOOL );
        }
        else
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return 1;
    case sw11cSHEAR:
        if ( domwind.shear_tool ) {
          knwind.panning = domwind.panning = false;
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_SHEAR_TOOL );
        }
        else
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        return 1;
    case sw11PAN_ZOOM:
        knwind.panning = domwind.panning;
        if ( domwind.panning ) {
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_PANNING_TOOL );
          xge_Redraw ();
        }
        else
          xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        return 1;
    case sw11COORDINATES:
        knwind.display_coord = domwind.display_coord;
        return 1;
    case sw15STATUS:
        StatusLineOnOff ( win1 );
        return 1;
    case sw11dDOMCPOINTS:
    case sw11dDOMSURRPATCHES:
    case sw11dDOMPATCHES1:
    case sw11dDOMPATCHES2:
    case sw11dDOMNUMBERS:
        xge_SetClipping ( domwind.er );
        domwind.er->redraw ( domwind.er, true );
        return 1;
    case sw11a1FIRST:
    case sw11a2FIRST:
        sw_opt_2 = !sw_opt_1;  /* do not break */
    case sw11a1SECOND:
    case sw11a2SECOND:
        sw_opt_1 = !sw_opt_2;
        if ( sw_opt_1 )
          menu14alist = menu14a1list;
        else
          menu14alist = menu14a2list;
        xge_SetMenuWidgets ( menu4, menu14alist, false );
        SetStatusLineText ( win1, "", false );
        xge_Redraw ();
        return 1;
    case sw11a1COONS:
        if ( options1.coons )
          options1.bezier = options1.spline = false;
        else
          options1.bezier = true;
        InitConstraintFrame ( 1 );
        goto inval_surf1;
    case sw11a1BEZIER:
        if ( options1.bezier )
          options1.coons = options1.spline = false;
        else
          options1.spline = true;
        InitConstraintFrame ( 1 );
        goto inval_surf1;
    case sw11a1SPLINE:
        if ( options1.spline )
          options1.coons = options1.bezier = false;
        else
          options1.coons = true;
        InitConstraintFrame ( 1 );
        goto inval_surf1;
    case sw11a1RESTRICTED:
    case sw11a1ALTCENTRE:
    case sw11a1GAUSSLEGENDRE:
        goto inval_surf1;
    case sw11a2COONS:
        if ( options2.coons )
          options2.bezier = options2.spline = false;
        else
          options2.bezier = true;
        InitConstraintFrame ( 2 );
        goto inval_surf2;
    case sw11a2BEZIER:
        if ( options2.bezier )
          options2.coons = options2.spline = false;
        else
          options2.spline = true;
        InitConstraintFrame ( 2 );
        goto inval_surf2;
    case sw11a2SPLINE:
        if ( options2.spline )
          options2.coons = options2.bezier = false;
        else
          options2.coons = true;
        InitConstraintFrame ( 2 );
        goto inval_surf2;
    case sw11a2RESTRICTED:
    case sw11a2ALTCENTRE:
    case sw11a2GAUSSLEGENDRE:
        goto inval_surf2;
    case sw11a1LINEAR:
        goto inval_surf1;
    case sw11a1QUASIG2:
        SetStatusLineText ( win1, "", true );
        if ( options1.order != 1 )
          Option1SetOrder ( 1 );
inval_surf1:
        UpdateDomain1 ();
        ConfigureConstraintWidgets ( false );
        if ( view_surf_1 ) {
          view_surf_1 = swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
        }
        xge_RedrawAll ();
        return 1;
    case sw11a2LINEAR:
        goto inval_surf2;
    case sw11a2QUASIG2:
        SetStatusLineText ( win1, "", true );
        if ( options2.order != 1 )
          Option2SetOrder ( 1 );
inval_surf2:
        UpdateDomain2 ();
        ConfigureConstraintWidgets ( false );
        if ( view_surf_2 ) {
          view_surf_2 = swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( false );
        }
        xge_RedrawAll ();
        return 1;
    default:
        return 0;
      }

  case xgemsg_SLIDEBAR_COMMAND:
      switch ( er->id ) {
    case sl11bPARAM0:
    case sl11bPARAM1:
    case sl11bPARAM2:
    case sl11bPARAM3:
        InvalFinalSurfaces ();
        InitGHDomainNetd ( hole_k, domcvect, 3, domcparam, domain_cp );
        ProjectDomainNet ();
        UpdateDomains ();
        view_surf_1 = view_surf_2 = false;
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( true );
        xge_RedrawAll ();
        return 1;
    case sl11a1QUASIG2CONST:
        options1.q2coeff = xge_LogSlidebarValued ( Q2COEFF_MIN, Q2COEFF_MAX,
                                                   options1.slp );
        NotifyDoubleNumber ( win1, options1.q2coeff, false );
        g1h_DestroyQ2PrivateDatad ( domain1 );
        options1.spl_basis_valid = false;
        if ( options1.quasiG2 ) {
          InvalFinalSurface1 ();
          if ( view_surf_1 )
            view_surf_1 = swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( true );
          xge_RedrawAll ();
        }
        return 1;
    case sl11a2QUASIG2CONST:
        options2.q2coeff = xge_LogSlidebarValued ( Q2COEFF_MIN, Q2COEFF_MAX,
                                                   options2.slp );
        NotifyDoubleNumber ( win1, options2.q2coeff, false );
        g1h_DestroyQ2PrivateDatad ( domain2 );
        options2.spl_basis_valid = false;
        if ( options2.quasiG2 ) {
          InvalFinalSurface2 ();
          if ( view_surf_2 )
            view_surf_2 = swind_picture = false;
          if ( RenderingIsOn )
            BreakRendering ( true );
          xge_RedrawAll ();
        }
        return 1;
    default:
        break;
      }
      return 0;

  case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
    case intw11a1ORDER:
        if ( Option1SetOrder ( key ) ) {
          options1.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints1 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a1NK:
        if ( Option1SetNK ( key ) ) {
          options1.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints1 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a1M1:
        if ( Option1SetM1 ( key ) ) {
          options1.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints1 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a1M2:
        if ( Option1SetM2 ( key ) )
          options1.constr_matrix_valid = false;
          xge_RedrawAll ();
        return 1;
    case intw11a2ORDER:
        if ( Option2SetOrder ( key ) ) {
          options2.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints2 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a2NK:
        if ( Option2SetNK ( key ) ) {
          options2.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints2 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a2M1:
        if ( Option2SetM1 ( key ) ) {
          options2.constr_matrix_valid = false;
          ConfigureConstraintWidgets ( constraints2 );
          xge_RedrawAll ();
        }
        return 1;
    case intw11a2M2:
        if ( Option2SetM2 ( key ) )
          options2.constr_matrix_valid = false;
          xge_RedrawAll ();
        return 1;
    case intw11bHOLE_SIDES:
        return SetupHoleSides ( key );
    default:
        break;
      }
      return 0;

  case xgemsg_2DWIN_RESIZE:
  case xgemsg_2DWIN_PROJCHANGE:
      ProjectDomainNet ();
      xge_SetClipping ( domwind.er );
      domwind.er->redraw ( domwind.er, true );
      return 1;

  case xgemsg_2DWIN_PICK_POINT:
      return FindNearestDomainCPoint ( x, y, xge_MINDIST );

  case xgemsg_2DWIN_MOVE_POINT:
      SetDomainCPoint ( x, y );
      view_surf_1 = view_surf_2 = swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( true );
      xge_RedrawAll ();
      return 1;

  case xgemsg_2DWIN_SELECT_POINTS:
      SelectDomainCPoints ();
      xge_SetClipping ( domwind.er );
      domwind.er->redraw ( domwind.er, true );
      return 1;

  case xgemsg_2DWIN_UNSELECT_POINTS:
      UnselectDomainCPoints ();
      xge_SetClipping ( domwind.er );
      domwind.er->redraw ( domwind.er, true );
      return 1;

  case xgemsg_2DWIN_CHANGE_TRANS:
      Notify2DTransChange ( win1, er->data0 );
      return 1;

  case xgemsg_2DWIN_SAVE_POINTS:
      SaveDomainCPoints ();
      return 1;

  case xgemsg_2DWIN_TRANSFORM_POINTS:
      Notify2DTrans ( win1, er->data0 );
      TransformDomainCPoints ();
      if ( !view_surf_1 && !view_surf_2 ) {
        xge_SetClipping ( domwind.er );
        domwind.er->redraw ( domwind.er, true );
      }
      else {
        view_surf_1 = view_surf_2 = swind_picture = false;
      if ( RenderingIsOn )
        BreakRendering ( false );
        xge_RedrawAll ();
      }
      return 1;

  case xgemsg_2DWIN_FIND_REFBBOX:
      FindDomainBoundingBox ( &domwind.RefBBox );
      xge_2DwindInitProjection ( &domwind, domwind.RefBBox.x0, domwind.RefBBox.x1,
                                 domwind.RefBBox.y0, domwind.RefBBox.y1 );
      ProjectDomainNet ();
      xge_SetClipping ( domwind.er );
      domwind.er->redraw ( domwind.er, true );
      return 1;

  case xgemsg_2DWIN_ERROR:
      return 0;

  case xgemsg_KNOTWIN_CHANGE_KNOT:
      if ( er == knwind.er ) {
        UpdateDomains ();
        InitConstraintFrame ( 1 );
        InitConstraintFrame ( 2 );
        xge_SetClipping ( domwind.er );
        domwind.er->redraw ( domwind.er, true );
        swind_picture = false;
        if ( RenderingIsOn )
          BreakRendering ( false );
        xge_SetClipping ( swind.fww.er );
        swind.fww.er->redraw ( swind.fww.er, true );
        return 1;
      }
      else
        return 0;

  default:
      return 0;
    }
  }
  else {
    switch ( msg ) {
  case xgemsg_RESIZE:
      ResizeWinStatus ( win1 );
      if ( key )
        xge_Redraw ();
      return 1;

  case xgemsg_KEY:
      return ProcessKey ( key );

  case xgemsg_IDLE_COMMAND:
      ProcessIdleCommand ( key, x, y );
      return 1;

  default:
      return 0;
    }
  }
  return 0;
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

