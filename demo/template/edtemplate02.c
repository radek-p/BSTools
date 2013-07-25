
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
#include "xgedit.h"

#include "spltemplate.h"
#include "edtempwidgets.h"
#include "edtemplate.h"


int ProcessKey ( int key )
{
  switch ( key ) {
case 'm':
    xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
    return 1;
case 'M':
    xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
    return 1;
case 'q':  case 'Q':
    xge_done = 1;
    return 1;
default:
    return 0;
  }
} /*ProcessKey*/

int Win0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
    switch ( msg ) {
  case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
    case btn00FILE:
        xge_AddPopup ( popup00 );
        xge_GrabFocus ( popup00, true );
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
        xge_SetupFileList ( &filelist, ".", file_filter );
        xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
        xge_AddPopup ( popup01 );
        xge_GrabFocus ( popup01, true );
        return 1;
    case btnP00SAVE:
        xge_RemovePopup ( true );
        if ( FilenameCorrect ( filename ) ) {
          return 1;
        }
        else goto open_save_as;
    case btnP00SAVEAS:
        xge_RemovePopup ( true );
open_save_as:
        getcwd ( current_directory, MAX_PATH_LGT+1 );
        xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
        xge_AddPopup ( popup02 );
        xge_GrabFocus ( popup02, true );
        return 1;
    case btnP00EXIT:
        xge_done = 1;
        return 1;
    case btnP01OPEN:
        xge_RemovePopup ( true );
        OpenFile ( filename );
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
    default:
        return 0;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
    case sw01aMARK_UNMARK:
        if ( swind.selecting_mode ) {
          swind.panning = swind.moving_tool = swind.scaling_tool =
          swind.rotating_tool = false;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        }
        xge_3DwindResetGeomWidgets ( &swind );
        xge_Redraw ();
        return 1;
    case sw01aMOVE:
        if ( swind.moving_tool ) {
          swind.panning = swind.selecting_mode = swind.scaling_tool =
          swind.rotating_tool = false;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_MOVING_TOOL );
        }
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_3DwindResetGeomWidgets ( &swind );
        xge_Redraw ();
        return 1;
    case sw01aSCALE:
        if ( swind.scaling_tool ) {
          swind.panning = swind.selecting_mode = swind.moving_tool =
          swind.rotating_tool = false;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_SCALING_TOOL );
        }
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_3DwindResetGeomWidgets ( &swind );
        xge_Redraw ();
        return 1;
    case sw01aROTATE:
        if ( swind.rotating_tool ) {
          swind.panning = swind.selecting_mode = swind.moving_tool =
          swind.scaling_tool = false;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_ROTATING_TOOL );
        }
        else
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        xge_3DwindResetGeomWidgets ( &swind );
        xge_Redraw ();
        return 1;
    case sw01PAN_ZOOM:
        if ( swind.panning ) {
          swind.selecting_mode = swind.moving_tool = swind.scaling_tool =
          swind.rotating_tool = false;
          xge_3DwindEnableGeomWidget ( &swind, xge_3DWIN_NO_TOOL );
        }
        xge_3DwindResetGeomWidgets ( &swind );
        xge_Redraw ();
        return 1;
    case sw01COORDINATES:
        return 1;
    case sw02STATUS:
        StatusLineOnOff ( win0 );
        return 1;
    default:
        return 0;
      }

  case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
    case lbP01DIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupFileList ( &filelist, ".", file_filter );
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup01 );
          popup01->redraw ( popup01, true );
        }
        break;
    case lbP01FILELIST:
        break;
    case lbP02DIRLIST:
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]]) ) {
          xge_SetupDirList ( &dirlist, ".", NULL, current_directory );
          getcwd ( current_directory, MAX_PATH_LGT+1 );
          xge_SetClipping ( popup02 );
          popup02->redraw ( popup02, true );
        }
        break;
    default:
        break;
      }
      return 1;

  case xgemsg_3DWIN_RESIZE:
      return 0;

  case xgemsg_3DWIN_PROJCHANGE:
      return 0;

  case xgemsg_3DWIN_PICK_POINT:
      return 0;

  case xgemsg_3DWIN_MOVE_POINT:
      return 0;

  case xgemsg_3DWIN_SELECT_POINTS:
      return 0;

  case xgemsg_3DWIN_UNSELECT_POINTS:
      return 0;

  case xgemsg_3DWIN_CHANGE_TRANS:
      Notify3DTransChange ( win0, er->data1 );
      return 1;

  case xgemsg_3DWIN_SAVE_POINTS:
      return 0;

  case xgemsg_3DWIN_TRANSFORM_POINTS:
      Notify3DTrans ( win0, er->data1 );
      return 1;

  case xgemsg_3DWIN_FIND_REFBBOX:
      return 0;

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
      if ( key )
        xge_Redraw ();
      return 1;

  case xgemsg_KEY:
      return ProcessKey ( key );

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
    default:
        return 0;
      }
      break;

  case xgemsg_SWITCH_COMMAND:
      switch ( er->id ) {
    case sw11aMARK_UNMARK:
      if ( domwind.selecting_mode ) {
        domwind.panning = domwind.moving_tool = domwind.scaling_tool =
        domwind.rotating_tool = false;
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
      }
      return 1;
    case sw11aMOVE:
      if ( domwind.moving_tool ) {
        domwind.panning = domwind.selecting_mode = domwind.scaling_tool =
        domwind.rotating_tool = false;
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_MOVING_TOOL );
      }
      else
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
      return 1;
    case sw11aSCALE:
      if ( domwind.scaling_tool ) {
        domwind.panning = domwind.selecting_mode = domwind.moving_tool =
        domwind.rotating_tool = false;
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_SCALING_TOOL );
      }
      else
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
      return 1;
    case sw11aROTATE:
      if ( domwind.rotating_tool ) {
        domwind.panning = domwind.selecting_mode = domwind.moving_tool =
        domwind.scaling_tool = false;
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_ROTATING_TOOL );
      }
      else
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
      return 1;
    case sw11PAN_ZOOM:
      if ( domwind.panning ) {
        domwind.selecting_mode = domwind.moving_tool =
        domwind.scaling_tool = domwind.rotating_tool = false;
        xge_2DwindEnableGeomWidget ( &domwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
      }
      return 1;
    case sw11COORDINATES:
      return 1;
    case sw15STATUS:
        StatusLineOnOff ( win1 );
        return 1;
    default:
        return 0;
      }

  case xgemsg_2DWIN_RESIZE:
      return 0;

  case xgemsg_2DWIN_PROJCHANGE:
      return 0;

  case xgemsg_2DWIN_PICK_POINT:
      return 0;

  case xgemsg_2DWIN_MOVE_POINT:
      return 0;

  case xgemsg_2DWIN_SELECT_POINTS:
      return 0;

  case xgemsg_2DWIN_UNSELECT_POINTS:
      return 0;

  case xgemsg_2DWIN_CHANGE_TRANS:
      Notify2DTransChange ( win1, er->data0 );
      return 1;

  case xgemsg_2DWIN_SAVE_POINTS:
      return 0;

  case xgemsg_2DWIN_TRANSFORM_POINTS:
      Notify2DTrans ( win1, er->data0 );
      return 1;

  case xgemsg_2DWIN_FIND_REFBBOX:
      return 0;

  case xgemsg_2DWIN_ERROR:
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

