
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


xge_widget *InitSide00Menu ( xge_widget *prev )
{
  xge_widget *w, *menu;

  w = xge_NewSwitch ( win0, NULL, swM01MARK, 100, 16, 0, 22, txtMarkUnmark,
                      &swwin0mark );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_0, 16, 16,  0, 42, NULL, &markbits[0] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_1, 16, 16, 18, 42, NULL, &markbits[1] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_2, 16, 16, 36, 42, NULL, &markbits[2] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_3, 16, 16, 54, 42, NULL, &markbits[3] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_4, 42, 16, 72, 42, txtBit, &markbits[4] );
  w = xge_NewSwitch ( win0, w, swM01TRANSLATE, 100, 16, 0, 62, txtTranslate,
                      &swwin0translate );
  w = xge_NewSwitch ( win0, w, swM01SCALE, 100, 16, 0, 82, txtScale,
                      &swwin0scale );
  w = xge_NewSwitch ( win0, w, swM01ROTATE, 100, 16, 0, 102, txtRotate,
                      &swwin0rotate );
  w = xge_NewSwitch ( win0, w, swM01SHEAR, 100, 16, 0, 122, txtShear,
                      &swwin0shear );
  w = xge_NewSwitch ( win0, w, swM01PANZOOM, 100, 16, 0, 142, txtPanZoom,
                      &swwin0panzoom );
  w = xge_NewSwitch ( win0, w, swM01COORDINATES, 100, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swwin0coordinates );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00widgets = side00awidgets = w;

  w = xge_NewSwitch ( win0, NULL, swM01MARK, 100, 16, 0, 22, txtMarkUnmark,
                      &swwin0mark );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_0, 16, 16,  0, 42, NULL, &markbits[0] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_1, 16, 16, 18, 42, NULL, &markbits[1] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_2, 16, 16, 36, 42, NULL, &markbits[2] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_3, 16, 16, 54, 42, NULL, &markbits[3] );
  w = xge_NewSwitch ( win0, w, swM01MKBIT_4, 42, 16, 72, 42, txtBit, &markbits[4] );
  w = xge_NewSwitch ( win0, w, swM01TRANSLATE, 100, 16, 0, 62, txtTranslate,
                      &swwin0translate );
  w = xge_NewSwitch ( win0, w, swM01SCALE, 100, 16, 0, 82, txtScale,
                      &swwin0scale );
  w = xge_NewSwitch ( win0, w, swM01ROTATE, 100, 16, 0, 102, txtRotate,
                      &swwin0rotate );
  w = xge_NewSwitch ( win0, w, swM01SHEAR, 100, 16, 0, 122, txtShear,
                      &swwin0shear );
  w = xge_NewSwitch ( win0, w, swM01PANZOOM, 100, 16, 0, 142, txtPanZoom,
                      &swwin0panzoom );
  w = xge_NewSwitch ( win0, w, swM01SELECT_VERTEX, 100, 16, 0, 162,
                      txtSelectVertex, &sw_bsm_selectvertex );
  w = xge_NewSwitch ( win0, w, swM01SELECT_EDGE, 100, 16, 0, 182,
                      txtSelectEdge, &sw_bsm_selectedge );
  w = xge_NewSwitch ( win0, w, swM01COORDINATES, 100, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swwin0coordinates );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00dwidgets = w;

  InitSide00bMenu ();
  InitSide00eMenu ();
  menu = xge_NewMenu ( win0, prev, SIDEMENU0, SIDEMENUWIDTH0,
                       xge_HEIGHT-TOPMENUHEIGHT, 0, TOPMENUHEIGHT,
                       side00widgets );
    return menu;
} /*InitSide00Menu*/

boolean ChangeSide00MenuWidth ( short h )
{
  boolean wide, result;

  if ( side00menu->data1 == side00awidgets )
    wide = false;
  else if ( side00menu->data1 == side00bwidgets )
    wide = false;
  else if ( side00menu->data1 == side00cwidgets )
    wide = false;
  else if ( side00menu->data1 == side00dwidgets )
    wide = false;
  else if ( side00menu->data1 == side00ewidgets )
    wide = side00esw.contents->h > h-40;
  else if ( side00menu->data1 == side00fwidgets )
    wide = false;
  else
    wide = false;
    
  result = wide != side00menu_wide;
  side00menu_wide = wide;
  return result;
} /*ChangeSide00MenuWidth*/

void SetupWin0StatusLine ( void )
{
  if ( win0statusline )
    command0[0] = 0;
  ResizeWindow0 ( xge_current_width, xge_current_height );
} /*SetupWin0StatusLine*/

void SetupWin0CommandLine ( void )
{
  if ( win0commandline ) {
    command0[0] = 0;
    command0_editor.start = command0_editor.pos = 0;
    command0_editor.er->state = xgestate_NOTHING;
  }
  ResizeWindow0 ( xge_current_width, xge_current_height );
} /*SetupWin0CommandLine*/

void SelectGeomTool ( char tool, boolean *sw )
{
  boolean on;
  on = *sw;

  swwin0mark = swwin0translate = swwin0scale = swwin0rotate =
  swwin0shear = swwin0panzoom =
  sw_bsm_selectvertex = sw_bsm_selectedge = false;
  if ( on ) {
    if ( geom00win == geom00win2D )
      xge_2DwindEnableGeomWidget ( &g00win2D, tool );
    else if ( geom00win == geom00win3D )
      xge_3DwindEnableGeomWidget ( &g00win3D, tool );
    *sw = true;
  }
  else {
    if ( geom00win == geom00win2D )
      xge_2DwindEnableGeomWidget ( &g00win2D, xge_2DWIN_NO_TOOL );
    else if ( geom00win == geom00win3D )
      xge_3DwindEnableGeomWidget ( &g00win3D, xge_3DWIN_NO_TOOL );
  }
  xge_SetClipping ( side00menu );
  side00menu->redraw ( side00menu, true );
  xge_SetClipping ( geom00menu );
  geom00menu->redraw ( geom00menu, true );
} /*SelectGeomTool*/

void SetToolSwitches ( char tool, char coord )
{
  char sw_bsm;

  sw_bsm = sw_bsm_selectvertex;
  swwin0mark = swwin0translate = swwin0scale = swwin0rotate =
  swwin0shear = swwin0panzoom = sw_bsm_selectvertex =
  sw_bsm_selectedge = false;
  switch ( tool ) {
case xge_3DWIN_SELECTING_TOOL:
    swwin0mark = true;
    break;
case xge_3DWIN_SPECIAL_SELECTING_TOOL:
    sw_bsm_selectvertex = sw_bsm;
    sw_bsm_selectedge = !sw_bsm;
    break;
case xge_3DWIN_MOVING_TOOL:
    swwin0translate = true;
    break;
case xge_3DWIN_SCALING_TOOL:
    swwin0scale = true;
    break;
case xge_3DWIN_ROTATING_TOOL:
    swwin0rotate = true;
    break;
case xge_3DWIN_SHEAR_TOOL:
    swwin0shear = true;
    break;
case xge_3DWIN_PANNING_TOOL:
    swwin0panzoom = true;
    break;
default:
    break;
  }
  swwin0coordinates = coord;
} /*SetToolSwitches*/

int Side00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM01MKBIT_0:
  case swM01MKBIT_1:
  case swM01MKBIT_2:
  case swM01MKBIT_3:
  case swM01MKBIT_4:
      GeomObjectSetMarkingMask ( markbits );
      xge_SetClipping ( geom00menu );
      geom00menu->redraw ( geom00menu, true );
      return 1;
  case swM01MARK:
      SelectGeomTool ( xge_3DWIN_SELECTING_TOOL, &swwin0mark );
      return 1;
  case swM01TRANSLATE:
      SelectGeomTool ( xge_3DWIN_MOVING_TOOL, &swwin0translate );
      return 1;
  case swM01SCALE:
      SelectGeomTool ( xge_3DWIN_SCALING_TOOL, &swwin0scale );
      return 1;
  case swM01ROTATE:
      SelectGeomTool ( xge_3DWIN_ROTATING_TOOL, &swwin0rotate );
      return 1;
  case swM01SHEAR:
      SelectGeomTool ( xge_3DWIN_SHEAR_TOOL, &swwin0shear );
      return 1;
  case swM01PANZOOM:
      SelectGeomTool ( xge_3DWIN_PANNING_TOOL, &swwin0panzoom );
      return 1;
  case swM01SELECT_VERTEX:
      SelectGeomTool ( xge_3DWIN_SPECIAL_SELECTING_TOOL, &sw_bsm_selectvertex );
      return 1;
  case swM01SELECT_EDGE:
      SelectGeomTool ( xge_3DWIN_SPECIAL_SELECTING_TOOL, &sw_bsm_selectedge );
      return 1;
  case swM01COORDINATES:
      g00win2D.display_coord = g00win3D.display_coord = swwin0coordinates;
      return 1;
  case swM01STATUS:
      SetupWin0StatusLine ();
      return 1;
  case swM01COMMAND:
      SetupWin0CommandLine ();
      return 1;
    }

default:
    return 0;
  }
} /*Side00MenuCallBack*/

