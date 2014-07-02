
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


static boolean sw2D = false, sw3D = true, swRational = false;


xge_widget *InitPopup11 ( void )
{
  xge_widget *w, *menu;

        /* this popup menu has 6 alternative widget lists,       */
        /* specific for the objects to be added. Before          */
        /* opening it, the appropriate list has to be assigned   */
        /* to the popup menu, and the address of the appropriate */
        /* widget (the first in each list) has to be assigned    */
        /* to the object name editor                             */
          /* widgets for a Bezier curve */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBezierCurve );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BEZC, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_2D, 60, 16, 20+50, 40+64, txt2D, &sw2D );
  w = xge_NewSwitch ( win1, w, swP11_3D, 60, 16, 20+168, 40+64, txt3D, &sw3D );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bezc =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );
          /* widgets for a Bezier patch */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBezierPatch );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BEZP, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_2D, 60, 16, 20+50, 40+64, txt2D, &sw2D );
  w = xge_NewSwitch ( win1, w, swP11_3D, 60, 16, 20+168, 40+64, txt3D, &sw3D );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bezp =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );
          /* widgets for a B-spline curve */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBSplineCurve );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BSC, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_2D, 60, 16, 20+50, 40+64, txt2D, &sw2D );
  w = xge_NewSwitch ( win1, w, swP11_3D, 60, 16, 20+168, 40+64, txt3D, &sw3D );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bsc =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );
          /* widgets for a B-spline patch */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBSplinePatch );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BSP, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_2D, 60, 16, 20+50, 40+64, txt2D, &sw2D );
  w = xge_NewSwitch ( win1, w, swP11_3D, 60, 16, 20+168, 40+64, txt3D, &sw3D );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bsp =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );
          /* widgets for a B-spline mesh */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBSplineMesh );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BSM, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_2D, 60, 16, 20+50, 40+64, txt2D, &sw2D );
  w = xge_NewSwitch ( win1, w, swP11_3D, 60, 16, 20+168, 40+64, txt3D, &sw3D );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bsm =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );
          /* widgets for a B-spline hole */  
  w = xge_NewTextWidget ( win1, NULL, 0, 18, 16, 20+10, 40+12, txtNew );
  w = xge_NewTextWidget ( win1, w, 0, 88, 16, 20+34, 40+12, txtBSplineHole );
  w = xge_NewTextWidget ( win1, w, 0, 40, 16, 20+10, 40+36, txtName );
  w = xge_NewButton ( win1, w, btnP11ADD_BSM, 58, 19, 20+78, 40+130-30, txtAdd );
  w = xge_NewButton ( win1, w, btnP11CANCEL, 58, 19, 20+224, 40+130-32, txtCancel );
  w = xge_NewSwitch ( win1, w, swP11_RATIONAL, 66, 16, 20+286, 40+64,
                      txtRational, &swRational );
  popup11wdg_bsh =
  w = xge_NewStringEd ( win1, w, txtedP11OBJNAME, 300, 19, 20+10+40, 40+36,
                        MAX_NAME_LENGTH, objectname, &objname_editor1 );

  menu = xge_NewFMenu ( win1, NULL, POPUP11, 360, 130, 20, 40, NULL );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup11*/

void SetupPopup11Bezc ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bezc, false );
  objname_editor1.er = popup11wdg_bezc;
} /*SetupPopup11Bezc*/

void SetupPopup11Bezp ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bezp, false );
  objname_editor1.er = popup11wdg_bezp;
} /*SetupPopup11Bezp*/

void SetupPopup11Bsc ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bsc, false );
  objname_editor1.er = popup11wdg_bsc;
} /*SetupPopup11Bsc*/

void SetupPopup11Bsp ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bsp, false );
  objname_editor1.er = popup11wdg_bsp;
} /*SetupPopup11Bsp*/

void SetupPopup11Bsm ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bsm, false );
  objname_editor1.er = popup11wdg_bsm;
} /*SetupPopup11Bsm*/

void SetupPopup11Bsh ( void )
{
  xge_SetMenuWidgets ( popup11, popup11wdg_bsh, false );
  objname_editor1.er = popup11wdg_bsh;
} /*SetupPopup11Bsh*/

void AddANewObject ( int obj_type )
{
  geom_object *go;
  char        dim;
  Box3d       box;

  if ( sw2D )
    dim = 2;
  else
    dim = 3;
  switch ( obj_type ) {
case GO_BEZIER_CURVE:
    go = GeomObjectAddBezierCurve ( objectname, dim, swRational );
    break;
case GO_BEZIER_PATCH:
    go = GeomObjectAddBezierPatch ( objectname, dim, swRational );
    break;
case GO_BSPLINE_CURVE:
    go = GeomObjectAddBSplineCurve ( objectname, dim, swRational );
    break;
case GO_BSPLINE_PATCH:
    go = GeomObjectAddBSplinePatch ( objectname, dim, swRational );
    break;
case GO_BSPLINE_MESH:
    go = GeomObjectAddBSplineMesh ( objectname, dim, swRational );
    break;
case GO_BSPLINE_HOLE:
    go = GeomObjectAddBSplineHole ( objectname, swRational );
    break;
default:
    return;
  }
  if ( go ) {
    objectname[0] = 0;
    objname_editor1.start = objname_editor1.pos = 0;
    SetupObjectSpecificMenus ( go );
    if ( geom00win == geom00win2D ) {
      if ( GeomObjectFindBoundingBox ( 2, &box ) ) {
        g00win2D.RefBBox.x0 = box.x0;  g00win2D.RefBBox.x1 = box.x1;
        g00win2D.RefBBox.y0 = box.y0;  g00win2D.RefBBox.y1 = box.y1;
        xge_2DwindSetupProjection ( &g00win2D );
        CallBack ( geom00win, xgemsg_2DWIN_PROJCHANGE, 0, 0, 0 );
      }
    }
    else {
      if ( GeomObjectFindBoundingBox ( 3, &g00win3D.RefBBox ) ) {
        xge_3DwindSetupParProj ( &g00win3D, &g00win3D.RefBBox );
        memcpy ( &g00win3D.PerspBBox, &g00win3D.RefBBox, sizeof(Box3d) );
        xge_3DwindSetupPerspProj ( &g00win3D, false );
      }
    }
    GeomObjectDisplayInfoText ( go );
    xge_RedrawAll ();
  }
} /*AddANewObject*/

int Popup11CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP11ADD_BEZC:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BEZIER_CURVE );
      return 1;
  case btnP11ADD_BEZP:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BEZIER_PATCH );
      return 1;
  case btnP11ADD_BSC:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BSPLINE_CURVE );
      return 1;
  case btnP11ADD_BSP:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BSPLINE_PATCH );
      return 1;
  case btnP11ADD_BSM:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BSPLINE_MESH );
      return 1;
  case btnP11ADD_BSH:
      xge_RemovePopup ( true );
      AddANewObject ( GO_BSPLINE_HOLE );
      return 1;
  case btnP11CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swP11_2D:
      sw3D = !sw2D;
      xge_SetClipping ( popup11 );
      popup11->redraw ( popup11, true );
      return 1;
  case swP11_3D:
      sw2D = !sw3D;
      xge_SetClipping ( popup11 );
      popup11->redraw ( popup11, true );
      return 1;
  case swP11_RATIONAL:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case txtedP11OBJNAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case txtedP11OBJNAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup11CallBack*/

