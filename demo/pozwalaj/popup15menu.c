
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
#include "pkrender.h"
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


xge_widget *InitPopup15 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewButton ( win1, NULL, btnP15GENERAL, 76, 19, 0,  79, txtGeneral );
  if ( w ) w->state = xgestate_BUTTON_COMBO_1;
  w = xge_NewButton ( win1, w, btnP15SWEPT,      76, 19, 0,  98, txtSwept );
  if ( w ) w->state = /*xgestate_BUTTON_COMBO_1*/ xgestate_BUTTON_INACTIVE;
  w = xge_NewButton ( win1, w, btnP15SPHERICAL,  76, 19, 0, 117, txtSpherical );
  if ( w ) w->state = xgestate_BUTTON_COMBO_1;
  w = xge_NewButton ( win1, w, btnP15LOFTED,     76, 19, 0, 136, txtLofted );
  if ( w ) w->state = /*xgestate_BUTTON_COMBO_1*/ xgestate_BUTTON_INACTIVE;
  w = xge_NewButton ( win1, w, btnP15BLENDINGG1, 76, 19, 0, 155, txtBlendingG1 );
  if ( w ) w->state = xgestate_BUTTON_COMBO_1;
  w = xge_NewButton ( win1, w, btnP15BLENDINGG2, 76, 19, 0, 174, txtBlendingG2 );
  if ( w ) w->state = xgestate_BUTTON_COMBO_1;
  menu = xge_NewFMenu ( win1, NULL, POPUP15,     76, 115, 0, 78, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup15*/

int Popup15CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplinePatch *obj;

  if ( current_go->obj_type != GO_BSPLINE_PATCH )
    return 0;
  obj = (GO_BSplinePatch*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP15GENERAL:
      xge_RemovePopup ( false );
      BSPatchSetupGeneralOptions ( obj );
      return 1;
  case btnP15SWEPT:
      xge_RemovePopup ( true );
/*      bsp_type_button->data0 = txtSwept; */
      xge_DisplayErrorMessage ( ErrorMsgNotImplemented, 0 );
      return 1;
  case btnP15SPHERICAL:
      xge_RemovePopup ( true );
      BSPatchSetupSphericalOptions ( obj );
      return 1;
  case btnP15LOFTED:
      xge_RemovePopup ( true );
/*      bsp_type_button->data0 = txtLofted; */
      xge_DisplayErrorMessage ( ErrorMsgNotImplemented, 0 );
      return 1;
  case btnP15BLENDINGG1:
      xge_RemovePopup ( false );
      BSPatchSetupBlG1Options ( obj );
      return 1;
  case btnP15BLENDINGG2:
      xge_RemovePopup ( false );
      BSPatchSetupBlG2Options ( obj );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup14CallBack*/

