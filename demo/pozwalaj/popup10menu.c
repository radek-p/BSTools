
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


xge_widget *InitPopup10 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewButton ( win1, NULL, btnP10BEZCURVE, 88, 19, 220+2, 50+22, txtBezierCurve );
  w = xge_NewButton ( win1, w, btnP10BEZPATCH, 88, 19, 220+2, 50+42, txtBezierPatch );
  w = xge_NewButton ( win1, w, btnP10BSCURVE, 88, 19, 220+2, 50+62, txtBSplineCurve );
  w = xge_NewButton ( win1, w, btnP10BSPATCH, 88, 19, 220+2, 50+82, txtBSplinePatch );
  w = xge_NewButton ( win1, w, btnP10BSMESH, 88, 19, 220+2, 50+102, txtBSplineMesh );
  w = xge_NewButton ( win1, w, btnP10BSHOLE, 88, 19, 220+2, 50+122, txtBSplineHole );
  if ( w ) w->state = xgestate_BUTTON_INACTIVE;
  menu = xge_NewFMenu ( win1, NULL, POPUP10, 92, 123, 220+0, 50+20, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup10*/

int Popup10CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP10BEZCURVE:
      xge_RemovePopups ( true );
      SetupPopup11Bezc ();
      OpenPopup ( popup11, true );
      return 1;
  case btnP10BEZPATCH:
      xge_RemovePopups ( true );
      SetupPopup11Bezp ();
      OpenPopup ( popup11, true );
      return 1;
  case btnP10BSCURVE:
      xge_RemovePopups ( true );
      SetupPopup11Bsc ();
      OpenPopup ( popup11, true );
      return 1;
  case btnP10BSPATCH:
      xge_RemovePopups ( true );
      SetupPopup11Bsp ();
      OpenPopup ( popup11, true );
      return 1;
  case btnP10BSMESH:
      xge_RemovePopups ( true );
      SetupPopup11Bsm ();
      OpenPopup ( popup11, true );
      return 1;
  case btnP10BSHOLE:
      xge_RemovePopups ( true );
      SetupPopup11Bsh ();
      OpenPopup ( popup11, true );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup10CallBack*/

