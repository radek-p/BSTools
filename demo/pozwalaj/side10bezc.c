
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


void InitSide10Menu_Bezc ( void )
{
  xge_widget *w;

        /* widgets specific for Bezier curves */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBezierCurve );
  w = xge_NewStringEd ( win1, w, textedM1BEZC_NAME, 109, 19, 0, 40,
                        MAX_NAME_LENGTH, objectname, &bezc_name_ed );
  w = xge_NewIntWidget ( win1, w, intwM1BEZC_DEG, 109, 19, 0, 60,
                         0, MAX_DEGREE+1, &intw_bezcdeg, txtDegree, &degree );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezc = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, btnM1BEZC_VIEW_CURVE, 109, 16, 0, 22,
                      txtCurve, &sw_view_curve );
  w = xge_NewSwitch ( win1, w, btnM1BEZC_VIEW_CPOLY, 109, 16, 0, 42,
                      txtControlPolygon, &sw_view_cpoly );
  w = xge_NewButton ( win1, w, btnM1BEZC_COLOUR, 60, 18, 0, 62, txtColour );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezc_view = w;
          /* data */
  w = xge_NewSwitch ( win1, NULL, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezc_data = w;
} /*InitSide10Menu_Bezc*/

void SetupBezierCurveWidgets ( GO_BezierCurve *obj )
{
  InitNameEditor ( &bezc_name_ed, obj->me.name );
  degree = obj->degree;
  sw_view_curve = obj->view_curve;
  sw_view_cpoly = obj->view_cpoly;
  SetGeomWin10Empty ();
} /*SetupBezierCurveWidgets*/

int Side10MenuBezcCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BezierCurve *obj;

  if ( current_go->obj_type != GO_BEZIER_CURVE )
    return 0;
  obj = (GO_BezierCurve*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BEZC_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case btnM1BEZC_VIEW_CURVE:
      obj->view_curve = sw_view_curve;
      RedrawGeom00Win ();
      return 1;
  case btnM1BEZC_VIEW_CPOLY:
      obj->view_cpoly = sw_view_cpoly;
      RedrawGeom00Win ();
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BEZC_DEG:
      if ( GeomObjectBezierCurveSetDegree ( obj, key ) ) {
        degree = obj->degree;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BEZC_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BEZC_NAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side10MenuBezcCallBack*/

