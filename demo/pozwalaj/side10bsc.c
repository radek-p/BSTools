
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


void InitSide10Menu_BSc ( void )
{
  xge_widget *w;

        /* widgets specific for B-spline curves */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBSplineCurve );
  w = xge_NewStringEd ( win1, w, textedM1BSC_NAME, 109, 19, 0, 40,
                        MAX_NAME_LENGTH, objectname, &bsc_name_ed );
  w = xge_NewIntWidget ( win1, w, intwM1BSC_DEG, 109, 19, 0, 60,
                         0, MAX_DEGREE+1, &intw_bscdeg, txtDegree, &degree );
  w = xge_NewButton ( win1, w, btnM1BSC_UNIFORM, 76, 19, 0, 80, txtUniform );
  w = xge_NewButton ( win1, w, btnM1BSC_REFINE, 76, 19, 0, 100, txtRefine );
  w = xge_NewSwitch ( win1, w, swM1BSC_CLOSED, 109, 16, 0, 120,
                      txtClosed, &bsc_sw_closed );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, btnM1BSC_VIEW_CURVE, 109, 16, 0, 22,
                      txtCurve, &sw_view_curve );
  w = xge_NewSwitch ( win1, w, btnM1BSC_VIEW_CPOLY, 109, 16, 0, 42,
                      txtControlPolygon, &sw_view_cpoly );
  w = xge_NewButton ( win1, w, btnM1BSC_COLOUR, 60, 18, 0, 62, txtColour );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_view = w;
          /* data */
  w = xge_NewSwitch ( win1, NULL, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bsc_data = w;

} /*InitSide10Menu_BSc*/

void SetupBSplineCurveWidgets ( GO_BSplineCurve *obj )
{
  InitNameEditor ( &bsc_name_ed, obj->me.name );
  degree = obj->degree;
  sw_view_curve = obj->view_curve;
  sw_view_cpoly = obj->view_cpoly;
  bsc_sw_closed = obj->closed;
  Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                        obj->knots, obj->closed );
  SetGeomWin10Knotw ( obj );
} /*SetupBSplineCurveWidgets*/

int Side10MenuBscCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BSplineCurve *obj;

  if ( current_go->obj_type != GO_BSPLINE_CURVE )
    return 0;
  obj = (GO_BSplineCurve*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BSC_UNIFORM:
      if ( GeomObjectBSplineCurveSetUniformKnots ( obj, true ) ) {
        xge_KnotWindFindMapping ( &g10knotwin );
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSC_REFINE:
      if ( GeomObjectBSplineCurveRefine ( obj ) ) {
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_KnotWindFindMapping ( &g10knotwin );
        xge_RedrawAll ();
      }
      return 1;
  case btnM1BSC_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case btnM1BSC_VIEW_CURVE:
      obj->view_curve = sw_view_curve;
      RedrawGeom00Win ();
      return 1;
  case btnM1BSC_VIEW_CPOLY:
      obj->view_cpoly = sw_view_cpoly;
      RedrawGeom00Win ();
      return 1;
  case swM1BSC_CLOSED:
      if ( GeomObjectBSplineCurveSetClosed ( obj, bsc_sw_closed ) ) {
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_RedrawAll ();
      }
      else if ( bsc_sw_closed )
        xge_DisplayErrorMessage ( ErrorMsgCannotClose, 0 );
      bsc_sw_closed = obj->closed;
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BSC_DEG:
      if ( GeomObjectBSplineCurveSetDegree ( obj, key ) ) {
        degree = obj->degree;
        Geom10winKNSetKnots ( &g10knotwin, obj->degree, obj->lastknot,
                              obj->knots, obj->closed );
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BSC_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BSC_NAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side10MenuBscCallBack*/

