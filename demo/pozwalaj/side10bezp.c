
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


void InitSide10Menu_Bezp ( void )
{
  xge_widget *w;

        /* widgets specific for Bezier patches */
          /* edit */
  w = xge_NewTextWidget ( win1, NULL, 0, 109, 19, 0, 20, txtBezierPatch );
  w = xge_NewStringEd ( win1, w, textedM1BEZP_NAME, 109, 19, 0, 40,
                        MAX_NAME_LENGTH, objectname, &bezp_name_ed );
  w = xge_NewIntWidget ( win1, w, intwM1BEZP_DEGU, 76, 19, 0, 60,
                         0, MAX_DEGREE+1, &intw_bezpdegu, txtDegreeU, &degreeu );
  w = xge_NewIntWidget ( win1, w, intwM1BEZP_DEGV, 76, 19, 0, 80,
                         0, MAX_DEGREE+1, &intw_bezpdegv, txtDegreeV, &degreev );
  w = xge_NewButton ( win1, w, btnM1BEZP_FLIP, 76, 19, 0, 100, txtFlip );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezp_edit = w;
          /* view */
  w = xge_NewSwitch ( win1, NULL, swM1BEZP_VIEW_SURF, 109, 16, 0, 22,
                      txtSurface, &sw_view_surf );
  w = xge_NewSwitch ( win1, w, swM1BEZP_VIEW_CNET, 109, 16, 0, 42,
                      txtControlNet, &sw_view_cnet );
  w = xge_NewIntWidget ( win1, w, intwM1BEZP_DENSITY_U, 109, 19, 0, 62,
                         1, MAX_PNET_DENSITY, &intwdensityu, txtDensityU,
                         &density_u );
  w = xge_NewIntWidget ( win1, w, intwM1BEZP_DENSITY_V, 109, 19, 0, 82,
                         1, MAX_PNET_DENSITY, &intwdensityv, txtDensityV,
                         &density_v );
  w = xge_NewButton ( win1, w, btnM1BEZP_COLOUR, 60, 18, 0, 102, txtColour );
  w = xge_NewSwitch ( win1, w, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezp_view = w;
          /* data */
  w = xge_NewSwitch ( win1, NULL, swM11STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win1statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win1, w, swM11COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win1commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side10wdg_bezp_data = w;
} /*InitSide10Menu_Bezp*/

boolean ChangeSide10MenuWidth_Bezp ( short h )
{
  boolean result;

  result = side10menu_wide;
  side10menu_wide = false;
  return result;
} /*ChangeSide10MenuWidth_Bezp*/

void SetupBezierPatchWidgets ( GO_BezierPatch *obj )
{
  InitNameEditor ( &bezp_name_ed, obj->me.name );
  degreeu = obj->degree_u;
  degreev = obj->degree_v;
  intwdensityu.title = &txtDensityU[0];
  intwdensityv.title = &txtDensityV[0];
  density_u = obj->dens_u;
  density_v = obj->dens_v;
  sw_view_surf = obj->view_surf;
  sw_view_cnet = obj->view_cnet;
  SetGeomWin10Empty ();
} /*SetupBezierPatchWidgets*/

int Side10MenuBezpCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  GO_BezierPatch *obj;

  if ( current_go->obj_type != GO_BEZIER_PATCH )
    return 0;
  obj = (GO_BezierPatch*)current_go;
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM1BEZP_FLIP:
      if ( GeomObjectBezierPatchFlipUV ( obj ) ) {
        SetupBezierPatchWidgets ( obj );
        xge_Redraw ();
      }
      return 1;
  case btnM1BEZP_COLOUR:
      memcpy ( colour_rgb, obj->me.colour, 3*sizeof(double) );
      OpenPopup ( popup13, false );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM1BEZP_VIEW_SURF:
      obj->view_surf = sw_view_surf;
      RedrawGeom00Win ();
      return 1;
  case swM1BEZP_VIEW_CNET:
      obj->view_cnet = sw_view_cnet;
      RedrawGeom00Win ();
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM1BEZP_DEGU:
      if ( GeomObjectBezierPatchSetDegreeU ( obj, key ) ) {
        degreeu = obj->degree_u;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BEZP_DEGV:
      if ( GeomObjectBezierPatchSetDegreeV ( obj, key ) ) {
        degreev = obj->degree_v;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BEZP_DENSITY_U:
      if ( GeomObjectBezierPatchSetDensityU ( obj, key ) ) {
        density_u = obj->dens_u;
        xge_RedrawAll ();
      }
      return 1;
  case intwM1BEZP_DENSITY_V:
      if ( GeomObjectBezierPatchSetDensityV ( obj, key ) ) {
        density_v = obj->dens_v;
        xge_RedrawAll ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case textedM1BEZP_NAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM1BEZP_NAME:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side10MenuBezpCallBack*/

