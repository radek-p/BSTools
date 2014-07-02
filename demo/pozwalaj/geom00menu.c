
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
#include "render.h"
#include "edlight.h"


xge_widget *InitGeom00Menu ( xge_widget *prev )
{
  xge_widget *menu;

  geom00win3D = xge_New3Dwind ( win0, NULL, GEOMWIN3D0,
                                xge_WIDTH-SIDEMENUWIDTH0,
                                xge_HEIGHT-TOPMENUHEIGHT,
                                SIDEMENUWIDTH0, TOPMENUHEIGHT,
                                &g00win3D, DrawG00win3Dpar, DrawG00win3Dpersp );

  geom00win2D = xge_New2Dwind ( win0, NULL, GEOMWIN2D0,
                                xge_WIDTH-SIDEMENUWIDTH0,
                                xge_HEIGHT-TOPMENUHEIGHT,
                                SIDEMENUWIDTH0, TOPMENUHEIGHT,
                                &g00win2D, DrawG00win2D );

  geom00win = geom00win3D;
  menu = xge_NewMenu ( win0, prev, GEOMMENU0, xge_WIDTH-SIDEMENUWIDTH0,
                       xge_HEIGHT-TOPMENUHEIGHT, SIDEMENUWIDTH0, TOPMENUHEIGHT,
                       geom00win );
  return menu;
} /*InitGeom00Menu*/

void DrawG00win3Dpar ( xge_widget *er, boolean onscreen )
{
  int id;

  id = er->id & 0x03;  /* it must be 0, 1 or 2 */
  glXWaitX ();
  xgle_DrawGeomWinBackground ( er, 0 );
  xgle_SetIdentMapping ( er );
  xgle_3DwindDrawCursorPos ( &g00win3D, id, xge_xx, xge_yy );
  xgle_SetGLCamerad ( &g00win3D.CPos[id] );
  glClear ( GL_DEPTH_BUFFER_BIT );
  glEnable ( GL_DEPTH_TEST );
  DisplaySpecial3DElements ( &g00win3D );
  GeomObjectDisplayActive ( 3 );
  glDisable ( GL_DEPTH_TEST );
  glXWaitGL ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );
  xge_DrawGeomWinSelectionRect ( er, &g00win3D.selection_rect );
  xge_3DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawG00win3Dpar*/

static void FindFogRange ( CameraRecd *CPos, Box3d *bbox, float *zmin, float *zmax )
{
  point3d p, q;
  float   d, e;

  SetPoint3d ( &p, 0.5*(bbox->x0+bbox->x1), 0.5*(bbox->y0+bbox->y1),
                   0.5*(bbox->z0+bbox->z1 ) );
  CameraProjectPoint3d ( CPos, &p, &q );
  d = bbox->x1-bbox->x0;
  e = bbox->y1-bbox->y0;
  d = max ( d, e );
  e = bbox->z1-bbox->z0;
  d = max ( d, e );
  *zmin = q.z-0.25*d;
  *zmax = q.z+1.75*d;
} /*FindFogRange*/

void DrawG00win3Dpersp ( xge_widget *er, boolean onscreen )
{
  int   id;
  float zmin, zmax;

  id = er->id & 0x03;  /* it must be 3 */
  if ( rendered_picture ) {
    XPutImage ( xgedisplay, xgepixmap, xgegc, rendimage,
                er->x, er->y, er->x, er->y, er->w, er->h );
    if ( onscreen )
      xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
  }
  else {
    glXWaitX ();
    xgle_DrawGeomWinBackground ( er, 0 );
    xgle_SetIdentMapping ( er );
    xgle_3DwindDrawCursorPos ( &g00win3D, id, xge_xx, xge_yy );
    xgle_SetGLCamerad ( &g00win3D.CPos[id] );
    glClear ( GL_DEPTH_BUFFER_BIT );
    glEnable ( GL_DEPTH_TEST );

    glEnable ( GL_FOG );
    glFogi ( GL_FOG_MODE, GL_LINEAR );
    FindFogRange ( &g00win3D.CPos[id], &g00win3D.PerspBBox, &zmin, &zmax );
    glFogf ( GL_FOG_START, zmin );
    glFogf ( GL_FOG_END, zmax );
    glHint ( GL_FOG_HINT, GL_DONT_CARE );
    glFogfv ( GL_FOG_COLOR, xglec_Blue6 );
    
    DisplaySpecial3DElements ( &g00win3D );
    GeomObjectDisplayActive ( 3 );
    glDisable ( GL_FOG );
    glDisable ( GL_DEPTH_TEST );
    glXWaitGL ();
    xgleCopyGLRect ( er->w, er->h, er->x, er->y );
    xge_DrawGeomWinSelectionRect ( er, &g00win3D.selection_rect );
    xge_DrawGeomWinFrame ( er, onscreen );
  }
} /*DrawG00win3Dpersp*/

void Geom00Win3DShowTransformation ( xge_3Dwind *ww )
{
  char text[256], txt1[32];
  vector3d *sf;

  text[0] = 0;
  switch ( ww->current_tool ) {
case xge_3DWIN_MOVING_TOOL:
    sf = &ww->trans_params;
    sprintf ( text, "moving: %7.3f, %7.3f, %7.3f", sf->x, sf->y, sf->z );
    break;
case xge_3DWIN_SCALING_TOOL:
    sf = &ww->scaling_factors;
    sprintf ( text, "scaling: %7.3f, %7.3f, %7.3f", sf->x, sf->y, sf->z );
    break;
case xge_3DWIN_ROTATING_TOOL:
    pkv_RadToDegreeStr ( ww->trans_params.x, txt1 );
    sprintf ( text, "rotating: %s (%8.5f)", txt1, ww->trans_params.x );
    break;
case xge_3DWIN_SHEAR_TOOL:
    switch ( ww->tool_mode ) {
  case 1:
  case 2:
      sprintf ( text, "shearing x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[0].x,
                ww->shear_axis[0].y, ww->shear_axis[0].z );
      break;
  case 3:
  case 4:
      sprintf ( text, "shearing x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[1].x,
                ww->shear_axis[1].y, ww->shear_axis[1].z );
      break;
  case 5:
  case 6:
      sprintf ( text, "shearing x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[2].x,
                ww->shear_axis[2].y, ww->shear_axis[2].z );
      break;
  default:
      break;
    }
    break;
default:
    break;
  }
  SetStatusText ( text, true );
} /*Geom00Win3DShowTransformation*/

void Geom00Win3DShowTransOrigin ( xge_3Dwind *ww )
{
  char text[256];

  text[0] = 0;
  switch ( ww->current_tool ) {
case xge_3DWIN_MOVING_TOOL:
    break;
case xge_3DWIN_SCALING_TOOL:
    sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->scaling_centre.x,
              ww->scaling_centre.y, ww->scaling_centre.z );
    break;
case xge_3DWIN_ROTATING_TOOL:
    sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->rotating_centre.x,
              ww->rotating_centre.y, ww->rotating_centre.z );
    break;
case xge_3DWIN_SHEAR_TOOL:
    switch ( ww->tool_mode ) {
  case 0:
      sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_centre.x,
                ww->shear_centre.y, ww->shear_centre.z );
      break;
  case 1:
  case 2:
      sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[0].x,
                ww->shear_axis[0].y, ww->shear_axis[0].z );
      break;
  case 3:
  case 4:
      sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[1].x,
                ww->shear_axis[1].y, ww->shear_axis[1].z );
      break;
  case 5:
  case 6:
      sprintf ( text, "x: %7.3f, y: %7.3f, z: %7.3f", ww->shear_axis[2].x,
                ww->shear_axis[2].y, ww->shear_axis[2].z );
      break;
  default:
      break;
    }
    break;
default:
    break;
  }
  SetStatusText ( text, true );
} /*Geom00Win3DShowTransOrigin*/

void Geom00Win3DShowCameraPos ( CameraRecd *CPos )
{
  char s[MAX_COMMAND_LGT];

  sprintf ( s,
    "Camera: w = %d, h = %d, pos = (%f,%f,%f), psi = %f, theta = %f, phi = %f, f = %f",
    CPos->width, CPos->height, CPos->position.x, CPos->position.y, CPos->position.z,
    CPos->psi, CPos->theta, CPos->phi, CPos->vd.persp.f );
  SetStatusText ( s, true );
} /*Geom00Win3DShowCameraPos*/

void Geom00WinRestartRendering ( void )
{
  if ( rendered_picture ) {
    rendered_picture = false;
    if ( g00win3D.fww.zoomwin == -1 || g00win3D.fww.zoomwin == 3 ) {
      RendEnterCamerad ( &g00win3D.CPos[3], g00win3D.fww.win[3] );
      xge_SetClipping ( g00win3D.fww.win[3] );
      g00win3D.fww.win[3]->redraw ( g00win3D.fww.win[3], false );
      XGetSubImage ( xgedisplay, xgepixmap, g00win3D.cwin[3]->x, g00win3D.cwin[3]->y,
                     g00win3D.cwin[3]->w, g00win3D.cwin[3]->h, 0xFFFFFFFF, ZPixmap,
                     rendimage, g00win3D.cwin[3]->x, g00win3D.cwin[3]->y );
      renderbtn0->data0 = renderbtn1->data0 = renderbtn2->data0 = txtInterrupt;
      if ( !RenderingIsOn ) {
        xge_PostIdleCommand ( IDLE_COMMAND_RENDER_CONT, 0, 0 );
        RenderingIsOn = true;
      }
      RendRestart ();
      rendered_picture = true;
    }
  }
} /*Geom00WinRestartRendering*/

int Geom00Win3DCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int        id;
  xge_3Dwind *ww;
  CameraRecd *CPos;
  int        spdimen, cpdimen;
  double     *pc;
  boolean    found, drawsidemenu;

  ww = er->data0;
  id = er->id & 0x03;
  CPos = &ww->CPos[id];
  switch ( msg ) {
case xgemsg_3DWIN_RESIZE:
case xgemsg_3DWIN_PROJCHANGE:
    drawsidemenu = false;
    if ( editing_camera ) {
      SetupCameraParamWidgets ();
      drawsidemenu = true;
    }
    if ( rendered_picture ) {
      if ( ww->fww.zoomwin != -1 && ww->fww.zoomwin != 3 ) {
        renderbtn0->data0 = renderbtn1->data0 = renderbtn2->data0 = txtRender;
        rendered_picture = false;
        drawsidemenu = true;
      }
      else
        Geom00WinRestartRendering ();
    }
    xge_SetClipping ( er );
    er->redraw ( er, true );
    if ( drawsidemenu ) {
      xge_SetClipping ( side00menu );
      side00menu->redraw ( side00menu, true );
    }
    return 1;

case xgemsg_3DWIN_PICK_POINT:
    if ( editing_shapefunc )
      return EdShapeFuncFindNearestPoint ( CPos, x, y );
    else if ( editing_lights )
      return EdLightFindNearestPoint ( CPos, x, y );
    else {
      found = GeomObjectFindCPoint ( 3, CPos, x, y );
      if ( found ) {
        if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                       &spdimen, &cpdimen, &pc ) )
          BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
      }
      return found;
    }

case xgemsg_3DWIN_MOVE_POINT:
    if ( editing_shapefunc )
      EdShapeFuncSetPoint ( CPos, x, y );
    else if ( editing_lights )
      EdLightSetPoint ( CPos, x, y );
    else {
      GeomObjectSetCPoint ( CPos, x, y );
      if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                     &spdimen, &cpdimen, &pc ) )
        BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
    }
    xge_SetClipping ( ww->fww.er );
    rendered_picture = false;
    ww->fww.er->redraw ( ww->fww.er, true );
    return 1;

case xgemsg_3DWIN_SELECT_POINTS:
    if ( editing_shapefunc || editing_lights )
      ;
    else if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        if ( GeomObjectSelectPoint ( 3, CPos, x, y ) ) {
          bsm_vertex_num0 = ((GO_BSplineMesh*)current_go)->current_vertex[0];
          bsm_vertex_num1 = ((GO_BSplineMesh*)current_go)->current_vertex[1];
          xge_RedrawAll ();
        }
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        if ( GeomObjectSelectEdge ( 3, CPos, x, y ) ) {
          bsm_edge_num0 = ((GO_BSplineMesh*)current_go)->current_edge[0];
          bsm_edge_num1 = ((GO_BSplineMesh*)current_go)->current_edge[1];
          xge_RedrawAll ();
        }
        return false;
      }
    }
    GeomObjectMarkCPoints ( 3, CPos, &ww->selection_rect, MARK_CP_SELECT );
    xge_SetClipping ( ww->fww.er );
    ww->fww.er->redraw ( ww->fww.er, true );
    return 1;

case xgemsg_3DWIN_UNSELECT_POINTS:
    if ( editing_shapefunc || editing_lights )
      ;
    else
      GeomObjectMarkCPoints ( 3, CPos, &ww->selection_rect, MARK_CP_UNSELECT );
    xge_SetClipping ( ww->fww.er );
    ww->fww.er->redraw ( ww->fww.er, true );
    return 1;

case xgemsg_3DWIN_TGSELECT_POINTS:
    if ( editing_shapefunc || editing_lights )
      ;
    else
      GeomObjectMarkCPoints ( 3, CPos, &ww->selection_rect, MARK_CP_TGSELECT );
    xge_SetClipping ( ww->fww.er );
    ww->fww.er->redraw ( ww->fww.er, true );
    return 1;

case xgemsg_3DWIN_SPECIAL_SELECT:
    if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        if ( GeomObjectSelectPoint ( 3, CPos, x, y ) ) {
          bsm_vertex_num0 = ((GO_BSplineMesh*)current_go)->current_vertex[0];
          bsm_vertex_num1 = ((GO_BSplineMesh*)current_go)->current_vertex[1];
          if ( GeomObjectGetPointCoord ( current_go, bsm_vertex_num0,
                                         &spdimen, &cpdimen, &pc ) )
            BottomDisplayPoint ( win1, spdimen, cpdimen,
                                 bsm_vertex_num0, pc, false );
          xge_RedrawAll ();
        }
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        if ( GeomObjectSelectEdge ( 3, CPos, x, y ) ) {
          bsm_edge_num0 = ((GO_BSplineMesh*)current_go)->current_edge[0];
          bsm_edge_num1 = ((GO_BSplineMesh*)current_go)->current_edge[1];
          xge_RedrawAll ();
        }
        return false;
      }
    }
    return false;

case xgemsg_3DWIN_SPECIAL_UNSELECT:
    if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        GeomObjectUnselectPoint ( 3 );
        bsm_vertex_num0 = bsm_vertex_num1 = -1;
        xge_RedrawAll ();
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        GeomObjectUnselectEdge ( 3 );
        bsm_edge_num0 = bsm_edge_num1 = -1;
        xge_RedrawAll ();
        return false;
      }
    }
    return false;

case xgemsg_3DWIN_CHANGE_TRANS:
    ww = er->data1;
    Geom00Win3DShowTransOrigin ( ww );
    return 1;

case xgemsg_3DWIN_SAVE_POINTS:
    if ( editing_pretrans )
      GeomObjectSavePretransformation ( current_go );
    else if ( editing_shapefunc )
      EdShapeFuncSavePoints ();
    else if ( editing_lights )
      EdLightSavePoints ();
    else
      GeomObjectSaveCPoints ( 3 );
    return 1;

case xgemsg_3DWIN_TRANSFORM_POINTS:
    Geom00Win3DShowTransformation ( ww );
    if ( editing_pretrans ) {
      GeomObjectCompPretransformation ( current_go, &ww->gwtrans );
      GetPreTransformation ();
      xge_SetClipping ( popup14 );
      popup14->redraw ( popup14, true );
    }
    else if ( editing_shapefunc )
      EdShapeFuncTransformPoints ( &ww->gwtrans );
    else if ( editing_lights )
      EdLightTransformPoints ( &ww->gwtrans );
    else {
      GeomObjectTransformCPoints3D ( &ww->gwtrans, marking_mask );
      rendered_picture = false;
    }
    return 1;

case xgemsg_3DWIN_FIND_REFBBOX:
    if ( GeomObjectFindBoundingBox ( 3, &ww->RefBBox ) )
      return 1;
    else
      return 0;

case xgemsg_3DWIN_UNDO:
    GeomObjectUndoLastTransformation ( 3 );
    rendered_picture = false;
    return 1;

case xgemsg_3DWIN_KEY:
    switch ( key ) {
  case 'P': case 'p':
      Geom00Win3DShowCameraPos ( &ww->CPos[3] );
      return 1;
  default:
      return 0;
    }

case xgemsg_3DWIN_ERROR:
    return 0;

default:
    return 0;
  }
} /*Geom00Win3DCallBack*/

void DrawG00win2D ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;

  _2Dwin = er->data0;
  glXWaitX ();
  xgle_DrawGeomWinBackground ( er, 0 );
  if ( _2Dwin->display_coord && _2Dwin->inside )
    xgle_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );
  xgle_SetGLCamerad ( &_2Dwin->CPos );
  GeomObjectDisplayActive ( 2 );
  glXWaitGL ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );
  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*DrawG00win2D*/

void Geom00Win2DShowTransformation ( xge_2Dwind *ww )
{
  char text[256], txt1[32];
  vector2d *sf;

  text[0] = 0;
  switch ( ww->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    sf = &ww->trans_params;
    sprintf ( text, "moving: %7.3f, %7.3f", sf->x, sf->y );
    break;
case xge_2DWIN_SCALING_TOOL:
    sf = &ww->scaling_factors;
    sprintf ( text, "scaling: %7.3f, %7.3f", sf->x, sf->y );
    break;
case xge_2DWIN_ROTATING_TOOL:
    pkv_RadToDegreeStr ( ww->trans_params.x, txt1 );
    sprintf ( text, "rotating: %s (%8.5f)", txt1, ww->trans_params.x );
    break;
case xge_2DWIN_SHEAR_TOOL:
    switch ( ww->tool_mode ) {
  case 1:
  case 2:
      sprintf ( text, "shearing: %7.3f, %7.3f",
                ww->shear_axis[0].x, ww->shear_axis[0].y );
      break;
  case 3:
  case 4:
      sprintf ( text, "shearing: %7.3f, %7.3f",
                ww->shear_axis[1].x, ww->shear_axis[1].y );
      break;
  default:
      break;
    }
    break;
default:
    break;
  }
  SetStatusText ( text, true );
} /*Geom00Win2DShowTransformation*/

void Geom00Win2DShowTransOrigin ( xge_2Dwind *ww )
{
  char text[256];

  text[0] = 0;
  switch ( ww->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    break;
case xge_2DWIN_SCALING_TOOL:
    sprintf ( text, "x: %7.3f, y: %7.3f",
              ww->scaling_centre.x, ww->scaling_centre.y );
    break;
case xge_2DWIN_ROTATING_TOOL:
    sprintf ( text, "x: %7.3f, y: %7.3f",
              ww->rotating_centre.x, ww->rotating_centre.y );
    break;
case xge_2DWIN_SHEAR_TOOL:
    switch ( ww->tool_mode ) {
  case 0:
      sprintf ( text, "x: %7.3f, y: %7.3f",
                ww->shear_centre.x, ww->shear_centre.y );
      break;
  case 1:
  case 2:
      sprintf ( text, "x: %7.3f, y: %7.3f",
                ww->shear_axis[0].x, ww->shear_axis[0].y );
      break;
  case 3:
  case 4:
      sprintf ( text, "x: %7.3f, y: %7.3f",
                ww->shear_axis[1].x, ww->shear_axis[1].y );
      break;
  default:
      break;
    }
    break;
default:
    break;
  }
  SetStatusText ( text, true );
} /*Geom00Win2DShowTransOrigin*/

int Geom00Win2DCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_2Dwind *ww;
  CameraRecd *CPos;
  Box3d      box;
  boolean    found;
  int        spdimen, cpdimen;
  double     *pc;

  ww = er->data0;
  CPos = &ww->CPos;
  switch ( msg ) {
case xgemsg_2DWIN_RESIZE:
case xgemsg_2DWIN_PROJCHANGE:
    return 0;

case xgemsg_2DWIN_PICK_POINT:
    found = GeomObjectFindCPoint ( 2, CPos, x, y );
    if ( found ) {
      if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                     &spdimen, &cpdimen, &pc ) )
        BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
    }
    return found;

case xgemsg_2DWIN_MOVE_POINT:
    GeomObjectSetCPoint ( CPos, x, y );
    if ( GeomObjectGetPointCoord ( currentp_go, current_point_ind,
                                   &spdimen, &cpdimen, &pc ) )
      BottomDisplayPoint ( win0, spdimen, cpdimen, current_point_ind, pc, true );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_SELECT_POINTS:
    if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        if ( GeomObjectSelectPoint ( 2, CPos, x, y ) ) {
          bsm_vertex_num0 = ((GO_BSplineMesh*)current_go)->current_vertex[0];
          bsm_vertex_num1 = ((GO_BSplineMesh*)current_go)->current_vertex[1];
          xge_RedrawAll ();
        }
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        if ( GeomObjectSelectEdge ( 2, CPos, x, y ) ) {
          bsm_edge_num0 = ((GO_BSplineMesh*)current_go)->current_edge[0];
          bsm_edge_num1 = ((GO_BSplineMesh*)current_go)->current_edge[1];
          xge_RedrawAll ();
        }
        return false;
      }
    }
    GeomObjectMarkCPoints ( 2, CPos, &ww->selection_rect, MARK_CP_SELECT );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_UNSELECT_POINTS:
    GeomObjectMarkCPoints ( 2, CPos, &ww->selection_rect, MARK_CP_UNSELECT );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_TGSELECT_POINTS:
    GeomObjectMarkCPoints ( 2, CPos, &ww->selection_rect, MARK_CP_TGSELECT );
    xge_SetClipping ( er );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_SPECIAL_SELECT:
    if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        if ( GeomObjectSelectPoint ( 2, CPos, x, y ) ) {
          bsm_vertex_num0 = ((GO_BSplineMesh*)current_go)->current_vertex[0];
          bsm_vertex_num1 = ((GO_BSplineMesh*)current_go)->current_vertex[1];
          if ( GeomObjectGetPointCoord ( current_go, bsm_vertex_num0,
                                         &spdimen, &cpdimen, &pc ) )
            BottomDisplayPoint ( win1, spdimen, cpdimen,
                                 bsm_vertex_num0, pc, false );
          xge_RedrawAll ();
        }
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        if ( GeomObjectSelectEdge ( 2, CPos, x, y ) ) {
          bsm_edge_num0 = ((GO_BSplineMesh*)current_go)->current_edge[0];
          bsm_edge_num1 = ((GO_BSplineMesh*)current_go)->current_edge[1];
          xge_RedrawAll ();
        }
        return false;
      }
    }
    return 1;

case xgemsg_2DWIN_SPECIAL_UNSELECT:
    if ( current_go->obj_type == GO_BSPLINE_MESH ) {
      if ( sw_bsm_selectvertex ) {
        GeomObjectUnselectPoint ( 2 );
        bsm_vertex_num0 = bsm_vertex_num1 = -1;
        xge_RedrawAll ();
        return false;
      }
      else if ( sw_bsm_selectedge ) {
        GeomObjectUnselectEdge ( 2 );
        bsm_edge_num0 = bsm_edge_num1 = -1;
        xge_RedrawAll ();
        return false;
      }
    }
    return 1;

case xgemsg_2DWIN_CHANGE_TRANS:
    Geom00Win2DShowTransOrigin ( ww );
    return 1;

case xgemsg_2DWIN_SAVE_POINTS:
    GeomObjectSaveCPoints ( 2 );
    return 1;

case xgemsg_2DWIN_TRANSFORM_POINTS:
    Geom00Win2DShowTransformation ( ww );
    GeomObjectTransformCPoints2D ( &ww->gwtrans, marking_mask );
    return 1;

case xgemsg_2DWIN_FIND_REFBBOX:
    if ( GeomObjectFindBoundingBox ( 2, &box ) ) {
      ww->RefBBox.x0 = box.x0;  ww->RefBBox.x1 = box.x1;
      ww->RefBBox.y0 = box.y0;  ww->RefBBox.y1 = box.y1;
      return 1;
    }
    else
      return 0;

case xgemsg_2DWIN_UNDO:
    GeomObjectUndoLastTransformation ( 3 );
    er->redraw ( er, true );
    return 1;

case xgemsg_2DWIN_ERROR:
    return 0;

default:
    return 0;
  }
} /*Geom00Win2DCallBack*/

void SetWin002D ( void )
{
  if ( geom00win != geom00win2D ) {
    geom00win = geom00win2D;
    xge_SetMenuWidgets ( geom00menu, geom00win2D, false );
    g00win2D.display_coord = swwin0coordinates;
    SetToolSwitches ( g00win2D.current_tool, g00win2D.display_coord );
    geom00win2D->msgproc ( geom00win2D, xgemsg_RESIZE, 0,
                           geom00menu->w, geom00menu->h );
  }
} /*SetWin002D*/

void SetWin003D ( void )
{
  if ( geom00win != geom00win3D ) {
    rendered_picture = false;
    geom00win = geom00win3D;
    g00win3D.display_coord = swwin0coordinates;
    xge_SetMenuWidgets ( geom00menu, geom00win3D, false );
    SetToolSwitches ( g00win3D.current_tool, g00win3D.display_coord );
    geom00win3D->msgproc ( geom00win3D, xgemsg_RESIZE, 0,
                           geom00menu->w, geom00menu->h );
  }
} /*SetWin003D*/

int Geom00MenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
default:
    return 0;
  }
} /*Geom00MenuCallBack*/

