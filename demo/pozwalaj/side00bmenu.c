
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
#include <sys/times.h>

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
#include "render.h"
#include "edlight.h"


static clock_t tick;
static long    ticksps;

void InitSide00bMenu ( void )
{
  xge_widget *w;

        /* rendering control widgets */
  w = renderbtn0 = xge_NewButton ( win0, NULL, btnM01bRENDER_INTERRUPT,
                                   58, 19, 0, 21, txtRender );
  w = xge_NewSwitch ( win0, w, swM01bPICTURE_ON, 16, 16, 60, 22, NULL,
                      &rendered_picture );
  w = xge_NewButton ( win0, w, btnM01bCAMERA, 58, 19, 0, 41, txtCamera );
  w = xge_NewTextWidget ( win0, w, 0, 12, 16, 0, 63, txtC);
  w = xge_NewTextWidget ( win0, w, 0, 89, 16, 20, 63, txtD_shape_func );
  w = xge_NewSwitch ( win0, w, swM01bGAUSSIAN_C, 16, 16, 0, 81, NULL,
                      &swGaussian_c );
  w = xge_NewSwitch ( win0, w, swM01bMEAN_C, 16, 16, 0, 101, NULL,
                      &swMean_c );
  w = xge_NewSwitch ( win0, w, swM01bLAMBERTISO_C, 16, 16, 0, 121, NULL,
                      &swLambert_c );
  w = xge_NewSwitch ( win0, w, swM01bREFLECTION_C, 16, 16, 0, 141, NULL,
                      &swReflection_c );
  w = xge_NewSwitch ( win0, w, swM01bHIGHLIGHT_C, 16, 16, 0, 161, NULL,
                      &swHighlight_c );
  w = xge_NewSwitch ( win0, w, swM01bSECTIONS_C, 16, 16, 0, 181, NULL,
                      &swSections_c );
  w = xge_NewSwitch ( win0, w, swM01bGAUSSIAN_D, 89, 16, 20, 81, txtGaussian,
                      &swGaussian_d );
  w = xge_NewSwitch ( win0, w, swM01bMEAN_D, 89, 16, 20, 101, txtMean,
                      &swMean_d );
  w = xge_NewSwitch ( win0, w, swM01bLAMBERTISO_D, 89, 16, 20, 121, txtIsophotes,
                      &swLambert_d );
  w = xge_NewSwitch ( win0, w, swM01bREFLECTION_D, 89, 16, 20, 141, txtReflection,
                      &swReflection_d );
  w = xge_NewSwitch ( win0, w, swM01bHIGHLIGHT_D, 89, 16, 20, 161, txtHighlight,
                      &swHighlight_d );
  w = xge_NewSwitch ( win0, w, swM01bSECTIONS_D, 89, 16, 20, 181, txtSections,
                      &swSections_d );
  w = xge_NewSlidebar2d ( win0, w, slM01bREND_CFRANGE, 109, 10, 0, 201,
                          render_cfrange );
  w = xge_NewSlidebard ( win0, w, slM01bREND_DFSF, 109, 10, 0, 215,
                         &render_dfsf );
  w = xge_NewSwitch ( win0, w, swM10bREFLECTION_FRAME, 109, 16, 0, 229,
                      txtReflectionFrame, &sw_reflection_frame );
  w = xge_NewSwitch ( win0, w, swM01bHIGHLIGHT_FRAME, 109, 16, 0, 249,
                      txtHighlightFrame, &sw_highlight_frame );
  w = xge_NewSwitch ( win0, w, swM01bSECTIONS_FRAME, 109, 16, 0, 269,
                      txtSectionsNormal, &sw_sections_frame );
        /* status & command line on/off */
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00bwidgets = w;

        /* light editing widgets */
  w = renderbtn1 = xge_NewButton ( win0, NULL, btnM01bRENDER_INTERRUPT,
                                   58, 19, 0, 21, txtRender );
  w = xge_NewSwitch ( win0, w, swM01bPICTURE_ON, 16, 16, 60, 22, NULL,
                      &rendered_picture );
  w = xge_NewButton ( win0, w, btnM01bSHAPE_FUNCTION,
                      58, 19, 0, 41, txtShapeFunction );
  w = xge_NewSwitch ( win0, w, swM01bLIGHT0DIR, 89, 16, 0, 62, txtLight0,
                      &sw_edit_light[0] );
  w = xge_NewSlidebard ( win0, w, slM01bLIGHT0INT, 109, 10, 0, 80,
                         &render_light_int[0] );
  w = xge_NewSwitch ( win0, w, swM01bLIGHT1DIR, 89, 16, 0, 92, txtLight1,
                      &sw_edit_light[1] );
  w = xge_NewSlidebard ( win0, w, slM01bLIGHT1INT, 109, 10, 0, 110,
                         &render_light_int[1] );
  w = xge_NewSwitch ( win0, w, swM01bLIGHT2DIR, 89, 16, 0, 122, txtLight2,
                      &sw_edit_light[2] );
  w = xge_NewSlidebard ( win0, w, slM01bLIGHT2INT, 109, 10, 0, 140,
                         &render_light_int[2] );
  w = xge_NewSwitch ( win0, w, swM01bLIGHT3DIR, 89, 16, 0, 152, txtLight3,
                      &sw_edit_light[3] );
  w = xge_NewSlidebard ( win0, w, slM01bLIGHT3INT, 109, 10, 0, 170,
                         &render_light_int[3] );
  w = xge_NewTextWidget ( win0, w, 0, 109, 16, 0, 182, txtAmbient );
  w = xge_NewSlidebard ( win0, w, 0, 109, 10, 0, 200, &render_light_int[4] );
  w = xge_NewSwitch ( win0, w, swM01bSHADOWS, 89, 16, 0, 216, txtShadows,
                      &swShadows );
  w = xge_NewSwitch ( win0, w, swM01bANTIALIAS, 89, 16, 0, 236, txtAntialias,
                      &swAntialias );
  w = xge_NewIntWidget ( win0, w, intwM01bRENDERINGTHREADS, 109, 19, 0, 256,
                         1, MAX_PTHREADS, &side00bnpthreads,
                         txtNPThreads, &rendering_npthreads );
        /* status & command line on/off */
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00cwidgets = w;

        /* camera parameters */
  w = renderbtn2 = xge_NewButton ( win0, NULL, btnM01bRENDER_INTERRUPT,
                                   58, 19, 0, 21, txtRender );
  w = xge_NewSwitch ( win0, w, swM01bPICTURE_ON, 16, 16, 60, 22, NULL,
                      &rendered_picture );
  w = xge_NewButton ( win0, w, btnM01bLIGHT,
                      58, 19, 0, 41, txtLight );
  w = xge_NewTextWidget ( win0, w, 0, 100, 15, 3, 62, txtPosition );
  w = xge_NewTextWidget ( win0, w, 0, 12, 16, 12, 83, txtX );
  w = xge_NewTextWidget ( win0, w, 0, 12, 16, 12, 103, txtY );
  w = xge_NewTextWidget ( win0, w, 0, 12, 16, 12, 123, txtZ );
  w = xge_NewStringEd ( win0, w, intwM01b_POS_X, 81, 19, 24, 82,
                        MAX_PARAM_LGT, side00fparam_str[0], &side00fparam_ed[0] );
  w = xge_NewStringEd ( win0, w, intwM01b_POS_Y, 81, 19, 24, 102,
                        MAX_PARAM_LGT, side00fparam_str[1], &side00fparam_ed[1] );
  w = xge_NewStringEd ( win0, w, intwM01b_POS_Z, 81, 19, 24, 122,
                        MAX_PARAM_LGT, side00fparam_str[2], &side00fparam_ed[2] );
  w = xge_NewDiald ( win0, w, dialM01b_PSI, 32, 32, 2, 148, NULL, &side00fpsi );
  w = xge_NewDiald ( win0, w, dialM01b_THETA, 32, 32, 2, 192, NULL, &side00ftheta );
  w = xge_NewDiald ( win0, w, dialM01b_PHI, 32, 32, 2, 238, NULL, &side00fphi );
  w = xge_NewStringEd ( win0, w, intwM01b_PSI, 69, 19, 35, 155,
                        MAX_PARAM_LGT, side00fparam_str[3], &side00fparam_ed[3] );
  w = xge_NewStringEd ( win0, w, intwM01b_THETA, 69, 19, 35, 199,
                        MAX_PARAM_LGT, side00fparam_str[4], &side00fparam_ed[4] );
  w = xge_NewStringEd ( win0, w, intwM01b_PHI, 69, 19, 35, 245,
                        MAX_PARAM_LGT, side00fparam_str[5], &side00fparam_ed[5] );
  w = xge_NewTextWidget ( win0, w, 0, 12, 16, 12, 279, txtF );
  w = xge_NewStringEd ( win0, w, intwM01b_F, 81, 19, 24, 280,
                        MAX_PARAM_LGT, side00fparam_str[6], &side00fparam_ed[6] );
        /* status & command line on/off */
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00fwidgets = w;
} /*InitSide00bMenu*/

boolean InitRendering ( void )
{
  xge_SetWindow ( win0 );
  RendReset ();
  RendEnterCamerad ( &g00win3D.CPos[3], g00win3D.fww.win[3] );
  RendEnterLightsd ( RENDERING_NLIGHTS, render_light_dir, render_light_int );
  RendEnterReflectionLinesFramed ( edshapef_reflection_frame );
  RendEnterHighlightLinesFramed ( edshapef_highlight_frame );
  RendEnterSectionPlanesNormald ( &edshapef_sectiondir );
  GeomObjectOutputToRenderer3D ( false );
  if ( RendBegin () ) {
    if ( !rendered_picture ) {
      XGetSubImage ( xgedisplay, xgepixmap, g00win3D.cwin[3]->x, g00win3D.cwin[3]->y,
                     g00win3D.cwin[3]->w, g00win3D.cwin[3]->h, 0xFFFFFFFF, ZPixmap,
                     rendimage, g00win3D.cwin[3]->x, g00win3D.cwin[3]->y );
      rendered_picture = true;
    }
    xge_ReleaseFocus ( xge_null_widget );
    xge_Redraw ();
    return true;
  }
  else {
    xge_ReleaseFocus ( xge_null_widget );
    rendered_picture = false;
    return false;
  }
} /*InitRendering*/

void StartRendering ( void )
{
  struct tms tt;

  if ( geom00win != geom00win3D ) {
    xge_DisplayErrorMessage ( ErrorMsgCannotRender, -1 );
    return;
  }
  if ( g00win3D.fww.zoomwin != -1 && g00win3D.fww.zoomwin != 3 ) {
    g00win3D.fww.zoomwin = -1;
    xge_CompSizeFourWW ( g00win3D.fww.er, 1 );
    xge_3DwindUpdatePerspProj ( &g00win3D );
    rendered_picture = false;
    g00win3D.fww.er->redraw ( g00win3D.fww.er, true );
    rendered_picture = true;
  }
  renderbtn0->data0 = renderbtn1->data0 = renderbtn2->data0 = txtInterrupt;
  ticksps = sysconf ( _SC_CLK_TCK );
  tick = times ( &tt );
  xge_GrabFocus ( xge_null_widget, true );
  xge_SetWindowCursor ( win0, xgeCURSOR_WATCH );
  xge_SetWindowCursor ( win1, xgeCURSOR_WATCH );
  xge_PostIdleCommand ( IDLE_COMMAND_RENDER_START, 0, 0 );
  xge_SetWindow ( win0 );
  xge_Redraw ();  /* there may be a popup */
} /*StartRendering*/

void StopRendering ( void )
{
  RenderingIsOn = false;  /* in case of interrupt */
  renderbtn0->data0 = renderbtn1->data0 = renderbtn2->data0 = txtRender;
  xge_SetWindow ( win0 );
  xge_Redraw ();  /* there may be a popup */
} /*StopRendering*/

void ContinueRendering ( void )
{
  clock_t    tick1;
  struct tms tt;
  xge_widget *er;

  if ( !rendered_picture ) {
    StopRendering ();
  }
  else {
    RenderLine ();
    tick1 = times ( &tt );
    if ( (tick1-tick) >= ticksps || !RenderingIsOn ) {
      er = g00win3D.fww.win[3];
      xge_SetClipping ( er );
      er->redraw ( er, false );
      xge_RedrawPopups ();
      xge_SetClipping ( er );
      xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
      tick = tick1;
    }
    if ( RenderingIsOn )
      xge_PostIdleCommand ( IDLE_COMMAND_RENDER_CONT, 0, 0 );
    else
      StopRendering ();
  }
} /*ContinueRendering*/

void SetupCameraParamWidgets ( void )
{
  CameraRecd *CPos;

  CPos = &g00win3D.CPos[3];
  sprintf ( side00fparam_str[0], "%f", CPos->position.x );
  sprintf ( side00fparam_str[1], "%f", CPos->position.y );
  sprintf ( side00fparam_str[2], "%f", CPos->position.z );
  sprintf ( side00fparam_str[3], "%f", CPos->psi );
  sprintf ( side00fparam_str[4], "%f", CPos->theta );
  sprintf ( side00fparam_str[5], "%f", CPos->phi );
  sprintf ( side00fparam_str[6], "%f", CPos->vd.persp.f );
  side00fpsi   = CPos->psi;
  side00ftheta = CPos->theta;
  side00fphi   = CPos->phi;
  editing_camera = true;
} /*SetupCameraParamWidgets*/

void UpdateCameraPosition ( CameraRecd *CPos, char which, double ang )
{
  double   r;
  vector3d v;

  SubtractPoints3d ( &CPos->g_centre, &CPos->position, &v );
  r = sqrt ( DotProduct3d ( &v, &v ) );
  CPos->position = CPos->g_centre;
  switch ( which ) {
case 0:
    CPos->psi = ang;
    break;
case 1:
    CPos->theta = ang;
    break;
case 2:
    CPos->phi = ang;
    break;
  }
  SetVector3d ( &v, 0.0, 0.0, -r );
  CameraMoveCd ( CPos, &v );
  SetupCameraParamWidgets ();
} /*UpdateCameraPosition*/

int Side00bMenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  CameraRecd *CPos;
  double     a;

  CPos = &g00win3D.CPos[3];
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM01bRENDER_INTERRUPT:
      if ( RenderingIsOn ) {
        StopRendering ();
        rendered_picture = false;
      }
      else
        StartRendering ();
      return 1;
  case btnM01bLIGHT:
      xge_SetMenuWidgets ( side00menu, side00cwidgets, editing_lights );
      if ( !editing_lights ) {
        picture_lights = editing_lights = true;
        editing_shapefunc = editing_camera = false;
        xge_Redraw ();
      }
      return 1;
  case btnM01bSHAPE_FUNCTION:
      xge_SetMenuWidgets ( side00menu, side00bwidgets, editing_shapefunc );
      if ( !editing_shapefunc ) {
        picture_lights = editing_lights = editing_camera = false;
        editing_shapefunc = true;
        xge_Redraw ();
      }
      return 1;
  case btnM01bCAMERA:
      SetupCameraParamWidgets ();
      xge_SetMenuWidgets ( side00menu, side00fwidgets, editing_camera );
      if ( !editing_camera ) {
        editing_camera = true;
        picture_lights = editing_lights = editing_shapefunc = false;
        xge_Redraw ();
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swM01bGAUSSIAN_C:
      if ( swGaussian_c ) {
        swMean_c = swLambert_c = swReflection_c =
        swHighlight_c = swSections_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bMEAN_C:
      if ( swMean_c ) {
        swGaussian_c = swLambert_c = swReflection_c =
        swHighlight_c = swSections_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bLAMBERTISO_C:
      if ( swLambert_c ) {
        swGaussian_c = swMean_c = swReflection_c =
        swHighlight_c = swSections_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bREFLECTION_C:
      if ( swReflection_c ) {
        swGaussian_c = swMean_c = swLambert_c =
        swHighlight_c = swSections_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bHIGHLIGHT_C:
      if ( swHighlight_c ) {
        swGaussian_c = swMean_c = swLambert_c =
        swReflection_c = swSections_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bSECTIONS_C:
      if ( swSections_c ) {
        swGaussian_c = swMean_c = swLambert_c =
        swReflection_c = swHighlight_c = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bGAUSSIAN_D:
      if ( swGaussian_d ) {
        swMean_d = swLambert_d = swReflection_d =
        swHighlight_d = swSections_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bMEAN_D:
      if ( swMean_d ) {
        swGaussian_d = swLambert_d = swReflection_d =
        swHighlight_d = swSections_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bLAMBERTISO_D:
      if ( swLambert_d ) {
        swGaussian_d = swMean_d = swReflection_d =
        swHighlight_d = swSections_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bREFLECTION_D:
      if ( swReflection_d ) {
        swGaussian_d = swMean_d = swLambert_d =
        swHighlight_d = swSections_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bHIGHLIGHT_D:
      if ( swHighlight_d ) {
        swGaussian_d = swMean_d = swLambert_d =
        swReflection_d = swSections_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM01bSECTIONS_D:
      if ( swSections_d ) {
        swGaussian_d = swMean_d = swLambert_d =
        swReflection_d = swHighlight_d = false;
        xge_SetClipping ( side00menu );
        side00menu->redraw ( side00menu, true );
      }
      return 1;
  case swM10bREFLECTION_FRAME:
  case swM01bHIGHLIGHT_FRAME:
  case swM01bSECTIONS_FRAME:
      geom00menu->redraw ( geom00menu, true );
      return 1;
  case swM01bSHADOWS:
      return 1;
  case swM01bANTIALIAS:
      return 1;
  case swM01bLIGHT0DIR:
  case swM01bLIGHT1DIR:
  case swM01bLIGHT2DIR:
  case swM01bLIGHT3DIR:
      xge_Redraw ();
      return 1;
  case swM01bPICTURE_ON:
      if ( rendered_picture ) {
        Geom00WinRestartRendering ();
        xge_Redraw ();
      }
      else if ( RenderingIsOn )
        StopRendering ();
      else
        xge_Redraw ();
      return 1;
  default:
      return 0;
    }

case xgemsg_INT_WIDGET_COMMAND:
    switch ( er->id ) {
  case intwM01bRENDERINGTHREADS:
      if ( key >= 1 && key <= MAX_PTHREADS )
        rendering_npthreads = ncpu > 1 ? key : 1;
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    switch ( er->id ) {
  case slM01bREND_DFSF:
      return 0;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR2_COMMAND:
    switch ( er->id ) {
  case slM01bREND_CFRANGE:
      return 0;
  default:
      return 0;
    }

case xgemsg_DIAL_COMMAND:
    switch ( er->id ) {
  case dialM01b_PSI:
      UpdateCameraPosition ( CPos, 0, side00fpsi );
      goto redraw_it;
  case dialM01b_THETA:
      UpdateCameraPosition ( CPos, 1, side00ftheta );
      goto redraw_it;
  case dialM01b_PHI:
      UpdateCameraPosition ( CPos, 2, side00fphi );
      goto redraw_it;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case intwM01b_POS_X:
      if ( sscanf ( side00fparam_str[0], "%lf", &a ) == 1 ) {
        CPos->position.x = a;
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_POS_Y:
      if ( sscanf ( side00fparam_str[1], "%lf", &a ) == 1 ) {
        CPos->position.y = a;
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_POS_Z:
      if ( sscanf ( side00fparam_str[2], "%lf", &a ) == 1 ) {
        CPos->position.z = a;
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_PSI:
      if ( sscanf ( side00fparam_str[3], "%lf", &a ) == 1 ) {
        if ( a < -PI ) a = -PI;
        else if ( a > PI ) a = PI;
        sprintf ( side00fparam_str[3], "%f", a );
        side00fpsi = a;
        UpdateCameraPosition ( CPos, 0, side00fpsi );
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_THETA:
      if ( sscanf ( side00fparam_str[4], "%lf", &a ) == 1 ) {
        if ( a < -PI ) a = -PI;
        else if ( a > PI ) a = PI;
        sprintf ( side00fparam_str[4], "%f", a );
        side00ftheta = a;
        UpdateCameraPosition ( CPos, 1, side00ftheta );
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_PHI:
      if ( sscanf ( side00fparam_str[5], "%lf", &a ) == 1 ) {
        if ( a < -PI ) a = -PI;
        else if ( a > PI ) a = PI;
        sprintf ( side00fparam_str[5], "%f", a );
        side00fphi = a;
        UpdateCameraPosition ( CPos, 2, side00fphi );
        goto redraw_it;
      }
      else
        return 0;
  case intwM01b_F:
      if ( sscanf ( side00fparam_str[6], "%lf", &a ) == 1 ) {
        if ( a < xge_3DWIN_MIN_ZOOM )      a = xge_3DWIN_MIN_ZOOM;
        else if ( a > xge_3DWIN_MAX_ZOOM ) a = xge_3DWIN_MAX_ZOOM;
        sprintf ( side00fparam_str[6], "%f", a );
        CPos->vd.persp.f = a;
      }
      else
        return 0;
redraw_it:
      xge_SetClipping ( side00menu );
      side00menu->redraw ( side00menu, true );
      xge_SetClipping ( g00win3D.cwin[3] );
      if ( rendered_picture )
        Geom00WinRestartRendering ();
      else
        g00win3D.cwin[3]->redraw ( g00win3D.cwin[3], true );
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case intwM01b_POS_X:
  case intwM01b_POS_Y:
  case intwM01b_POS_Z:
  case intwM01b_PSI:
  case intwM01b_THETA:
  case intwM01b_PHI:
  case intwM01b_F:
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Side00bMenuCallBack*/

