
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "xgledit.h"

#include "linkage.h"
#include "palm.h"
#include "ludzik.h"
#include "anima.h"
#include "poruszaj.h"
#include "lighting.h"
#include "widgets.h"


/* ///////////////////////////////////////////////////////////////////////// */
xge_3Dwind    wind;
xge_KnotWind  knwind;

xge_widget    *topmenu0, *sidemenu0, *geommenu0, *bottommenu0,
              *popup0, *popup1, *popup2, *popup3;
xge_widget    *side0a, *side0b, *side0c;
xge_widget    *bottom0a, *bottom0b;

xge_listbox   dirlist1, filelist1, dirlist2, filelist2;
xge_string_ed filename_editor;
const char    file_filter[] = "*.key";
const char    file_ext[] = ".key";
char          initial_directory[MAX_PATH_LGT+1];
char          current_directory[MAX_PATH_LGT+1], current_dir[MAX_PATH_SHRT+1];
char          filename[MAX_FILENAME_LGT+1];

xge_scroll_widget scrollw0;

boolean sw_palmlod1 = false, sw_palmlod2 = false, sw_palmlod3 = true;
boolean sw_antialias = false;

/* ///////////////////////////////////////////////////////////////////////// */
static boolean animate = false;
static clock_t time_cnt = 0;
static double  time_start, time_now = 0.0;

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *TopMenu0Init ( xge_widget *prev )
{
  xge_widget *w, *m;

  w = xge_NewButton ( 0, NULL, btnT0FILE, 64, 18, 0, 0, txtFile );
  w = xge_NewButton ( 0, w, btnT0ARTICULATE, 64, 18, 66, 0, txtArticulate );
  w = xge_NewButton ( 0, w, btnT0OPTIONS, 64, 18, 132, 0, txtOptions );
  w = xge_NewButton ( 0, w, btnT0LIGHT, 64, 18, 198, 0, txtLight );
  w = xge_NewButton ( 0, w, btnT0ABOUT, 64, 18, xge_WIDTH-64, 0, txtAbout );
  xge_SetWidgetPositioning ( w, 1, -60, 0 );
  m = xge_NewMenu ( 0, prev, 0, xge_WIDTH, 20, 0, 0, w );
  return m;
} /*TopMenu0Init*/

int TopMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnT0FILE:
      xge_AddPopup ( popup0 );
      xge_GrabFocus ( popup0, true );
      return 1;
  case btnT0ARTICULATE:
      xge_SetMenuWidgets ( sidemenu0, side0a, true );
      return 1;
  case btnT0OPTIONS:
      xge_SetMenuWidgets ( sidemenu0, side0b, true );
      return 1;
  case btnT0LIGHT:
      xge_SetMenuWidgets ( sidemenu0, side0c, true );
      return 1;
  case btnT0ABOUT:
      xge_DisplayInfoMessage ( InfoMsg, -1 );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*TopMenu0CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *SideMenu0Init ( xge_widget *prev )
{
  xge_widget *w, *scm, *m;
  int        i;

        /* side menu a - articulation parameters */
  w = NULL;
  for ( i = 0; i < N_PARAMS; i++ ) {
    w = xge_NewSlidebard ( 0, w, slS0aARTPARAM0+i, 109, 10, 0, 2+15*i,
                           &art_param[i] );
    xge_SetWidgetPositioning ( w, 0, 2, 2+15*i );
  }
  scm = xge_NewMenu ( 0, NULL, 0, 111, 15*N_PARAMS, 0, 0, w );
  w = xge_NewScrollWidget ( 0, NULL, 0, SIDE_MENU_WIDTH, xge_HEIGHT-62, 0, 20,
                            &scrollw0, scm );

  w = xge_NewSwitch ( 0, w, swS0aPANZOOM, 109, 16, 0, xge_HEIGHT-16,
                       txtPanZoom, &wind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( 0, w, swS0aCOORDS, 109, 16, 0, xge_HEIGHT-36,
                       txtCoordinates, &wind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  side0a = w;

        /* side menu b - options */
  w = xge_NewSwitch ( 0, NULL, swS0bPALMLOD1, 16, 16, 0, 20, NULL, &sw_palmlod1 );
  w = xge_NewSwitch ( 0, w, swS0bPALMLOD2, 16, 16, 20, 20, NULL, &sw_palmlod2 );
  w = xge_NewSwitch ( 0, w, swS0bPALMLOD3, 70, 16, 40, 20, txtPalmLOD, &sw_palmlod3 );
  w = xge_NewSwitch ( 0, w, swS0bSPLHANDS, 90, 16, 0, 40,
                      txtSplineHands, &sw_spline_hands );
  w = xge_NewSwitch ( 0, w, swS0bSPLNET, 90, 16, 0, 60,
                      txtDisplayNets, &sw_draw_bsnets );
  w = xge_NewSwitch ( 0, w, swS0bANTIALIAS, 90, 16, 0, 90,
                      txtAntialias, &sw_antialias );
  w = xge_NewSwitch ( 0, w, swS0aPANZOOM, 109, 16, 0, xge_HEIGHT-16,
                       txtPanZoom, &wind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( 0, w, swS0aCOORDS, 109, 16, 0, xge_HEIGHT-36,
                       txtCoordinates, &wind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  side0b = w;

        /* side menu c - light parameters */
  w = xge_NewDialf ( 0, NULL, diS0cANG1, 32, 32, 0, 20,
                     NULL, &HorizAngle );
  w = xge_NewDialf ( 0, w, diS0cANG2, 32, 32, 36, 20,
                     NULL, &VertAngle );
  w = xge_NewSwitch ( 0, w, swS0cUSESPEC, 80, 16, 0, 56, txtSpecular, &useSpecular );

  w = xge_NewSwitch ( 0, w, swS0aPANZOOM, 109, 16, 0, xge_HEIGHT-16,
                       txtPanZoom, &wind.panning );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( 0, w, swS0aCOORDS, 109, 16, 0, xge_HEIGHT-36,
                       txtCoordinates, &wind.display_coord );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  side0c = w;


  m = xge_NewMenu ( 0, prev, SIDEMENU0a,
                    SIDE_MENU_WIDTH, xge_HEIGHT-TOP_MENU_HEIGHT,
                    0, TOP_MENU_HEIGHT, side0a );
  return m;
} /*SideMenu0Init*/

int SideMenu0aCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swS0aPANZOOM:
      if ( wind.panning )
        xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_PANNING_TOOL );
      else
        xge_3DwindEnableGeomWidget ( &wind, xge_3DWIN_NO_TOOL );
      knwind.panning = wind.panning;
      xge_Redraw ();
      return 1;
  case swS0aCOORDS:
      knwind.display_coord = wind.display_coord;
      return 1;
  default:
      return 0;
    }

case xgemsg_SLIDEBAR_COMMAND:
    if ( er->id >= slS0aARTPARAM0 && er->id < slS0aARTPARAM0+N_PARAMS ) {
      xge_SetClipping ( geommenu0 );
      geommenu0->redraw ( geommenu0, true );
      return 1;
    }
    else
      return 0;

default:
    return 0;
  }
} /*SideMenu0aCallBack*/

int SideMenu0bCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swS0bPALMLOD1:
      if ( sw_palmlod1 )
        { sw_palmlod2 = sw_palmlod3 = false;  palm_lod = 1; }
      else
        { sw_palmlod2 = true; sw_palmlod3 = false;  palm_lod = 2; }
      goto redraw;
  case swS0bPALMLOD2:
      if ( sw_palmlod2 )
        { sw_palmlod1 = sw_palmlod3 = false;  palm_lod = 2; }
      else
        { sw_palmlod3 = true; sw_palmlod1 = false;  palm_lod = 3; }
      goto redraw;
  case swS0bPALMLOD3:
      if ( sw_palmlod3 )
        { sw_palmlod1 = sw_palmlod2 = false;  palm_lod = 3; }
      else
        { sw_palmlod1 = true; sw_palmlod2 = false;  palm_lod = 1; }
      goto redraw;
  case swS0bSPLHANDS:
      goto redraw;
  case swS0bSPLNET:
      goto redraw;
  case swS0bANTIALIAS:
      if ( _xgle_no_accum )
        sw_antialias = false;
redraw:
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*SideMenu0bCallBack*/

int SideMenu0cCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swS0cUSESPEC:
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

case xgemsg_DIAL_COMMAND:
    switch ( er->id ) {
  case diS0cANG1:
  case diS0cANG2:
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*SideMenu0cCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawWinContents ( void )
{
  glPushMatrix ();
  glTranslated ( 0.0, 0.28, 0.0 );
  glScaled ( 0.2, 0.2, 0.2 );
  glEnable ( GL_DEPTH_TEST );
  RenderLightPosition ();
  SetupLightToDisplay ();
  DisplayCharacter ();
  glDisable ( GL_DEPTH_TEST );
  glPopMatrix ();
} /*DrawWinContents*/

void RysujOkno ( xge_widget *er, boolean onscreen )
{
#define TTOL 0.05
#define NTR  9
  int     id, i;
  clock_t ct;
  static double pshx[NTR] = {-0.4, -0.4, -0.4, 0.0, 0.0, 0.0, 0.4, 0.4, 0.4 };
  static double pshy[NTR] = {-0.4, 0.0, 0.4, -0.4, 0.0, 0.4, -0.4, 0.0, 0.4 };

  if ( animate ) {
    if ( AnimaCurrentTime ( &ct ) - time_now > TTOL )
      Animate ( false );
  }
  glXWaitX ();
  id = er->id & 0x03;

  if ( sw_antialias && !_xgle_no_accum ) {
    glClear ( GL_ACCUM_BUFFER_BIT );
    for ( i = 0; i < NTR; i++ ) {
      xgle_DrawGeomWinBackground ( er, GL_DEPTH_BUFFER_BIT );
      xgle_SetIdentMapping ( er );
      xgle_3DwindDrawCursorPos ( &wind, id, xge_xx, xge_yy );
      xgle_SetGLaccCamerad ( &wind.CPos[id], pshx[i], pshy[i],
                             0.0, 0.0, wind.CPos[id].zmin );
      DrawWinContents ();
      if ( i == 0 )
        glAccum ( GL_LOAD, 1.0/(double)NTR );
      else
        glAccum ( GL_ACCUM, 1.0/(double)NTR );
    }
    glAccum ( GL_RETURN, 1.0 );
  }
  else {
    xgle_DrawGeomWinBackground ( er, GL_DEPTH_BUFFER_BIT );
    xgle_SetIdentMapping ( er );
    xgle_3DwindDrawCursorPos ( &wind, id, xge_xx, xge_yy );
    xgle_SetGLCamerad ( &wind.CPos[id] );
    DrawWinContents ();
  }

  glXWaitGL ();
  glFlush ();
  xgleCopyGLRect ( er->w, er->h, er->x, er->y );

  xge_DrawGeomWinSelectionRect ( er, &wind.selection_rect );
  xge_3DwinfDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
#undef TTOL
} /*RysujOkno*/

xge_widget *GeomMenu0Init ( xge_widget *prev )
{
  xge_widget *m;

  m = xge_New3Dwind ( 0, prev, GEOMMENU0, xge_WIDTH-SIDE_MENU_WIDTH,
                      xge_HEIGHT-TOP_MENU_HEIGHT-BOTTOM_MENU_HEIGHT,
                      SIDE_MENU_WIDTH, TOP_MENU_HEIGHT, &wind, RysujOkno, RysujOkno );
  CameraSetDepthRanged ( &wind.CPos[0], -10.0, 10.0 );
  CameraSetDepthRanged ( &wind.CPos[1], -10.0, 10.0 );
  CameraSetDepthRanged ( &wind.CPos[2], -10.0, 10.0 );
  CameraSetDepthRanged ( &wind.CPos[3],  10.0, 40.0 );
  return m;
} /*GeomMenu0Init*/

int GeomMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
default:
    return 0;
  }
} /*GeomMenu0CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void MyDrawKnotWind ( xge_widget *er, boolean onscreen )
{
  xge_KnotWind *knw;
  short      x, y0, y1;

  knw = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( knw->display_coord && knw->xx >= 0 )
    xge_KnotWindDrawCursorPos ( knw );
  xge_KnotWindDrawAxis ( knw );

  y0 = (short)(er->y+2);
  y1 = (short)(er->y+er->h-2);
  if ( animate )
    xgeSetForeground ( xgec_Yellow );
  else
    xgeSetForeground ( xgec_Green3 );
  x = xge_KnotWindMapKnot ( knw, time_now );
  xgeDrawLine ( x, y0, x, y1 );

  xge_KnotWindDrawKnots ( knw );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*MyDrawKnotWind*/

void MyDrawButton ( xge_widget *er, boolean onscreen )
{
  XPoint p[3];
  short  xc, yc;

  xge_DrawButton ( er, false );
  xgeSetForeground ( xgec_White );
  xc = er->x + er->w/2;
  yc = er->y + er->h/2;
  switch ( er->id ) {
case btnB0PREVKEY:
    xgeDrawLine ( xc-5, yc, xc+5, yc );
    xgeDrawLine ( xc-5, yc, xc-1, yc-4 );
    xgeDrawLine ( xc-5, yc, xc-1, yc+4 );
    break;
case btnB0NEXTKEY:
    xgeDrawLine ( xc-5, yc, xc+5, yc );
    xgeDrawLine ( xc+5, yc, xc+1, yc-4 );
    xgeDrawLine ( xc+5, yc, xc+1, yc+4 );
    break;
case btnB0PLAYSTOP:
    if ( animate ) {
      xgeFillRectangle ( 2, 10, xc-4, yc-5 );
      xgeFillRectangle ( 2, 10, xc+2, yc-5 );
    }
    else {
      p[0].x = p[2].x = xc-5;  p[1].x = xc+5;
      p[0].y = yc-5;  p[1].y = yc;  p[2].y = yc+5;
      xgeFillPolygon ( Convex, 3, p );
    }
    break;
default:
    break;
  }
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*MyDrawButton*/

xge_widget *BottomMenu0Init ( xge_widget *prev )
{
  xge_widget *w, *m;

  w = xge_NewButton ( 0, NULL, btnB0PLAY, 58, 19,
                      SIDE_MENU_WIDTH, xge_HEIGHT-19, txtPlay );
  xge_SetWidgetPositioning ( w, 2, 0, -19 );
  w = xge_NewButton ( 0, w, btnB0PREVKEY, 32, 19,
                      SIDE_MENU_WIDTH+60, xge_HEIGHT-19, NULL );
  w->redraw = MyDrawButton;
  xge_SetWidgetPositioning ( w, 2, 60, -19 );
  w = xge_NewButton ( 0, w, btnB0SETPOSE, 32, 19,
                      SIDE_MENU_WIDTH+94, xge_HEIGHT-19, txtSet );
  w->redraw = MyDrawButton;
  xge_SetWidgetPositioning ( w, 2, 94, -19 );
  w = xge_NewButton ( 0, w, btnB0NEXTKEY, 32, 19,
                      SIDE_MENU_WIDTH+128, xge_HEIGHT-19, NULL );
  w->redraw = MyDrawButton;
  xge_SetWidgetPositioning ( w, 2, 128, -19 );
  bottom0a = w;

  w = xge_NewButton ( 0, NULL, btnB0EDIT, 58, 19,
                      SIDE_MENU_WIDTH, xge_HEIGHT-19, txtEdit );
  xge_SetWidgetPositioning ( w, 2, 0, -19 );
  w = xge_NewButton ( 0, w, btnB0PLAYSTOP, 32, 19,
                      SIDE_MENU_WIDTH+60, xge_HEIGHT-19, NULL );
  w->redraw = MyDrawButton;
  xge_SetWidgetPositioning ( w, 2, 60, -19 );
  w = xge_NewSwitch ( 0, w, swB0PERIODIC, 66, 16,
                      SIDE_MENU_WIDTH+94, xge_HEIGHT-18, txtPeriodic, &sw_periodic );
  xge_SetWidgetPositioning ( w, 2, 94, -18 );
  bottom0b = w;

  knwind.er = w =
      xge_NewKnotWind ( 0, bottom0a, knwB0KNOTWIN, xge_WIDTH-SIDE_MENU_WIDTH, 20,
                        SIDE_MENU_WIDTH, xge_HEIGHT-BOTTOM_MENU_HEIGHT+1,
                        &knwind, MAX_KEY_POSES, ic_knots );
  w->redraw = MyDrawKnotWind;
  knwind.degree = 1;
  knwind.lastknot = GetNPoses () + 1;
  m = xge_NewMenu ( 0, prev, BOTTOMMENU0,
                    xge_WIDTH-SIDE_MENU_WIDTH, BOTTOM_MENU_HEIGHT,
                    SIDE_MENU_WIDTH, xge_HEIGHT-BOTTOM_MENU_HEIGHT, w );
  return m;
} /*BottomMenu0Init*/

void AnimaOn ( void )
{
  struct tms tt;

  time_cnt = times ( &tt );
  time_start = time_now;
  xge_PostIdleCommand ( cmdANIMATE, 0, 0 );
} /*AnimaOn*/

double AnimaCurrentTime ( clock_t *ct )
{
#define TIME_INTERVAL 1.0
  struct tms tt;
  clock_t    cnt;

  *ct = cnt = times ( &tt );
  return time_start + ((double)(cnt-time_cnt))/
                       (TIME_INTERVAL*sysconf(_SC_CLK_TCK));
#undef TIME_INTERVAL
} /*AnimaCurrentTime*/

boolean Animate ( boolean nudge )
{
  int     nposes;
  clock_t ct;

  if ( animate ) {
/*printf ( "*" );*/
    time_now = AnimaCurrentTime ( &ct );
    nposes = GetNPoses ();
    while ( time_now > ic_knots[nposes] ) {
      time_start = time_now += ic_knots[1]-ic_knots[nposes];
      time_cnt = ct;
    }
    if ( GetPose ( time_now, art_param ) ) {
      if ( nudge )
        xge_PostIdleCommand ( cmdANIMATE, 0, 0 );
      return true;
    }
    else
      animate = false;
  }
  return false;
} /*Animate*/

int BottomMenu0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int i, nkposes;

  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnB0PLAY:
      knwind.locked = true;
      knwind.er->prev = bottom0b;
      bottom0b->next = knwind.er;
      xge_SetMenuWidgets ( bottommenu0, knwind.er, true );
      return 1;
  case btnB0PREVKEY:
      nkposes = GetNPoses ();
      for ( i = nkposes; i > 1; i-- )
        if ( time_now > ic_knots[i] )
          break;
      time_now = ic_knots[i];
      GetPose ( time_now, art_param );
      xge_Redraw ();
      return 1;
  case btnB0SETPOSE:
      if ( EnterKeyPose ( time_now, art_param ) ) {
        knwind.lastknot = GetNPoses () + 1;
        xge_SetClipping ( bottommenu0 );
        bottommenu0->redraw ( bottommenu0, true );
      }
      return 1;
  case btnB0NEXTKEY:
      nkposes = GetNPoses ();
      for ( i = 1; i <= nkposes; i++ )
        if ( time_now < ic_knots[i] )
          break;
      time_now = ic_knots[i];
      GetPose ( time_now, art_param );
      xge_Redraw ();
      return 1;
  case btnB0EDIT:
      knwind.locked = animate = false;
      knwind.er->prev = bottom0a;
      bottom0a->next = knwind.er;
      xge_SetMenuWidgets ( bottommenu0, knwind.er, true );
      return 1;
  case btnB0PLAYSTOP:
      animate = !animate;
      xge_SetClipping ( er );
      er->redraw ( er, true );
      if ( animate )
        AnimaOn ();
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swB0PERIODIC:
      spline_ok = false;
      if ( animate ) {
        GetPose ( time_now, art_param );
        animate = false;
      }
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

case xgemsg_KNOTWIN_CHANGE_KNOT:
    if ( er->id == knwB0KNOTWIN ) {
/* printf ( "change\n" ); */
      if ( MoveKeyKnot ( knwind.current_knot-1 ) ) {
        time_now = ic_knots[knwind.current_knot];
        xge_SetClipping ( geommenu0 );
        geommenu0->redraw ( geommenu0, true );
        return 1;
      }
      else
        return 0;
    }
    else
      return 0;

case xgemsg_KNOTWIN_INSERT_KNOT:
    if ( er->id == knwB0KNOTWIN ) {
/* printf ( "insert\n" ); */
      if ( EnterKeyPose ( knwind.newknot, art_param ) ) {
        time_now = knwind.newknot;
        knwind.lastknot = GetNPoses () + 1;
        xge_SetClipping ( geommenu0 );
        geommenu0->redraw ( geommenu0, true );
        return 1;
      }
      else
        return 0;
    }
    else
      return 0;

case xgemsg_KNOTWIN_REMOVE_KNOT:
    if ( er->id == knwB0KNOTWIN ) {
/* printf ( "remove\n" ); */
      if ( DeleteKeyPose ( knwind.current_knot-1 ) ) {
        knwind.lastknot = GetNPoses () + 1;
        xge_SetClipping ( geommenu0 );
        geommenu0->redraw ( geommenu0, true );
        return 1;
      }
      else
        return 0;
    }
    else
      return 0;

case xgemsg_KNOTWIN_CHANGE_MAPPING:
    if ( er->id == knwB0KNOTWIN ) {
      return 1;
    }
    else
      return 0;

case xgemsg_KNOTWIN_MCLICK:
    if ( er->id == knwB0KNOTWIN ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
        if ( key & xgemouse_LBUTTON_CHANGE )
          xge_GrabFocus ( er, false );
        goto point_time;
      }
      else if ( er == xge_GetFocusWidget ( 0 ) )
        xge_ReleaseFocus ( er );
      return 1;
    }
    else
      return 0;

case xgemsg_KNOTWIN_MMOVE:
    if ( er->id == knwB0KNOTWIN ) {
      if ( key & xgemouse_LBUTTON_DOWN ) {
point_time:
        animate = false;
        time_now = xge_KnotWindUnmapKnot ( &knwind, x );
        GetPose ( time_now, art_param );
        xge_Redraw ();
      }
      else if ( er == xge_GetFocusWidget ( 0 ) )
        xge_ReleaseFocus ( er );
      return 1;
    }
    else
      return 0;

case xgemsg_KNOTWIN_ERROR:
    if ( er->id == knwB0KNOTWIN ) {
      return 1;
    }
    else
      return 0;

default:
    return 0;
  }
} /*BottomMenu0CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *Popup0Init ( void )
{
  xge_widget *w, *m;

  w = xge_NewButton ( 0, NULL, btnP0OPEN, 62, 19, 2, 22, txtOpen );
  w = xge_NewButton ( 0, w, btnP0SAVE, 62, 19, 2, 42, txtSave );
  w = xge_NewButton ( 0, w, btnP0SAVEAS, 62, 19, 2, 62, txtSaveAs );
  w = xge_NewButton ( 0, w, btnP0EXIT, 62, 19, 2, 82, txtExit );
  m = xge_NewFMenu ( 0, NULL, POPUP0, 66, 83, 0, 20, w );
  if ( m )
    m->msgproc = xge_PopupMenuMsg;
  return m;
} /*Popup0Init*/

void SetCurrentDir ( void )
{
  int l;

  getcwd ( current_directory, MAX_PATH_LGT+1 );
  l = strlen ( current_directory );
  if ( l <= MAX_PATH_SHRT )
    strcpy ( current_dir, current_directory );
  else {
    strcpy ( current_dir, "..." );
    strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
  }
} /*SetCurrentDir*/

boolean FilenameCorrect ( char *fn )
{
        /* !!!TU POWINNO BYC GRUNTOWNE SPRAWDZANIE POPRAWNOSCI NAZWY!!! */
  return fn[0] != 0;
} /*FilenameCorrect*/

void OpenSaveAsPopup ( void )
{
  SetCurrentDir ();
  xge_SetupFileList ( &filelist2, ".", file_filter, false );
  xge_SetupDirList ( &dirlist2, ".", NULL, false, NULL );
  xge_AddPopup ( popup2 );
  xge_GrabFocus ( popup2, true );
} /*OpenSaveAsPopup*/

int Popup0CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP0OPEN:
      SetCurrentDir ();
      xge_SetupFileList ( &filelist1, ".", file_filter, false );
      xge_SetupDirList ( &dirlist1, ".", NULL, false, NULL );
      xge_RemovePopup ( true );
      xge_AddPopup ( popup1 );
      xge_GrabFocus ( popup1, true );
      return 1;
  case btnP0SAVE:
      if ( FilenameCorrect ( filename ) ) {
           /* tu wstawiamy wywolanie procedury pisania danych w pliku */

           /* okno popup mozemy zostawic, jesli blad */
        xge_RemovePopup ( true );
        return 1;
      }
      else {
        xge_RemovePopup ( true );
        OpenSaveAsPopup ();
      }
      return 1;
  case btnP0SAVEAS:
      xge_RemovePopup ( true );
      OpenSaveAsPopup ();
      return 1;
  case btnP0EXIT:
      xge_RemovePopup ( true );
      xge_AddPopup ( popup3 );
      xge_GrabFocus ( popup3, true );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup0CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *Popup1Init ( void )
{
  xge_widget *w, *m;

  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir );
  w = xge_NewButton ( 0, w, btnP1OPEN, 58, 19, 20+91, 40+180-30, txtOpen );
  w = xge_NewButton ( 0, w, btnP1CANCEL, 58, 19, 20+251, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lbP1DIRLIST, 180, 99, 20+10, 40+40, &dirlist1 );
  w = xge_NewListBox ( 0, w, lbP1FILELIST, 180, 99, 220+10, 40+40, &filelist1 );
  m = xge_NewFMenu ( 0, NULL, POPUP1, 400, 180, 20, 40, w );
  if ( m )
    m->msgproc = xge_PopupMenuMsg;
  return m;
} /*Popup1Init*/

void Popup1ChangeDir ( void )
{
  if ( !chdir ( &dirlist1.itemstr[dirlist1.itemind[dirlist1.currentitem]] ) ) {
    xge_SetupFileList ( &filelist1, ".", file_filter, false );
    xge_SetupDirList ( &dirlist1, ".", NULL, false, current_directory );
    SetCurrentDir ();
    xge_SetClipping ( popup1 );
    popup1->redraw ( popup1, true );
  }
} /*Popup1ChangeDir*/

boolean Popup1OpenFile ( void )
{
    /* procedura czytania danych - tymczasem rozwiazanie najprostsze */
  FILE   *f;
  int    i, j, nkp, d;
  double kpostime[MAX_KEY_POSES], pose[N_PARAMS];

  f = fopen ( filename, "r+" );
  if ( !f )
    return false;
  i = fscanf ( f, "%d %d", &d, &nkp );
  if ( i != 2 || d != N_PARAMS || nkp < 0 || nkp > MAX_KEY_POSES )
    goto failure;
        /* empty the key pose table */
  SetNArtParams ( N_PARAMS );

  for ( i = 0; i < nkp; i++ )
    if ( fscanf ( f, "%lf", &kpostime[i] ) != 1 )
      goto failure;
  for ( i = 0; i < nkp; i++ ) {
    for ( j = 0; j < N_PARAMS; j++ )
      if ( fscanf ( f, "%lf", &pose[j] ) != 1 )
        goto failure;
    if ( !EnterKeyPose ( kpostime[i], pose ) )
      goto failure;
  }
  fclose ( f );

  knwind.lastknot = GetNPoses () + 1;
  xge_KnotWindFindMapping ( &knwind );
  time_now = kpostime[0];
  GetPose ( time_now, art_param );
  return true;

failure:
  fclose ( f );
  return false;
} /*Popup1OpenFile*/

int Popup1CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP1OPEN:
      if ( xge_GetCurrentListBoxString ( &filelist1, filename ) ) {
        xge_RemovePopup ( false );
        Popup1OpenFile ();
        xge_Redraw ();
      }
      return 1;
  case btnP1CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lbP1DIRLIST:
      Popup1ChangeDir ();
      return 1;
  case lbP1FILELIST:
      if ( xge_GetCurrentListBoxString ( &filelist1, filename ) ) {
        xge_RemovePopup ( false );
        Popup1OpenFile ();
        xge_Redraw ();
      }
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup1CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *Popup2Init ( void )
{
  xge_widget *w, *m;

  w = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir );
  w = xge_NewButton ( 0, w, btnP2SAVE, 58, 19, 20+91, 40+180-30, txtSave );
  w = xge_NewButton ( 0, w, btnP2CANCEL, 58, 19, 20+251, 40+180-30, txtCancel );
  w = xge_NewListBox ( 0, w, lbP2DIRLIST, 180, 99, 20+10, 40+40, &dirlist2 );
  w = xge_NewListBox ( 0, w, lbP2FILELIST, 180, 67, 220+10, 40+72, &filelist2 );
  w = xge_NewTextWidget ( 0, w, 0, 120, 16, 220+10, 40+35, txtSaveAs );
  w = xge_NewStringEd ( 0, w, txtedP2FILENAME, 180, 19, 220+10, 40+50,
                        MAX_FILENAME_LGT, filename, &filename_editor );
  m = xge_NewFMenu ( 0, NULL, POPUP2, 400, 180, 20, 40, w );
  if ( m )
    m->msgproc = xge_PopupMenuMsg;
  return m;
} /*Popup2Init*/

void Popup2ChangeDir ( void )
{
  if ( !chdir ( &dirlist2.itemstr[dirlist2.itemind[dirlist2.currentitem]] ) ) {
    xge_SetupFileList ( &filelist2, ".", file_filter, false );
    xge_SetupDirList ( &dirlist2, ".", NULL, false, current_directory );
    SetCurrentDir ();
    xge_SetClipping ( popup2 );
    popup1->redraw ( popup2, true );
  }
} /*Popup2ChangeDir*/

boolean Popup2SaveFile ( void )
{
        /* tu procedura zapisywania danych na w pliku */
  FILE   *f;
  int    i, j, np;
  double pose[N_PARAMS];

  f = fopen ( filename, "w+" );
  if ( !f )
    return false;
  np = GetNPoses ();
  fprintf ( f, "%d %d\n", N_PARAMS, np );
  for ( i = 0; i < np; i++ )
    fprintf ( f, "%7.4f ", ic_knots[i+1] );
  fprintf ( f, "\n" );
  for ( i = 0; i < np; i++ ) {
    GetPose ( ic_knots[i+1], pose );
    for ( j = 0; j < N_PARAMS; j++ )
      fprintf ( f, "%7.4f", pose[j] );
    fprintf ( f, "\n" );
  }
  fclose ( f );
  return true;
} /*Popup2SaveFile*/

int Popup2CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP2SAVE:
      if ( FilenameCorrect ( filename ) ) {
        xge_RemovePopup ( true );
        Popup2SaveFile ();
      }
      else
        xge_DisplayErrorMessage ( ErrorMessageIncorrectFilename, -1 );
      return 1;
  case btnP2CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lbP2DIRLIST:
      Popup2ChangeDir ();
      return 1;
  case lbP2FILELIST:
      if ( xge_GetCurrentListBoxString ( &filelist2, filename ) ) {
        filename_editor.start = filename_editor.pos = 0;
        xge_SetClipping ( popup2 );
        popup2->redraw ( popup2, true );
      }
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    return 1;

case xgemsg_TEXT_EDIT_ENTER:
    return 1;

default:
    return 0;
  }
} /*Popup2CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *Popup3Init ( void )
{
  xge_widget *w, *m;

  w = xge_NewTextWidget ( 0, NULL, 0, 240, 16, 110, 70, txtMsgReconsider );
  w = xge_NewButton ( 0, w, btnP3EXIT, 58, 19, 110, 100, txtExit );
  w = xge_NewButton ( 0, w, btnP3SAVE, 58, 19, 210, 100, txtSave );
  w = xge_NewButton ( 0, w, btnP3CANCEL, 58, 19, 310, 100, txtCancel );
  m = xge_NewFMenu ( 0, NULL, POPUP3, 300, 70, 90, 60, w );
  if ( m )
    m->msgproc = xge_PopupMenuMsg;
  return m;
} /*Popup3Init*/

int Popup3CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP3EXIT:
      xge_done = 1;
      return 1;
  case btnP3SAVE:
      xge_RemovePopup ( true );
      OpenSaveAsPopup ();
      return 1;
  case btnP3CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup3CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
int NonWidgetCallBack ( int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_KEY:
    switch ( key ) {
  case 'Q': case 'q':
      xge_done = 1;
      return 1;
  case 'm':
      xgeResizeWindow ( xge_WIDTH, xge_HEIGHT );
      return 1;
  case 'M':
      xgeResizeWindow ( xge_MAX_WIDTH, xge_MAX_HEIGHT );
      return 1;
  default:
      return 0;
    }

case xgemsg_IDLE_COMMAND:
    switch ( key ) {
  case cmdANIMATE:
      Animate ( true );
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

case xgemsg_RESIZE:
    wind.fww.er->msgproc ( wind.fww.er, xgemsg_RESIZE, 0,
                 (short)(xge_current_width-SIDE_MENU_WIDTH),
                 (short)(xge_current_height-TOP_MENU_HEIGHT-BOTTOM_MENU_HEIGHT) );
    sidemenu0->msgproc ( sidemenu0, xgemsg_RESIZE, 0, SIDE_MENU_WIDTH,
                         (short)(xge_current_height-TOP_MENU_HEIGHT) );
    scrollw0.er->msgproc ( scrollw0.er, xgemsg_RESIZE, 0, SIDE_MENU_WIDTH,
                           (short)(xge_current_height-62) );
    topmenu0->msgproc ( topmenu0, xgemsg_RESIZE, 0, xge_current_width, 20 );
    bottommenu0->y = xge_current_height-BOTTOM_MENU_HEIGHT;
    knwind.er->y = bottommenu0->y+1;
    knwind.er->msgproc ( knwind.er, xgemsg_RESIZE, 0,
                         (short)(xge_current_width-SIDE_MENU_WIDTH), 20 );
    bottommenu0->msgproc ( bottommenu0, xgemsg_RESIZE, 0,
                           (short)(xge_current_width-SIDE_MENU_WIDTH),
                           BOTTOM_MENU_HEIGHT );
    xge_Redraw ();
    return 1;

default:
    return 0;
  }
} /*NonWidgetCallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
/*    printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
    switch ( er->id & MENUMASK ) {
  case TOPMENU0:
      return TopMenu0CallBack ( er, msg, key, x, y );
  case SIDEMENU0a:
      return SideMenu0aCallBack ( er, msg, key, x, y );
  case SIDEMENU0b:
      return SideMenu0bCallBack ( er, msg, key, x, y );
  case SIDEMENU0c:
      return SideMenu0cCallBack ( er, msg, key, x, y );
  case GEOMMENU0:
      return GeomMenu0CallBack ( er, msg, key, x, y );
  case BOTTOMMENU0:
      return BottomMenu0CallBack ( er, msg, key, x, y );
  case POPUP0:
      return Popup0CallBack ( er, msg, key, x, y );
  case POPUP1:
      return Popup1CallBack ( er, msg, key, x, y );
  case POPUP2:
      return Popup2CallBack ( er, msg, key, x, y );
  case POPUP3:
      return Popup3CallBack ( er, msg, key, x, y );
  default:
      return 0;
    }
  }
  else {
/*    printf ( "%d, %d, %d, %d\n", msg, key, x, y ); */
    return NonWidgetCallBack ( msg, key, x, y );
  }
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void init_edwin ( void )
{
  pkv_InitScratchMem ( 655360 );
  getcwd ( initial_directory, MAX_PATH_LGT+1 );

  InitCharacter ();
  if ( !SetNArtParams ( N_PARAMS ) ) {
    printf ( "Error: could not initialise animation tables." );
    exit ( 1 );
  }
  EnterKeyPose ( 0.0, art_param );
  EnterKeyPose ( 1.0, art_param );
  InitSpotlight ();

    /* kolejnosc tworzenia tych wihajstrow jest istotna */
  topmenu0 = TopMenu0Init ( NULL );
  sidemenu0 = SideMenu0Init ( topmenu0 );
  geommenu0 = GeomMenu0Init ( sidemenu0 );
  bottommenu0 = BottomMenu0Init ( geommenu0 );
  xge_SetWinEdRect ( bottommenu0 );
    /* tych juz nie */
  popup0 = Popup0Init ();
  popup1 = Popup1Init ();
  popup2 = Popup2Init ();
  popup3 = Popup3Init ();

  xge_Redraw ();
  xge_DisplayInfoMessage ( InfoMsg, -1 );
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

/* ///////////////////////////////////////////////////////////////////////// */
int main ( int argc, char *argv[] )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  xgle_Init ( argc, argv, CallBack, NULL,
              XGLE_WANTED, XGLE_WOULD_BE_NICE, XGLE_NOT_NEEDED );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xgle_Cleanup ();
  exit ( 0 );
} /*main*/

