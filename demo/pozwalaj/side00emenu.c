
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
#include "edlight.h"
#include "arrowbtn.h"


void InitSide00eMenu ( void )
{
  xge_widget *w, *ww;

  w = xge_NewTextWidget ( win0, NULL, 0, 100, 19, 0, 0, txtRefPoint );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 20, txtX );
  w = xge_NewStringEd ( win0, w, textedM01ePX, 93, 19, 15, 20,
                        MAX_PARAM_LGT, side00eparam_str[0], &side00eparam_ed[0] );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 40, txtY );
  w = xge_NewStringEd ( win0, w, textedM01ePY, 93, 19, 15, 40,
                        MAX_PARAM_LGT, side00eparam_str[1], &side00eparam_ed[1] );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 60, txtZ );
  w = xge_NewStringEd ( win0, w, textedM01ePZ, 93, 19, 15, 60,
                        MAX_PARAM_LGT, side00eparam_str[2], &side00eparam_ed[2] );
  w = xge_NewTextWidget ( win0, w, 0, 100, 19, 1, 80, txtRefVector );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 100, txtX );
  w = xge_NewStringEd ( win0, w, textedM01eVX, 93, 19, 15, 100,
                        MAX_PARAM_LGT, side00eparam_str[3], &side00eparam_ed[3] );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 120, txtY );
  w = xge_NewStringEd ( win0, w, textedM01eVY, 93, 19, 15, 120,
                        MAX_PARAM_LGT, side00eparam_str[4], &side00eparam_ed[4] );
  w = xge_NewTextWidget ( win0, w, 0, 12, 19, 1, 140, txtZ );
  w = xge_NewStringEd ( win0, w, textedM01eVZ, 93, 19, 15, 140,
                        MAX_PARAM_LGT, side00eparam_str[5], &side00eparam_ed[5] );
  w = xge_NewTextWidget ( win0, w, 0, 100, 19, 1, 160, txtRadius );
  w = xge_NewStringEd ( win0, w, textedM01eR, 93, 19, 15, 180,
                        MAX_PARAM_LGT, side00eparam_str[6], &side00eparam_ed[6] );
  w = xge_NewTextWidget ( win0, w, 0, 100, 19, 1, 200, txtAngleDeg );
  w = xge_NewStringEd ( win0, w, textedM01eANGLE, 93, 19, 15, 220,
                        MAX_PARAM_LGT, side00eparam_str[7], &side00eparam_ed[7] );
  w = xge_NewButton ( win0, w, btnM01eTRANSLATE, 72, 19, 0, 244, txtTranslate );
  w = xge_NewButton ( win0, w, btnM01eSCALE, 72, 19, 0, 264, txtScale );
  w = xge_NewButton ( win0, w, btnM01eROTATE, 72, 19, 0, 284, txtRotate );
  w = xge_NewTextWidget ( win0, w, 0, 100, 19, 1, 304, txtProject );
  w = xge_NewButton ( win0, w, btnM01ePROJ_LINE, 72, 19, 0, 324, txtOnLine );
  w = xge_NewButton ( win0, w, btnM01ePROJ_PLANE, 72, 19, 0, 344, txtOnPlane );
  w = xge_NewButton ( win0, w, btnM01ePROJ_PLANE_UP, 18, 19, 73, 344, NULL );
  w->redraw = DrawButtonArrowUp;
  w = xge_NewButton ( win0, w, btnM01ePROJ_PLANE_DOWN, 18, 19, 92, 344, NULL );
  w->redraw = DrawButtonArrowDown;
  w = xge_NewButton ( win0, w, btnM01ePROJ_SPHERE, 72, 19, 0, 364, txtOnSphere );
  w = xge_NewButton ( win0, w, btnM01ePROJ_SPHERE_UP, 18, 19, 73, 364, NULL );
  w->redraw = DrawButtonArrowUp;
  w = xge_NewButton ( win0, w, btnM01ePROJ_SPHERE_DOWN, 18, 19, 92, 364, NULL );
  w->redraw = DrawButtonArrowDown;
  w = xge_NewButton ( win0, w, btnM01ePROJ_CYLINDER, 72, 19, 0, 384, txtOnCylinder );
  w = xge_NewButton ( win0, w, btnM01ePROJ_CYLINDER_UP, 18, 19, 73, 384, NULL );
  w->redraw = DrawButtonArrowUp;
  w = xge_NewButton ( win0, w, btnM01ePROJ_CYLINDER_DOWN, 18, 19,92, 384, NULL );
  w->redraw = DrawButtonArrowDown;
  for ( ww = w; ww; ww = ww->prev )
    xge_SetWidgetPositioning ( ww, 0, ww->x, ww->y );
  side00econtents = xge_NewMenu ( win0, NULL, scwM01CONTENTS,
                                  SIDEMENUWIDTH0, 407, 0, 20, w );
  side00escroll = xge_NewScrollWidget ( win0, NULL, scwM01SCROLL,
                                  SIDEMENUWIDTH0, xge_HEIGHT-TOPMENUHEIGHT-40, 0, 20,
                                  &side00esw, side00econtents );

        /* status & command line on/off */
  w = xge_NewSwitch ( win0, side00escroll, swM01COORDINATES, 100, 16, 0, xge_HEIGHT-36,
                      txtCoordinates, &swwin0coordinates );
  xge_SetWidgetPositioning ( w, 2, 0, -36 );
  w = xge_NewSwitch ( win0, w, swM01STATUS, 16, 16, 0, xge_HEIGHT-16,
                      txtNull, &win0statusline );
  xge_SetWidgetPositioning ( w, 2, 0, -16 );
  w = xge_NewSwitch ( win0, w, swM01COMMAND, 16, 16, 20, xge_HEIGHT-16,
                      txtNull, &win0commandline );
  xge_SetWidgetPositioning ( w, 2, 20, -16 );
  side00ewidgets = w;
} /*InitSide00eMenu*/

boolean VerifyNumParam ( char *txt, double *param )
{
  double x;

  if ( sscanf ( txt, "%le", &x ) == 1 ) {
    *param = x;
    return true;
  }
  else
    return false;
} /*VerifyNumParam*/

boolean VerifyAngleParam ( char *txt, double *param )
{
  double x;

  if ( pkv_DegreeStrToRad ( txt, &x ) ) {
    *param = x;
    return true;
  }
  else
    return false;
} /*VerifyAngleParam*/

void EscapeNumParam ( char *txt, double param )
{
  sprintf ( txt, "%g", param );
} /*EscapeNumParam*/

void EscapeAngleParam ( char *txt, double param )
{
  pkv_RadToDegreeStr ( param, txt );
} /*EscapeAngleParam*/

int Side00eMenuCallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnM01eTRANSLATE:
      GeomObjectSpecial3DTranslate ( (vector3d*)&side00eparam[3] );
      xge_Redraw ();
      return 1;
  case btnM01eSCALE:
      GeomObjectSpecial3DScale ( (point3d*)&side00eparam[0],
                                 (vector3d*)&side00eparam[3] );
      xge_Redraw ();
      return 1;
  case btnM01eROTATE:
      GeomObjectSpecial3DRotate ( (point3d*)&side00eparam[0],   
                                  (vector3d*)&side00eparam[3],
                                  side00eparam[7] );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_LINE:
      GeomObjectProject3DPointsOnALine ( (point3d*)&side00eparam[0],
                                         (vector3d*)&side00eparam[3] );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_PLANE:
      GeomObjectProject3DPointsOnAPlane ( (point3d*)&side00eparam[0],
                                          (vector3d*)&side00eparam[3], 0 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_PLANE_UP:
      GeomObjectProject3DPointsOnAPlane ( (point3d*)&side00eparam[0],
                                          (vector3d*)&side00eparam[3], 2 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_PLANE_DOWN:
      GeomObjectProject3DPointsOnAPlane ( (point3d*)&side00eparam[0],
                                          (vector3d*)&side00eparam[3], 1 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_SPHERE:
      GeomObjectProject3DPointsOnASphere ( (point3d*)&side00eparam[0],
                                           side00eparam[6], 0 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_SPHERE_UP:
      GeomObjectProject3DPointsOnASphere ( (point3d*)&side00eparam[0],
                                           side00eparam[6], 2 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_SPHERE_DOWN:
      GeomObjectProject3DPointsOnASphere ( (point3d*)&side00eparam[0],
                                           side00eparam[6], 1 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_CYLINDER:
      GeomObjectProject3DPointsOnACylinder ( (point3d*)&side00eparam[0],
                                             (vector3d*)&side00eparam[3],
                                             side00eparam[6], 0 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_CYLINDER_UP:
      GeomObjectProject3DPointsOnACylinder ( (point3d*)&side00eparam[0],
                                             (vector3d*)&side00eparam[3],
                                             side00eparam[6], 2 );
      xge_Redraw ();
      return 1;
  case btnM01ePROJ_CYLINDER_DOWN:
      GeomObjectProject3DPointsOnACylinder ( (point3d*)&side00eparam[0],
                                             (vector3d*)&side00eparam[3],
                                             side00eparam[6], 1 );
      xge_Redraw ();
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case textedM01ePX:
      return VerifyNumParam ( side00eparam_str[0], &side00eparam[0] );
  case textedM01ePY:
      return VerifyNumParam ( side00eparam_str[1], &side00eparam[1] );
  case textedM01ePZ:
      return VerifyNumParam ( side00eparam_str[2], &side00eparam[2] );
  case textedM01eVX:
      return VerifyNumParam ( side00eparam_str[3], &side00eparam[3] );
  case textedM01eVY:
      return VerifyNumParam ( side00eparam_str[4], &side00eparam[4] );
  case textedM01eVZ:
      return VerifyNumParam ( side00eparam_str[5], &side00eparam[5] );
  case textedM01eR:
      return VerifyNumParam ( side00eparam_str[6], &side00eparam[6] );
  case textedM01eANGLE:
      return VerifyAngleParam ( side00eparam_str[7], &side00eparam[7] );
  default:
      return 1;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case textedM01ePX:
      EscapeNumParam ( side00eparam_str[0], side00eparam[0] );
      return 1;
  case textedM01ePY:
      EscapeNumParam ( side00eparam_str[1], side00eparam[1] );
      return 1;
  case textedM01ePZ:
      EscapeNumParam ( side00eparam_str[2], side00eparam[2] );
      return 1;
  case textedM01eVX:
      EscapeNumParam ( side00eparam_str[3], side00eparam[3] );
      return 1;
  case textedM01eVY:
      EscapeNumParam ( side00eparam_str[4], side00eparam[4] );
      return 1;
  case textedM01eVZ:
      EscapeNumParam ( side00eparam_str[5], side00eparam[5] );
      return 1;
  case textedM01eR:
      EscapeNumParam ( side00eparam_str[6], side00eparam[6] );
      return 1;
  case textedM01eANGLE:
      EscapeAngleParam ( side00eparam_str[7], side00eparam[7] );
      return 1;
  default:
      return 1;
    }

case xgemsg_RESIZE:
    side00escroll->msgproc ( side00escroll, msg, key, x, y-40 );
    return 1;

default:
    return 0;
  }
} /*Side00eMenuCallBack*/

void EnterRefLine ( point3d *p )
{
  int      i;

  SubtractPoints3d ( &p[1], &p[0], (vector3d*)&side00eparam[3] );
  memcpy ( &side00eparam[0], &p[0], sizeof(point3d) );
  for ( i = 0; i < 6; i++ )
    EscapeNumParam ( side00eparam_str[i], side00eparam[i] );
} /*EnterRefLine*/

