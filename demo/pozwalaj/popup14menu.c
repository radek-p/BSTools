
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


xge_widget *InitPopup14 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewTextWidget ( win1, NULL, 0, 144, 16, 1+10, 81+10,
                          txtPreTransformation );
  w = xge_NewStringEd ( win1, w, txtedP14_A11, 81, 19, 1+10, 81+30,
                        MAX_PARAM_LGT, pretrans_str[0], &pretrans_ed[0] );
  w = xge_NewStringEd ( win1, w, txtedP14_A12, 81, 19, 1+94, 81+30,
                        MAX_PARAM_LGT, pretrans_str[1], &pretrans_ed[1] );
  w = xge_NewStringEd ( win1, w, txtedP14_A13, 81, 19, 1+178, 81+30,
                        MAX_PARAM_LGT, pretrans_str[2], &pretrans_ed[2] );
  w = xge_NewStringEd ( win1, w, txtedP14_A14, 81, 19, 1+262, 81+30,
                        MAX_PARAM_LGT, pretrans_str[3], &pretrans_ed[3] );
  w = xge_NewStringEd ( win1, w, txtedP14_A21, 81, 19, 1+10, 81+52,
                        MAX_PARAM_LGT, pretrans_str[4], &pretrans_ed[4] );
  w = xge_NewStringEd ( win1, w, txtedP14_A22, 81, 19, 1+94, 81+52,
                        MAX_PARAM_LGT, pretrans_str[5], &pretrans_ed[5] );
  w = xge_NewStringEd ( win1, w, txtedP14_A23, 81, 19, 1+178, 81+52,
                        MAX_PARAM_LGT, pretrans_str[6], &pretrans_ed[6] );
  w = xge_NewStringEd ( win1, w, txtedP14_A24, 81, 19, 1+262, 81+52,
                        MAX_PARAM_LGT, pretrans_str[7], &pretrans_ed[7] );
  w = xge_NewStringEd ( win1, w, txtedP14_A31, 81, 19, 1+10, 81+74,
                        MAX_PARAM_LGT, pretrans_str[8], &pretrans_ed[8] );
  w = xge_NewStringEd ( win1, w, txtedP14_A32, 81, 19, 1+94, 81+74,
                        MAX_PARAM_LGT, pretrans_str[9], &pretrans_ed[9] );
  w = xge_NewStringEd ( win1, w, txtedP14_A33, 81, 19, 1+178, 81+74,
                        MAX_PARAM_LGT, pretrans_str[10], &pretrans_ed[10] );
  w = xge_NewStringEd ( win1, w, txtedP14_A34, 81, 19, 1+262, 81+74,
                        MAX_PARAM_LGT, pretrans_str[11], &pretrans_ed[11] );
  w = xge_NewButton ( win1, w, btnP14_RESET, 60, 19, 1+10, 81+96,
                      txtReset );
  w = xge_NewButton ( win1, w, btnP14_OK, 60, 19, 1+147, 81+120, txtOK );
  menu = xge_NewFMenu ( win1, NULL, POPUP14, 354, 150, 1, 81, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup14*/

void GetPreTransformation ( void )
{
  double  *a;
  int     i;

  a = &current_go->pretrans.U0.a11;
  for ( i = 0; i < 12; i++ ) {
    sprintf ( pretrans_str[i], "%g", a[i] );
    pretrans_ed[i].start = pretrans_ed[i].pos = 0;
  }
} /*GetPreTransformation*/

void CleanupPopup14 ( void )
{
  current_go->display_pretrans = false;
  editing_pretrans = false;
  xge_SetWindow ( win0 );
  xge_Redraw ();
  xge_SetWindow ( win1 );
} /*CleanupPopup14*/

int Popup14CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  trans3d tr;
  int     ci;
  double  *c;

  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP14_RESET:
      IdentTrans3d ( &tr );
      GeomObjectSetPretransformation ( current_go, &tr );
      GetPreTransformation ();
      xge_SetClipping ( popup14 );
      popup14->redraw ( popup14, true );
      xge_SetWindow ( win0 );
      xge_Redraw ();
      return 1;
  case btnP14_OK:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
case xgemsg_TEXT_EDIT_ENTER:
    ci = er->id - txtedP14_A11;
    if ( ci >= 0 && ci < 12 ) {
      c = &current_go->pretrans.U0.a11 + ci;
      if ( VerifyNumParam ( pretrans_str[ci], c ) ) {
        xge_SetClipping ( geom00menu );
        geom00menu->redraw ( geom00menu, true );
        xge_SetWindow ( win1 );
        return 1;
      }
      else
        return 0;
    }
    return 1;

case xgemsg_TEXT_EDIT_ESCAPE:
    ci = er->id - txtedP14_A11;
    if ( ci >= 0 && ci < 12 ) {
      c = &current_go->pretrans.U0.a11 + ci;
      EscapeNumParam ( pretrans_str[ci], *c );
    }
    return 1;

default:
    return 0;
  }
} /*Popup14CallBack*/

