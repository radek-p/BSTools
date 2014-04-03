
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


xge_widget *InitPopup01 ( void )
{
  xge_widget *w, *menu;

  w = xge_NewTextWidget ( win0, NULL, txtP01DIRSTR, 380, 16, 20+10, 40+10,
                          current_dir );
  w = xge_NewSwitch ( win0, w, swP01HIDDEN, 54, 16, 20+10, 40+30,
                      txtHidden, &showhiddenfiles );
  w = xge_NewButton ( win0, w, btnP01OPEN, 58, 19, 20+91, 40+232-30, txtOpen );
  w = xge_NewButton ( win0, w, btnP01CANCEL, 58, 19, 20+251, 40+232-30, txtCancel );
  w = xge_NewListBox ( win0, w, lbP01DIRLIST, 180, 131, 20+10, 40+60, &dirlist1 );
  w = xge_NewListBox ( win0, w, lbP01FILELIST, 180, 131, 220+10, 40+60, &filelist1 );
  menu = xge_NewFMenu ( win0, NULL, POPUP01, 400, 232, 20, 40, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup01*/

void Popup01ChangeDir ( void )
{
  if ( !chdir ( &dirlist1.itemstr[dirlist1.itemind[dirlist1.currentitem]] ) ) {
    PreparePopup01 ();
    xge_SetClipping ( popup01 );
    popup01->redraw ( popup01, true );
  }
} /*Popup01ChangeDir*/

void Popup01ChangeDirAlt ( short x )
{
  short l1, l2;

  l1 = strlen ( current_directory );
  l2 = strlen ( current_dir );
  x += l2-l1;
  while ( x < l1 && current_directory[x] != '/' )
    x++;
  if ( x < l1 ) {
    current_directory[x] = 0;
    if ( chdir ( current_directory ) )
      current_directory[x] = '/';
    else {
      PreparePopup01 ();
      xge_SetClipping ( popup01 );
      popup01->redraw ( popup01, true );
    }
  }
} /*Popup01ChangeDirAlt*/

void Popup01CameraReader ( void *usrdata, int ident, CameraRecd *camera )
{
  if ( !camera->parallel ) {  /* now only the perspective camera */
    g00win3D.CPos[4] = g00win3D.CPos[3];  /* save the last projection */
    g00win3D.CPos[3].position = camera->position;
    g00win3D.CPos[3].psi      = camera->psi;
    g00win3D.CPos[3].theta    = camera->theta;
    g00win3D.CPos[3].phi      = camera->phi;
    g00win3D.CPos[3].zmin     = camera->zmin;
    g00win3D.CPos[3].zmax     = camera->zmax;
    g00win3D.CPos[3].vd.persp.f    = camera->vd.persp.f;
    g00win3D.CPos[3].vd.persp.xi0  = camera->vd.persp.xi0;
    g00win3D.CPos[3].vd.persp.eta0 = camera->vd.persp.eta0;
    if ( !CameraSetMappingd ( &g00win3D.CPos[3] ) )
      g00win3D.CPos[3] = g00win3D.CPos[4];  /* restore if wrong data */
  }
} /*Popup01CameraReader*/

void Popup01OpenFile ( void )
{
  geom_object *save_current;

  save_current = current_go;
  if ( !GeomObjectReadFile ( filename, Popup01CameraReader ) )
    xge_DisplayErrorMessage ( ErrorMsgCannotOpen, 0 );
  if ( save_current )
    current_go = save_current;
  else
    current_go = first_go;
  SetupObjectSpecificMenus ( current_go );
  if ( xge_IsPopupOn ( popup12 ) )
    UpdateObjectNameList ();
  xge_RedrawAll ();
} /*Popup01OpenFile*/

int Popup01CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP01OPEN:
      if ( xge_GetCurrentListBoxString ( &filelist1, filename ) ) {
        xge_RemovePopup ( true );
        Popup01OpenFile ();
      }
      return 1;
  case btnP01CANCEL:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swP01HIDDEN:
      PreparePopup01 ();
      xge_SetClipping ( popup01 );
      popup01->redraw ( popup01, true );
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_WIDGET_CLICK:
    switch ( er->id ) {
  case txtP01DIRSTR:
      Popup01ChangeDirAlt ( x );
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lbP01DIRLIST:
      Popup01ChangeDir ();
      return 1;
  case lbP01FILELIST:
      if ( xge_GetCurrentListBoxString ( &filelist1, filename ) ) {
        xge_RemovePopup ( true );
        Popup01OpenFile ();
      }
      return 1;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup01CallBack*/

