
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
#include "render.h"

#define PARENT_SIDE
#include "pozwalajipc.h"
#undef PARENT_SIDE

static boolean swActive = false;
static xge_widget *wdgActive;

xge_widget *InitPopup12 ( void )
{
  xge_widget *w, *menu;

  memset ( &objectlist, 0, sizeof(xge_listbox) );
  w = xge_NewListBox ( win1, NULL, lbP12OBJLIST, 180, 99, 20+10, 40+10,
                       &objectlist );
  objectlist.bk0 = xgec_Blue6;
  objectlist.bk1 = xgec_Blue3;
  w = xge_NewButton ( win1, w, btnP12NEW, 58, 19, 20+200, 40+10, txtNew );
  w = xge_NewButton ( win1, w, btnP12COPY, 58, 19, 20+200, 40+34, txtCopy );
  w = xge_NewButton ( win1, w, btnP12DELETE, 58, 19, 20+200, 40+58, txtDelete );
  w = xge_NewButton ( win1, w, btnP12PURGE, 58, 19, 20+200, 40+82, txtPurge );
  wdgActive = w = xge_NewSwitch ( win1, w, swP12SHOW, 58, 16, 20+200, 40+106,
                                  txtActive, &swActive );
  w = xge_NewButton ( win1, w, btnP12OK, 58, 19, 20+105, 40+160-30, txtOK );
  menu = xge_NewFMenu ( win1, NULL, POPUP12, 268, 160, 20, 40, w );
  menu->msgproc = xge_PopupMenuMsg;
  return menu;
} /*InitPopup12*/

boolean SetupObjectNameList ( void )
{
  int         lsum, nit, lgt;
  geom_object *go;

  nit = GeomObjectNumber ();
  objectlist.itemind = malloc ( (nit+1)*sizeof(int) );
  if ( !objectlist.itemind ) {
    DeleteObjectNameList ();
    return false;
  }
  objectlist.currentitem = 0;
  for ( go = first_go, lsum = 0, nit = 0;  go;  go = go->next, nit ++ ) {
    lgt = strlen ( go->name ) + 1;
    objectlist.itemind[nit] = lsum;
    lsum += lgt;
    if ( go == current_go )
      objectlist.currentitem = nit;
  }
  objectlist.itemstr = malloc ( (lsum+1)*sizeof(char) );
  if ( !objectlist.itemstr ) {
    DeleteObjectNameList ();
    return false;
  }
  for ( go = first_go, nit = 0;  go;  go = go->next, nit ++ )
    strcpy ( &objectlist.itemstr[objectlist.itemind[nit]], go->name );
  objectlist.nitems = nit;
  if ( nit <= objectlist.dlistnpos )
    objectlist.fditem = 0;
  else if ( objectlist.currentitem < objectlist.fditem )
    objectlist.fditem = objectlist.currentitem;
  else if ( objectlist.currentitem-objectlist.fditem >= objectlist.dlistnpos )
    objectlist.fditem = objectlist.currentitem-objectlist.dlistnpos+1;
  if ( current_go )
    swActive = current_go->active;
  else
    swActive = false;
  return true;
} /*SetupObjectNameList*/

void DeleteObjectNameList ( void )
{
  if ( objectlist.itemind ) free ( objectlist.itemind );
  if ( objectlist.itemstr ) free ( objectlist.itemstr );
  objectlist.nitems = objectlist.fditem = objectlist.currentitem = 0;
} /*DeleteObjectNameList*/

void UpdateObjectNameList ( void )
{
  DeleteObjectNameList ();
  SetupObjectNameList ();
} /*UpdateObjectNameList*/

void CleanupPopup12 ( void )
{
  DeleteObjectNameList ();
  SetupObjectSpecificMenus ( current_go );
  if ( RenderingIsOn )
    StopRendering ();
  rendered_picture = false;
  xge_RedrawAll ();
} /*CleanupPopup12*/

int Popup12CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int sc;

  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP12NEW:
      OpenPopup ( popup10, false );
      return 1;
  case btnP12COPY:
      if ( GeomObjectCopyCurrent () ) {
        UpdateObjectNameList ();
        SetupObjectSpecificMenus ( current_go );
        xge_SetClipping ( popup12 );
        popup12->redraw ( popup12, true );
        xge_ResetClipping ();
      }
      return 1;
  case btnP12DELETE:
      GeomObjectDeleteCurrent ();
      UpdateObjectNameList ();
      SetupObjectSpecificMenus ( current_go );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  case btnP12PURGE:
      if ( ipc_state == ipcstate_CHILD_BUSY )
        IPCInterruptTheChild ();
      GeomObjectPurgeList ();
      UpdateObjectNameList ();
      SetupObjectSpecificMenus ( current_go );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  case btnP12OK:
      xge_RemovePopup ( true );
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case swP12SHOW:
      current_go->active = swActive;
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_VERIFY:
    switch ( er->id ) {
  case txtedP12OBJNAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ENTER:
    switch ( er->id ) {
  case txtedP12OBJNAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_EDIT_ESCAPE:
    switch ( er->id ) {
  case txtedP12OBJNAME:
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_SET:
    switch ( er->id ) {
  case lbP12OBJLIST:
      GeomObjectSelect ( objectlist.currentitem );
      swActive = current_go->active;
      switch ( current_go->cpdimen ) {
    case 2:
        SetWin002D ();
        break;
    case 3:
        SetWin003D ();
        if ( RenderingIsOn )
          StopRendering ();
        rendered_picture = false;
        break;
    default:
        break;
      }
      xge_RedrawAll ();
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lbP12OBJLIST:
      SetupObjectSpecificMenus ( current_go );
      if ( RenderingIsOn )
        StopRendering ();
      rendered_picture = false;
      xge_RedrawAll ();
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_EXCHANGE:
    switch ( er->id ) {
  case lbP12OBJLIST:
      if ( GeomObjectExchangeInTheList ( key == -1 ) ) {
        sc = objectlist.currentitem;
        UpdateObjectNameList ();
        objectlist.currentitem = sc;
        return 1;
      }
      else
        return 0;
  default:
      return 0;
    }

default:
    return 0;
  }
} /*Popup12CallBack*/

