
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdio.h>
#include <malloc.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "xgedit.h"

#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
void xge_ShortenString ( const char *s, char *buf, int maxlen )
{
  int len, l;
  char dots[] = "...";

  len = strlen ( s );
  if ( len <= maxlen ) {
    memcpy ( buf, s, len+1 );
  }
  else if (maxlen >= 5 ) {
    l = (maxlen-3)/2;
    memcpy ( buf, s, l );
    memcpy ( &buf[l], dots, 3 );
    memcpy ( &buf[l+3], &s[len-l], l );
    buf[maxlen] = 0;
  }
  else
    buf[0] = 0;
} /*xge_ShortenString*/

void xge_DrawListBox ( xge_widget *er, boolean onscreen )
{
  int         i, p0, p1;
  xge_listbox *lbox;
  char        buf[120];

  lbox = er->data0;
  xgeSetForeground ( xgec_Grey6 );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( lbox->nitems ) {
    p0 = lbox->fditem;
    p1 = p0 + lbox->dlistnpos;  p1 = min ( p1, lbox->nitems );
    for ( i = p0; i < p1; i++ ) {
      if ( i == lbox->currentitem )
        xgeSetForeground ( lbox->bk1 );
      else
        xgeSetForeground ( lbox->bk0 );
      xgeFillRectangle ( er->w-2, xge_LISTDIST-1,
                         er->x+1, er->y+2+(i-lbox->fditem)*xge_LISTDIST );
      xgeSetForeground ( xgec_White );
      xge_ShortenString ( &lbox->itemstr[lbox->itemind[i]], buf, lbox->maxitl );
      xgeDrawString ( buf, er->x+2, er->y+3+(i-lbox->fditem+1)*xge_LISTDIST-5 );
    }
  }

  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawListBox*/

boolean xge_ListBoxMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_listbox *lbox;
  int npos;

  lbox = er->data0;
  switch ( msg ) {
case xgemsg_ENTERING:
    if ( er->state != xgestate_LISTBOX_PICKING ) {
                   /* grab the keyboard input */
      xge_GrabFocus ( er, false );
      er->state = xgestate_LISTBOX_PICKING;
      xge_SetCurrentWindowCursor ( xgeCURSOR_HAND );
    }
    break;

case xgemsg_EXITING:
    if ( er->state == xgestate_LISTBOX_PICKING ) {
release_input:
                   /* release the keyboard input */
      xge_ReleaseFocus ( er );
      er->state = xgestate_NOTHING;
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    }
    break;

case xgemsg_MMOVE:
    if ( !xge_PointInRect ( er, x, y ) && er->state == xgestate_LISTBOX_PICKING )
      goto release_input;
    return false;

case xgemsg_MCLICK:
    if ( er->state != xgestate_LISTBOX_PICKING )
      return false;
    if ( !xge_PointInRect ( er, x, y ) )
      goto release_input;
    if ( (key & xgemouse_LBUTTON_DOWN) && (key && xgemouse_LBUTTON_CHANGE) ) {
      npos = (y-er->y-2)/xge_LISTDIST;
      npos = max ( npos, 0 );
      npos = min ( npos, lbox->dlistnpos );
      npos += lbox->fditem;
      if ( npos < lbox->nitems ) {
        if ( npos == lbox->currentitem ) {
          xge_callback ( er, xgemsg_LISTBOX_ITEM_PICK, 0, 0, 0 );
        }
        else {
          lbox->currentitem = (short)npos;
          goto notify;
        }
      }
    }
    else if ( (key & xgemouse_WHEELFW_DOWN) && (key && xgemouse_WHEELFW_CHANGE) )
      goto move_up;
    else if ( (key & xgemouse_WHEELBK_DOWN) && (key && xgemouse_WHEELBK_CHANGE) )
      goto move_down;
    break;

case xgemsg_KEY:
    switch ( key ) {
  case 0x0D:   /* Enter */
      if ( lbox->nitems )  /* bother only if the box is nonempty */
        xge_callback ( er, xgemsg_LISTBOX_ITEM_PICK, 0, 0, 0 );
      break;
  case 'X':    /* request to move the current item one position up */
      if ( lbox->nitems > 1 && lbox->currentitem > 0 ) {
        if ( xge_callback ( er, xgemsg_LISTBOX_EXCHANGE, -1, 0, 0 ) ) {
          lbox->currentitem --;
          goto notify;
        }
      }
      break;
  case 'x':    /* request to move the current item one position down */
      if ( lbox->nitems > 1 && lbox->currentitem < lbox->nitems-1 ) {
        if ( xge_callback ( er, xgemsg_LISTBOX_EXCHANGE, +1, 0, 0 ) ) {
          lbox->currentitem ++;
          goto notify;
        }
      }
      break;
  default:
      return false;
    }
    break;

case xgemsg_SPECIAL_KEY:
    switch ( key ) {
  case 0xFF95:  /* Home */
  case 0xFF50:  /* Home - notebook */
      if ( lbox->currentitem > 0 ) {
        lbox->currentitem = lbox->fditem = 0;
        goto notify;
      }
      break;
  case 0xFF9C:  /* End */
  case 0xFF57:  /* End - notebook */
      if ( lbox->currentitem < lbox->nitems-1 ) {
        lbox->currentitem = (short)(lbox->nitems-1);
        lbox->fditem = (short)(max ( lbox->fditem, lbox->nitems-lbox->dlistnpos));
        goto notify;
      }
      break;
  case 0xFF52:  /* arrow up */
move_up:
      if ( lbox->currentitem > 0 ) {
        lbox->currentitem --;
        if ( lbox->currentitem < lbox->fditem )
          lbox->fditem = lbox->currentitem;
        goto notify;
      }
      break;
  case 0xFF54:  /* arrow down */
move_down:
      if ( lbox->currentitem < lbox->nitems-1 ) {
        lbox->currentitem ++;
        if ( lbox->currentitem >= lbox->fditem+lbox->dlistnpos )
          lbox->fditem ++;
        goto notify;
      }
      break;
  case 0xFF9A:  /* PgUp - Vobis */
  case 0xFF55:  /* PgUp - notebook */
      if ( lbox->currentitem > 0 ) {
        lbox->currentitem = (short)(lbox->currentitem-lbox->dlistnpos);
        lbox->currentitem = (short)(max ( 0, lbox->currentitem ));
        lbox->fditem = (short)(lbox->fditem-lbox->dlistnpos);
        lbox->fditem = (short)(max( 0, lbox->fditem ));
        goto notify;
      }
      break;
  case 0xFF9B:  /* PgDn - Vobis */
  case 0xFF56:  /* PgDn - notebook */
      if ( lbox->currentitem < lbox->nitems-1 ) {
        npos = lbox->dlistnpos;
        if ( lbox->currentitem+npos > lbox->nitems-1 )
          npos = lbox->nitems-1-lbox->currentitem;
        lbox->currentitem = (short)(lbox->currentitem+npos);
        if ( lbox->currentitem-lbox->fditem >= lbox->dlistnpos )
          lbox->fditem = (short)(lbox->fditem+npos);
notify:
        xge_callback ( er, xgemsg_LISTBOX_ITEM_SET, 0, 0, 0 );
        xge_SetClipping ( er );
        er->redraw ( er, true );
      }
      break;
  default:
      return false;
    }
    break;

default:
    return false;
  }
  return true;
} /*xge_ListBoxMsg*/

xge_widget *xge_NewListBox ( char window_num, xge_widget *prev, int id,
                             short w, short h, short x, short y,
                             xge_listbox *listbox )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y, listbox, NULL,
                       xge_ListBoxMsg, xge_DrawListBox );
  if ( er ) {
    listbox->er = er;
    listbox->dlistnpos = (char)((h-3)/xge_LISTDIST);
    listbox->maxitl = (char)((w-4)/xge_CHAR_WIDTH);
    listbox->nitems = listbox->currentitem = listbox->fditem = 0;
    listbox->itemind = NULL;
    listbox->itemstr = NULL;
        /* assign default background colours */
    listbox->bk0 = xgec_Chocolate4;
    listbox->bk1 = xgec_DarkOrange2;
        /* other fields of this structure are undefined; */
        /* the application must take care of that. */
  }
  return er;
} /*xge_NewListBox*/

void xge_ClearListBox ( xge_listbox *lbox )
{
  if ( lbox->itemind ) {
    free ( lbox->itemind );
    lbox->itemind = NULL;
  }
  if ( lbox->itemstr ) {
    free ( lbox->itemstr );
    lbox->itemstr = NULL;
  }
  lbox->nitems = lbox->currentitem = lbox->fditem = 0;
} /*xge_ClearListBox*/

boolean xge_GetCurrentListBoxString ( xge_listbox *lbox, char *string )
{
  if ( lbox->currentitem >= 0 && lbox->currentitem < lbox->nitems ) {
    strcpy ( string, &lbox->itemstr[lbox->itemind[lbox->currentitem]] );
    return true;
  }
  else {
    string[0] = 0;
    return false;
  }
} /*xge_GetCurrentListBoxString*/

int xge_MoveInListBox ( xge_listbox *lbox, short amount )
{
  lbox->currentitem += amount;
  if ( lbox->currentitem < 0 )
    lbox->currentitem = 0;
  else if ( lbox->currentitem >= lbox->nitems )
    lbox->currentitem = lbox->nitems-1;
  if ( lbox->currentitem - lbox->fditem >= lbox->dlistnpos )
    lbox->fditem = lbox->currentitem-lbox->dlistnpos+1;
  else if ( lbox->currentitem < lbox->fditem )
    lbox->fditem = lbox->currentitem;
  return lbox->currentitem;
} /*xge_MoveInListBox*/

