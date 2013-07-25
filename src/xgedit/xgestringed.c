
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "xgedit.h"
#include "xgeprivate.h"

static boolean insert = true;
static boolean storebytes = false;

static void CorrectCursorPos ( xge_string_ed *ed )
{
  if ( ed->pos < 0 )
    ed->pos = 0;
  else if ( ed->pos > ed->maxlength )
    ed->pos = ed->maxlength;
  if ( ed->pos < ed->start )
    ed->start = ed->pos;
  else if ( ed->pos >= ed->start+ed->chdisp )
    ed->start = (short)(ed->pos-ed->chdisp+1);
} /*CorrectCursorPos*/

void xge_DrawStringEd ( xge_widget *er, boolean onscreen )
{
  xge_string_ed *ed;
  char          *text;
  char          buffer[xge_MAX_STRING_LENGTH+1];
  short         x, h, lgt;

  text = er->data0;
  ed   = er->data1;
  lgt  = (short)(strlen(text)-ed->start);
  lgt = (short)min(lgt, xge_MAX_STRING_LENGTH);
  if ( er->state == xgestate_TEXT_EDITING ) {
    xgeSetForeground ( xgec_Blue6 );
    xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
    CorrectCursorPos ( ed );
        /* draw the text cursor */
    xgeSetForeground ( xgec_DodgerBlue );
    x = (short)(er->x + 2 + (ed->pos-ed->start)*xge_CHAR_WIDTH);
    if ( insert )
      h = (xge_CHAR_HEIGHT+1) / 2;
    else
      h = xge_CHAR_HEIGHT / 4;
    xgeFillRectangle ( xge_CHAR_WIDTH, h, x, er->y+er->h-2-h );
  }
  else {
    xgeSetForeground ( xgec_Blue3 );
    xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  }

  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( lgt > 0 ) {
    memcpy ( buffer, &text[ed->start], lgt );
    buffer[lgt] = 0;
    xgeDrawString ( buffer, er->x+2, er->y+er->h-5 );
  }
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*xge_DrawStringEd*/

boolean xge_StringEdMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  xge_string_ed *ed;
  char          *text, *buf;
  short         i, j, lgt;
  int           nb;

  text = er->data0;
  ed   = er->data1;
  lgt  = (short)strlen(text);
  if ( er->state == xgestate_TEXT_EDITING ) {
    switch ( msg ) {
  case xgemsg_SPECIAL_KEY:
      switch ( key ) {
    case 0xFF51:  /* arrow left */
        ed->pos--;
        goto redraw_string;
    case 0xFF53:  /* arrow right */
        ed->pos++;
        goto redraw_string;
    case 0xFF95:  /* Home - skladak */
    case 0xFF50:  /* Home - notebook Toshiba */
        ed->pos = 0;
        goto redraw_string;
    case 0xFF9c:  /* End - skladak */
    case 0xFF57:  /* End - notebook Toshiba */
        ed->pos = lgt;
        goto redraw_string;
    case 0xFF9f:  /* Del - skladak */
        goto delete_char;
    case 0xFE20:  /* Shift+Tab - notebook Toshiba */
        break;
    case 0xFF63:  /* insert - notebook Toshiba */
    case 0xFF9e:  /* insert - Dell */
        insert = (boolean)(!insert);
        goto redraw_string;
    default:
/*printf ( "%x\n", key );*/
        break;
      }
      break;

  case xgemsg_KEY:
/*printf ( "+key: %x, %x\n", key, (int)xgekeysym );*/
      switch ( key ) {
  case 0x08:  /* Backspace */
        if ( ed->pos > lgt ) {
          ed->pos = lgt;
          goto redraw_string;
        }
        else if ( ed->pos > 0 ) {
          ed->pos --;
          goto delete_char;
        }
        else
          break;
    case 0x7F:  /* Del - notebook Toshiba */
delete_char:
        if ( ed->pos < lgt ) {
          for ( i = ed->pos;  (text[i] = text[i+1]);  i++ )
            ;
          goto redraw_string;
        }
        else
          break;
  case 0x09:  /* Tab */
  case 0x0D:
        goto exit_editing_mode;
  case 0x1B:  /* Esc */
        if ( xge_callback ( er, xgemsg_TEXT_EDIT_ESCAPE, 0, 0, 0 ) == 1 ) {
                /* stop editing, without alerting the application */
          er->state = xgestate_NOTHING;
          xge_ReleaseFocus ( er );
          xge_SetClipping ( er );
          er->redraw ( er, true );
          xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
        }
        else
          break;
  default:
        if ( key >= 0x20 && key <= 0x7F ) { /* enter a character */
          if ( ed->pos > lgt )
            ed->pos = lgt;
          if ( insert && lgt >= ed->maxlength && ed->pos < ed->maxlength )
            text[--lgt] = 0;
          if ( insert && ed->pos < lgt ) {
            for ( i = lgt; i >= ed->pos; i-- )
              text[i+1] = text[i];
            text[ed->pos++] = (char)key;
          }
          else if ( ed->pos < ed->maxlength ) {
            text[ed->pos++] = (char)key;
            if ( ed->pos == lgt+1 )
              text[ed->pos] = 0;
          }
redraw_string:
          CorrectCursorPos ( ed );
          xge_SetClipping ( er );
          er->redraw ( er, true );
          break;
        }
        break;
      }
      break;

  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_CHANGE) &&
           (key & xgemouse_LBUTTON_DOWN) ) {
        if ( !xge_PointInRect ( er, x, y ) ) {
exit_editing_mode:
          if ( xge_callback ( er, xgemsg_TEXT_EDIT_VERIFY, 0, 0, 0 ) ) {
            er->state = xgestate_NOTHING;
            xge_ReleaseFocus ( er );
            xge_SetClipping ( er );
            er->redraw ( er, true );
            xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
            xge_callback ( er, xgemsg_TEXT_EDIT_ENTER, 0, 0, 0 );
          }
          storebytes = false;
        }
        else {  /* move the text cursor position */
          ed->pos = (short)(ed->start + (x-er->x-2) / xge_CHAR_WIDTH);
          CorrectCursorPos ( ed );
          xge_SetClipping ( er );
          er->redraw ( er, true );
          storebytes = true;
        }
      }
      else if ( er->state == xgestate_TEXT_EDITING &&
                (key & xgemouse_LBUTTON_CHANGE) &&
                !(key & xgemouse_LBUTTON_DOWN) ) {
        nb = (short)(ed->start + (x-er->x-2) / xge_CHAR_WIDTH);
        if ( nb > lgt ) nb = lgt;
        nb -= ed->start;
        if ( storebytes && nb > ed->pos ) {
          XStoreBytes ( xgedisplay, &text[ed->pos], nb-ed->pos );
          storebytes = false;
        }
      }
      else if ( (key & xgemouse_RBUTTON_CHANGE) &&
                (key & xgemouse_RBUTTON_DOWN) ) {
        buf = XFetchBytes ( xgedisplay, &nb );
                /* type in the characters */
        if ( ed->pos > lgt )
          ed->pos = lgt;
        for ( i = 0;  i < nb && buf[i] >= 0x20 && !(buf[i] &0x80);  i++ ) {
          if ( insert ) {
            if ( lgt >= ed->maxlength && ed->pos < ed->maxlength )
              text[--lgt] = 0;
            if ( ed->pos < lgt ) {
              for ( j = lgt; j >= ed->pos; j-- )
                text[j+1] = text[j];
            }
            text[ed->pos++] = buf[i];
            lgt ++;
          }
          else {
            text[ed->pos++] = buf[i];
            if ( ed->pos == lgt+1 ) {
              text[ed->pos] = 0;
              lgt ++;
            }
          }
          if ( ed->pos >= ed->maxlength )
            break;
        }
        goto redraw_string;
      }
      break;

  default:
      break;
    }
  }
  else if ( er->state == xgestate_NOTHING ) {
    switch ( msg ) {
  case xgemsg_KEY:
/*printf ( "-key: %x, %x\n", key, (int)xgekeysym );*/
      if ( key == 0x0D )
        goto enter_editing_mode;
      break;

  case xgemsg_MCLICK:
      if ( (key & xgemouse_LBUTTON_CHANGE) &&
           (key & xgemouse_LBUTTON_DOWN) ) {
enter_editing_mode:
        er->state = xgestate_TEXT_EDITING;
        xge_GrabFocus ( er, false );
        xge_SetClipping ( er );
        er->redraw ( er, true );
        xge_SetCurrentWindowCursor ( xgeCURSOR_PENCIL );
        storebytes = false;
      }
      break;

  default:
      break;
    }
  }

  return true;
} /*xge_StringEdMsg*/

xge_widget *xge_NewStringEd ( char window_num, xge_widget *prev, int id,
                              short w, short h, short x, short y,
                              short maxlength, char *text, xge_string_ed *ed )
{
  xge_widget *er;

  er = xge_NewWidget ( window_num, prev, id, w, h, x, y, text, ed,
                       xge_StringEdMsg, xge_DrawStringEd );
  ed->er = er;
  ed->maxlength = maxlength;
  ed->chdisp = (short)((w-4) / xge_CHAR_WIDTH);
  ed->start = 0;
  ed->pos = 0;
  return er;
} /*xge_NewStringEd*/

