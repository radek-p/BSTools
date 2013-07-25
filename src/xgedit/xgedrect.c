
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2012                            */
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

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

#define TABSIZE 32

typedef struct edr_table {
    struct edr_table *next;
    xge_widget edr[TABSIZE];
  } edr_table, *edr_tablep;

static edr_tablep rect_list = NULL;
static int rect_num = TABSIZE;

boolean xge_InitRectAllocation ( void )
{
  rect_list = NULL;
  rect_num = TABSIZE;
  return true;
} /*xge_InitRectAllocation*/

xge_widget *xge_AllocEdRect ( void )
{
  edr_tablep  next;
  xge_widget *new_rect;

  if ( rect_num >= TABSIZE ) {
    if ( !(next = malloc ( sizeof(edr_table) )) )
      return NULL;
    next->next = rect_list;
    rect_list = next;
    rect_num = 0;
  }
  new_rect = &(rect_list->edr[rect_num++]);
  return new_rect;
} /*xge_AllocEdRect*/

void xge_FreeEdRectangles ( void )
{
  edr_tablep next;

  while ( (next = rect_list) ) {
    rect_list = next->next;
    free ( next );
  }
} /*xge_FreeEdRectangles*/

/* ///////////////////////////////////////////////////////////////////////// */
xge_widget *xge_NewWidget (
         char window_num, xge_widget *prev, int id,
         short w, short h, short x, short y,
         void *data0, void *data1,
         boolean (*msgproc)(xge_widget*, int, int, short, short),
         void (*redraw)(xge_widget*, boolean) )
{
  xge_widget *new_rect;

  if ( (new_rect = xge_AllocEdRect () ) ) {
    new_rect->id = id;
    new_rect->w = w;
    new_rect->h = h;
    new_rect->x = new_rect->xofs = x;
    new_rect->y = new_rect->yofs = y;
    new_rect->rpos = -1;      /* the default value - do not reposition */
    new_rect->window_num = window_num;
    new_rect->data0 = data0;
    new_rect->data1 = data1;
    new_rect->data2 = NULL;
    new_rect->msgproc = msgproc;
    new_rect->redraw = redraw;
    new_rect->next = NULL;
    new_rect->prev = prev;
    new_rect->up = NULL;
    new_rect->state = xgestate_NOTHING;  /* default, 0 */
    if ( prev )
      prev->next = new_rect;
  }
  return new_rect;
} /*xge_NewWidget*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_SetupEdRect ( char window_num, xge_widget *edr, int en, int id,
                       short w, short h, short x, short y,
                       boolean (*msgproc) ( xge_widget*, int, int, short, short ),
                       void (*redraw) ( xge_widget*, boolean ) )
{
  edr[en].id = id;
  edr[en].w = w;
  edr[en].h = h;
  edr[en].x = edr[en].xofs = x;
  edr[en].y = edr[en].yofs = y;
  edr[en].rpos = -1;
  edr[en].window_num = window_num;
  edr[en].msgproc = msgproc;
  edr[en].redraw = redraw;
} /*xge_SetupEdRect*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_SetWidgetPositioning ( xge_widget *edr,
                                char rpos, short xofs, short yofs )
{
  edr->rpos = rpos;
  edr->xofs = xofs;
  edr->yofs = yofs;
} /*xge_SetWidgetPositioning*/

