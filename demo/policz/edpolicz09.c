
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"

#define SWITCH_DIST 20

static int switch_nrows, switch_ncols;

static char switch_row_txt[MAX_CONSTR_SWITCH_ROWS][3] =
  {" 0", " 1"," 2"," 3"," 4"," 5"," 6"," 7",
   " 8"," 9","10","11","12","13","14","15"};

/* ///////////////////////////////////////////////////////////////////////// */
boolean ConstructConstraintSwitches ( void )
{
  int i;

  constr_sw[0] = xge_NewSwitch ( 0, NULL, sw01cCONSTR_SWITCH,
                             109, 16, 0, 0, txtCentralPoint, &constraint[0] );
  for ( i = 1; i < NUM_CONSTR_SWITCHES; i++ ) {
    constr_sw[i] = xge_NewSwitch ( 0, NULL, sw01cCONSTR_SWITCH+i,
                                   16, 16, 0, 0, NULL, &constraint[i] );
    if ( !constr_sw[i] )
      return false;
  }
  for ( i = 0; i < MAX_CONSTR_SWITCH_ROWS; i++ ) {
    constr_row_text[i] = xge_NewTextWidget ( 0, NULL, 0, 16, 16, 0, 0,
                                             switch_row_txt[i] );
    if ( !constr_row_text[i] )
      return false;
  }
  return true;
} /*ConstructConstraintSwitches*/

void ConfigureConstraintWidgets ( boolean reset )
{
  int        i, j, k, x, y, nconstr;
  xge_widget *prev, *v, *w;
  GHoptions  *opt;

        /* compute the number of available constraints */
  switch_nrows = hole_k;
  if ( constraints1 ) {
    opt = &options1;
    constr_surfno = 1;
  }
  else {
    opt = &options2;
    constr_surfno = 2;
  }
  if ( opt->order == 1 )
    switch_ncols = 2;
  else                /* it must be 2 */
    switch_ncols = 4;
  if ( opt->spline )
    switch_ncols += opt->nk*opt->m1;
  opt->nconstrsw = nconstr = switch_nrows*switch_ncols+1;
  if ( reset ) {
    opt->constr_type = 0;
    memset ( opt->constrsw, 0, nconstr*sizeof(boolean) );
  }
  memcpy ( constraint, opt->constrsw, nconstr*sizeof(boolean) );

        /* position the constraint switches in the scrolled menu */
  prev = constr_sw[0];
  prev->prev = NULL;
  xge_SetWidgetPositioning ( constr_sw[0], 0, 2, 2 );
  for ( i = 0, k = 1;  i < switch_nrows;  i++ ) {
    y = 2 + SWITCH_DIST*(i+1);
    constr_row_text[i]->prev = prev;
    xge_SetWidgetPositioning ( constr_row_text[i], 0, 5, y );
    prev = constr_row_text[i];
    for ( j = 0;  j < switch_ncols;  j++, k++ ) {
      x = 2+SWITCH_DIST*(j+1);
      constr_sw[k]->prev = prev;
      xge_SetWidgetPositioning ( constr_sw[k], 0, x, y );
      prev = constr_sw[k];
    }
  }  
  prev->next = NULL;
  for ( v = prev, w = v->prev;  w;  v = w, w = v->prev )
    w->next = v;
  menu01cscrolled->x = scroll_constr_sw.er->x;
  menu01cscrolled->y = scroll_constr_sw.er->y;
  menu01cscrolled->w = (switch_ncols+1)*SWITCH_DIST;
  menu01cscrolled->w = max ( menu01cscrolled->w, scroll_constr_sw.er->w );
  menu01cscrolled->h = (switch_nrows+1)*SWITCH_DIST;
  xge_SetMenuWidgets ( menu01cscrolled, prev, false );
  scroll_constr_sw.er->msgproc ( scroll_constr_sw.er, xgemsg_RESIZE, 0,
                        scroll_constr_sw.er->w, scroll_constr_sw.er->h );
  if ( !scroll_constr_sw.xslon && !scroll_constr_sw.yslon ) {
    for ( w = prev; w; w = w->prev ) {
      w->xofs -= 2;
      w->yofs -= 2;
    }
    xge_RepositionWidgets ( menu01cscrolled->w, menu01cscrolled->h,
                            menu01cscrolled->x, menu01cscrolled->y,
                            menu01cscrolled->data1 );
  }
} /*ConfigureConstraintWidgets*/

void SwitchTheConstraint ( int constrno )
{
  GHoptions *opt;
  boolean   *view_surf; 

  if ( constraints1 ) {
    opt = &options1;
    view_surf = &view_surf_1;
  }
  else {
    opt = &options2;
    view_surf = &view_surf_2;
  }
  opt->constrsw[constrno] = constraint[constrno];
  opt->constr_matrix_valid = false;
  CountTheConstraints ();
  opt->surf_valid = false;
  if ( *view_surf ) {
    if ( SurfFast ( opt ) ) {
      if ( constraints1 )
        UpdateFinalSurface1 ();
      else
        UpdateFinalSurface2 ();
    }
    else {
      *view_surf = false;
    }
    xge_Redraw ();
  }
} /*SwitchTheConstraint*/

