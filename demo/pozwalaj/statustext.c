
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

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

/* ///////////////////////////////////////////////////////////////////////// */
void SetStatusText ( char *text, boolean onscreen )
{
  if ( xge_CurrentWindow () == win0 ) {
    if ( text )
      strncpy ( status0, text, MAX_COMMAND_LGT );
    else
      status0[0] = 0;
    if ( onscreen ) {
      xge_SetClipping ( win0statl );
      win0statl->redraw ( win0statl, true );
    }
  }
  else if ( xge_CurrentWindow () == win1 ) {
    if ( text )
      strncpy ( status1, text, MAX_COMMAND_LGT );
    else
      status1[0] = 0;
    if ( onscreen ) {
      xge_SetClipping ( win1statl );
      win1statl->redraw ( win1statl, true );
    }
  }
} /*SetStatusText*/

void NotifyParam1 ( double param )
{
  char text[32];

  sprintf ( text, "%6.3f", param );
  SetStatusText ( text, true );
} /*NotifyParam1*/

void NotifyParam2 ( double param )
{
  char text[32];

  sprintf ( text, "%6.4f", param );
  SetStatusText ( text, true );
} /*NotifyParam2*/

/* ///////////////////////////////////////////////////////////////////////// */
static void DeleteChar ( char *s )
{
  while ( *s ) {
    *s = *(s+1);
    s ++;
  }
} /*DeleteChar*/

static char *FindChar ( char c, char *s )
{
  while ( *s && *s != c )
    s ++;
  return s;
} /*FindChar*/

static void WriteDoubleNumber ( double x, char *s )
{
#define EPS 1.0e-15
#define UX  1.0e+12
#define LX  1.0e-3
  double ax;
  char   *d, *e;
  int    l;

  ax = fabs ( x );
  if ( ax >= UX ) {
    sprintf ( s, "%15.12e", x );
    e = FindChar ( 'e', s ) +1;
    if ( *e == '+' ) e++;
    goto trunc_fract;
  }
  else if ( ax < EPS )
    memcpy ( s, "0.0\000", 4 );
  else if ( ax < LX ) {
    sprintf ( s, "%15.12e", x );
    e = FindChar ( 'e', s ) + 1;
    e = FindChar ( '-', e ) + 1;
trunc_fract:
    while ( *e == 0 && *(e+1) )
      DeleteChar ( e );
      e = FindChar ( 'e', s )-1;
      while ( *e == '0' )
        DeleteChar ( e-- );
  }
  else {
    sprintf ( s, "%14.12f", x );
    d = FindChar ( '.', s );
    if ( d-s < 14 )
      s[15] = 0;
    l = strlen ( s )-1;
    while ( s[l] == '0' && s[l-1] != '.' )
      s[l--] = 0;
  }
#undef LX
#undef UX
#undef EPS
} /*WriteDoubleNumber*/

void BottomDisplayPoint2D ( int win, int pn, point2d *p, boolean onscreen )
{
  char sx[20], sy[20], *s;

  if ( win == win0 )
    s = status0;
  else if ( win == win1 )
    s = status1;
  else
    return;
  if ( pn >= 0 ) {
    WriteDoubleNumber ( p->x, sx );
    WriteDoubleNumber ( p->y, sy );
    sprintf ( s, "%d: x = %s, y = %s", pn, sx, sy );
  }
  else
    s[0] = 0;
  if ( win == win0 ) {
    if ( win0statusline ) {
      xge_SetClipping ( win0statl );
      win0statl->redraw ( win0statl, onscreen );
    }
  }
  else if ( win == win1 ) {
    if ( win1statusline ) {
      xge_SetClipping ( win1statl );
      win1statl->redraw ( win1statl, onscreen );
    }
  }
} /*BottomDisplayPoint2D*/

void BottomDisplayPoint3D ( int win, int pn, point3d *p, boolean onscreen )
{
  char sx[20], sy[20], sz[20], *s;

  if ( win == win0 )
    s = status0;
  else if ( win == win1 )
    s = status1;
  else
    return;
  if ( pn >= 0 ) {
    WriteDoubleNumber ( p->x, sx );
    WriteDoubleNumber ( p->y, sy );
    WriteDoubleNumber ( p->z, sz );
    sprintf ( s, "%d: x = %s, y = %s, z = %s", pn, sx, sy, sz );
  }
  else
    s[0] = 0;
  if ( win == win0 ) {
    if ( win0statusline ) {
      xge_SetClipping ( win0statl );
      win0statl->redraw ( win0statl, onscreen );
    }
  }
  else if ( win == win1 ) {
    if ( win1statusline ) {
      xge_SetClipping ( win1statl );
      win1statl->redraw ( win1statl, onscreen );
    }
  }
} /*BottomDisplayPoint3D*/

void BottomDisplayPoint4D ( int win, int pn, point4d *p, boolean onscreen )
{
  char sx[20], sy[20], sz[20], sw[20], *s;

  if ( win == win0 )
    s = status0;
  else if ( win == win1 )
    s = status1;
  else
    return;
  if ( pn >= 0 ) {
    WriteDoubleNumber ( p->x, sx );
    WriteDoubleNumber ( p->y, sy );
    WriteDoubleNumber ( p->z, sz );
    WriteDoubleNumber ( p->w, sw );
    sprintf ( s, "%d: x = %s, y = %s, z = %s, sw = %s", pn, sx, sy, sz, sw );
  }
  else
    s[0] = 0;
  if ( win == win0 ) {
    if ( win0statusline ) {
      xge_SetClipping ( win0statl );
      win0statl->redraw ( win0statl, onscreen );
    }
  }
  else if ( win == win1 ) {
    if ( win1statusline ) {
      xge_SetClipping ( win1statl );
      win1statl->redraw ( win1statl, onscreen );
    }
  }
} /*BottomDisplayPoint4D*/

void BottomDisplayPoint ( int win, int spdimen, int cpdimen,
                          int pn, double *pc, boolean onscreen )
{
  switch ( cpdimen ) {
case 2:
    BottomDisplayPoint2D ( win, pn, (point2d*)pc, onscreen );
    break;
case 3:
    BottomDisplayPoint3D ( win, pn, (point3d*)pc, onscreen );
    break;
case 4:
    BottomDisplayPoint4D ( win, pn, (point4d*)pc, onscreen );
    break;
  }
} /*BottomDisplayPoint*/

