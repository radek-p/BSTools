
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"

#include "xgledit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
static void xgle_3DwinfDrawXYZ ( CameraRecf *CPos, char cc, char left,
                                 int xi, int eta, float xyz )
{
  char s[30];

  glColor3fv ( xglec_Green );
  switch ( cc ) {
case 0: sprintf ( s, "x=%5.3f", xyz );  break;
case 1: sprintf ( s, "y=%5.3f", xyz );  break;
case 2: sprintf ( s, "z=%5.3f", xyz );  break;
  }
  if ( left ) {
    xi = CPos->xmin+2;
    eta = max ( eta, CPos->ymin+11 );
  }
  else {
    xi = min ( xi+2, CPos->xmin+CPos->width-6*strlen(s) );
    eta = CPos->ymin+CPos->height-4;
  }
  xgleDrawString ( s, xi, eta );
} /*xgle_3DwinfDrawXYZ*/

void xgle_3DwinfDrawCursorPos ( xge_3Dwinf *_3Dwin,
                                int id, short x, short y )
{
  point3f p, q;
  char    zoomwin;

  id &= 0x03;
  if ( !_3Dwin->display_coord ||
       _3Dwin->CoordWin < 0 || _3Dwin->CoordWin > 2 ||
       id < 0 || id > 2 )
    return;
  xgle_SetIdentMapping ( _3Dwin->cwin[id] );
  zoomwin = _3Dwin->fww.zoomwin;
  glColor3fv ( xglec_Green5 );
  switch ( _3Dwin->CoordWin ) {
case 0:
    switch ( id ) {
  case 0:
      xgleDrawLine ( x, _3Dwin->CPos[id].ymin,
                     x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgleDrawLine ( _3Dwin->CPos[id].xmin, y,
                     _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      if ( zoomwin == 0 )
        xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
  case 1:   
      xgleDrawLine ( _3Dwin->CPos[id].xmin, y,
                     _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y );
      break;
  case 2:   
      xgleDrawLine ( x, _3Dwin->CPos[id].ymin,
                     x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
    }
    break;
case 1:   
    switch ( id ) {
  case 0:
      xgleDrawLine ( _3Dwin->CPos[id].xmin, y,
                    _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y );
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      break;
  case 1:   
      xgleDrawLine ( x, _3Dwin->CPos[id].ymin,
                     x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgleDrawLine ( _3Dwin->CPos[id].xmin, y,
                     _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 1, 0, x, y, p.y );
      if ( zoomwin == 1 )
        xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 2, 1, x, y, p.z );
      break;
  case 2:   
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[_3Dwin->CoordWin], &q, &p );
      CameraProjectPoint3f ( &_3Dwin->CPos[id], &p, &q );   
      xgleDrawLine ( _3Dwin->CPos[id].xmin, (int)q.y,
                     _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, (int)q.y );
      break;
    }
    break;
case 2:   
    switch ( id ) {
  case 0:
      xgleDrawLine ( x, _3Dwin->CPos[id].ymin,
                     x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      break;
  case 1:   
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[_3Dwin->CoordWin], &q, &p );
      CameraProjectPoint3f ( &_3Dwin->CPos[id], &p, &q );   
      xgleDrawLine ( (int)q.x, _3Dwin->CPos[id].ymin,
                     (int)q.x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      break;
  case 2:   
      xgleDrawLine ( x, _3Dwin->CPos[id].ymin,
                     x, _3Dwin->CPos[id].ymin+_3Dwin->CPos[id].height );
      xgleDrawLine ( _3Dwin->CPos[id].xmin, y,
                     _3Dwin->CPos[id].xmin+_3Dwin->CPos[id].width, y ); 
      SetPoint3f ( &q, (float)x, (float)y, 0.0 );
      CameraUnProjectPoint3f ( &_3Dwin->CPos[id], &q, &p );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 1, 1, x, y, p.y );
      xgle_3DwinfDrawXYZ ( &_3Dwin->CPos[id], 0, 0, x, y, p.x );
      break;
    }
    break;
  }
} /*xgle_3DwinfDrawCursorPos*/

