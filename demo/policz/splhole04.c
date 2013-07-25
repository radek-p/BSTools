
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>
#include <malloc.h>  
#include <string.h>     

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
#include "drawbezd.h"
#include "splhole.h"
#include "datagend.h"


point2s rlight_dir[4][R_NLIGHTS+1];

/* ////////////////////////////////////////////////////////////////////////// */
boolean FindNearestLightPoint ( int id, short x, short y, short mindist )
{
  short d, e;
  int   i, j;

  id &= 0x03;
  e = (short)(mindist+1);
  j = -1;
        /* first try the light directions */
  for ( i = 0; i < R_NLIGHTS; i++ )
    if ( edit_light_sw[i] ) {
      d = (short)(abs(rlight_dir[id][i].x-x)+abs(rlight_dir[id][i].y-y));
      if ( d < e ) {
        e = d;
        j = i;
      }
    }
        /* then the reflection lines frame */
        /* then the highlight lines frame */
        /* and the sections frame */

  swind.current_point = j;
  return e <= mindist;
} /*FindNearestLightPoint*/

void SetLightPoint ( int id, short x, short y )
{
  int     i;
  point3d p, q, r;

  id &= 0x03;
  swind_picture = false;
  i = swind.current_point;
  if ( i < R_NLIGHTS ) {
        /* set the light direction */
    p = swind.CPos[id].g_centre;
    SetPoint3d ( &q, p.x+light_dir[i].x, p.y+light_dir[i].y, p.z+light_dir[i].z );
    CameraProjectPoint3d ( &swind.CPos[id], &q, &r );
    r.x = (double)x;
    r.y = (double)y;
    CameraUnProjectPoint3d ( &swind.CPos[id], &r, &q );
    SubtractPoints3d ( &q, &p, &light_dir[i] );
    NormalizeVector3d ( &light_dir[i] );
  }
} /*SetLightPoint*/

void ProjectLightPoints ( int id )
{
  int     i;
  point3d p, q;

  id &= 0x03;
        /* project the light direction vectors */
  p = swind.CPos[id].g_centre;
  CameraProjectPoint3d2s ( &swind.CPos[id], &p, &rlight_dir[id][R_NLIGHTS] );
  for ( i = 0; i < R_NLIGHTS; i++ ) {
    SetPoint3d ( &q, p.x+light_dir[i].x, p.y+light_dir[i].y, p.z+light_dir[i].z );
    CameraProjectPoint3d2s ( &swind.CPos[id], &q, &rlight_dir[id][i] );
  }
        /* then the reflection lines frame */
        /* then the highlight lines frame */
        /* and the sections frame */
} /*ProjectLightPoints*/

void MyMarkPoint2s ( point2s *p )
{
  xgeFillRectangle ( 3, 3, (int)p->x-1, (int)p->y-1 );
} /*MyMarkPoint2s*/

void DrawLightPoints ( int id )
{
  int i;

  id &= 0x03;
  ProjectLightPoints ( id );
        /* draw the light direction vectors */
  for ( i = 0; i < R_NLIGHTS; i++ )
    if ( edit_light_sw[i] ) {
      xgeSetForeground ( xgec_Green );
      xgeDrawLine ( rlight_dir[id][R_NLIGHTS].x, rlight_dir[id][R_NLIGHTS].y,
                    rlight_dir[id][i].x, rlight_dir[id][i].y );
      xgeSetForeground ( xgec_Cyan );
      MyMarkPoint2s ( &rlight_dir[id][R_NLIGHTS] );
      xgeSetForeground ( xgec_Yellow );
      MyMarkPoint2s ( &rlight_dir[id][i] );
    }
        /* then the reflection lines frame */
        /* then the highlight lines frame */
        /* and the sections frame */
} /*DrawLightPoints*/

