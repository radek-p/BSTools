
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"

#include "xgedit.h"
#include "spl3d.h"
#include "ed3dswidgets.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean SetClosedUSurface ( void )
{
  int id;
  int i, k, pitch;

  if ( kwind.closed_u ) {
    sw_bind_blending = blending_mat_valid = false;
    if ( lastknot_u <= 3*degree_u ) {
      xge_DisplayErrorMessage ( ErrorMsgNotEnoughUKnots, -1 );
      kwind.closed_u = false;
      return false;
    }
    kwind.clcKu = lastknot_u-2*degree_u;
    kwind.clcTu = knots_u[kwind.clcKu+degree_u]-knots_u[degree_u];
    for ( i = 1; i < 2*degree_u; i++ )
      knots_u[i] = (double)(0.5*(knots_u[i]+knots_u[kwind.clcKu+i]-kwind.clcTu));
    for ( i = 1; i < 2*degree_u; i++ )
      knots_u[kwind.clcKu+i] = knots_u[i]+kwind.clcTu;
    knots_u[0] = knots_u[1];
    knots_u[lastknot_u] = knots_u[lastknot_u-1];
    pitch = lastknot_v-degree_v;
    k = kwind.clcKu*pitch;
    for ( i = 0; i < degree_u*pitch; i++ ) {
      MidPoint4d ( &cpoints[k+i], &cpoints[i], &cpoints[i] );
      cpoints[k+i] = cpoints[i];
    }
    ClearPointMarking (
        (lastknot_u-degree_u)*(lastknot_v-degree_v),
        mkpoints ); 
    for ( id = 0; id < 4; id++ )
      ProjectSurface ( id );
  }
  return true;
} /*SetClosedUSurface*/

boolean SetClosedVSurface ( void )
{
  int id;
  int i, j, k, pitch;

  if ( kwind.closed_v ) {
    sw_bind_blending = blending_mat_valid = false;
    if ( lastknot_v <= 3*degree_v ) {
      xge_DisplayErrorMessage ( ErrorMsgNotEnoughVKnots, -1 );
      kwind.closed_v = false;
      return false;
    }
    kwind.clcKv = lastknot_v-2*degree_v;
    kwind.clcTv = knots_v[kwind.clcKv+degree_v]-knots_v[degree_v];
    for ( i = 1; i < 2*degree_v; i++ )
      knots_v[i] = (double)(0.5*(knots_v[i]+knots_v[kwind.clcKv+i]-kwind.clcTv));
    for ( i = 1; i < 2*degree_v; i++ )
      knots_v[kwind.clcKv+i] = knots_v[i]+kwind.clcTv;
    knots_v[0] = knots_v[1];
    knots_v[lastknot_v] = knots_v[lastknot_v-1];
    pitch = lastknot_v-degree_v;
    for ( i = k = 0;  i < lastknot_u-degree_u;  i++, k += pitch )
      for ( j = 0; j < degree_v; j++ ) {
        MidPoint4d ( &cpoints[k+j+kwind.clcKv], &cpoints[k+j], &cpoints[k+j] );
        cpoints[k+j+kwind.clcKv] = cpoints[k+j];
      }
    ClearPointMarking (
        (lastknot_u-degree_u)*(lastknot_v-degree_v),
        mkpoints ); 
    for ( id = 0; id < 4; id++ )
      ProjectSurface ( id );
  }
  return true;
} /*SetClosedVSurface*/

