
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

/* ///////////////////////////////////////////////////////////////////////// */
void SaveControlPoints ( void )
{
  int ncp;

  ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
  memcpy ( savedcpoints, cpoints, ncp*sizeof(point4d) );
  if ( n_blending_constraints ) {
    ncp = (lastknot_v-degree_v)*n_blending_constraints;
    memcpy ( savedbl_constr_cp, blending_constr_cp, ncp*sizeof(point4d) );
  }
} /*SaveControlPoints*/

void TransformMarkedControlPoints ( trans3d *tr )
{
  int ncp, i, id;

  ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
  for ( i = 0; i < ncp; i++ )
    if ( mkpoints[i] & 0x01 )
      Trans3Point4d ( tr, &savedcpoints[i], &cpoints[i] );
  if ( n_blending_constraints ) {
    ncp = (lastknot_v-degree_v)*n_blending_constraints;
    for ( i = 0; i < ncp; i++ )
      if ( mkblending_cp[i] & 0x01 )
        Trans3Point4d ( tr, &savedbl_constr_cp[i], &blending_constr_cp[i] );
  }
  if ( sw_bind_blending ) {
    if ( sw_blending_g1 )
      sw_bind_blending = ConstructG1BlendingSurface ();
    else
      sw_bind_blending = ConstructG2BlendingSurface ();
  }
  for ( id = 0; id < 4; id++ ) {
    ProjectSurface ( id );
    ProjectG2BlendingConstrCPoly ( id );
  }
} /*TransformMarkedControlPoints*/

