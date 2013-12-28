
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
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
#include "xgedit.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "bsfile.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean FilenameCorrect ( char *filename )
{
  return (boolean)(filename[0] != 0);
} /*FilenameCorrect*/

boolean SaveBSPatch ( char *filename )
{
  int lfn, lex;

  lfn = strlen ( filename );
  lex = strlen ( file_ext );
  if ( strcmp ( file_ext, &filename[lfn-lex] ) )
    strcpy ( &filename[lfn], file_ext );
  if ( bsf_OpenOutputFile ( filename, false ) ) {
    bsf_WriteBSplinePatchd ( 3, 4, eqmer_nurbs, degree_u, lastknot_u, knots_u,
                             degree_v, lastknot_v, knots_v,
                             kwind.closed_u, kwind.closed_v,
                             4*(lastknot_v-degree_v), &cpoints[0].x,
                             mkpoints, NULL, NULL, NULL );
    bsf_CloseOutputFile ();
    return true;
  }
  else
    return false;
} /*SaveBSPatch*/

boolean ReadBSPatch ( char *filename )
{
  int     pitch, spdimen;
  boolean result, closed_u, closed_v, rational;

  if ( bsf_OpenInputFile ( filename ) ) {
    result = bsf_ReadBSplinePatch4d ( MAX_DEGREE, MAX_KNOTS*(MAX_DEGREE+1)-1,
                                MAX_CPOINTS,
                                &degree_u, &lastknot_u, knots_u,
                                &degree_v, &lastknot_v, knots_v,
                                &closed_u, &closed_v,
                                &pitch, cpoints, &spdimen, &rational,
                                mkpoints, NULL );
    bsf_CloseInputFile ();
    if ( !result )
      bsf_PrintErrorLocation ();
    else {
      SetKWindNKnots ();
      kwind.closed_u = closed_u;
      if ( kwind.closed_u ) {
        kwind.clcKu = lastknot_u-2*degree_u;
        kwind.clcTu = knots_u[lastknot_u-degree_u]-knots_u[degree_u];
      }
      kwind.closed_v = closed_v;
      if ( kwind.closed_v ) {
        kwind.clcKv = lastknot_v-2*degree_v;
        kwind.clcTv = knots_v[lastknot_v-degree_v]-knots_v[degree_v];
      }
      if ( closed_u )
        blending_opt_part[0] = 0;
      else
        blending_opt_part[0] = degree_u;
      blending_opt_part[1] = lastknot_u-2*degree_u-1;
      blending_opt_part[2] = degree_v;
      blending_opt_part[3] = lastknot_v-2*degree_v-1;
      UpdateBlendingRangeWidgets ();
    }
    return result;
  }
  else
    return false;
} /*ReadBSPatch*/

