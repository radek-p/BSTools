
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

boolean WritePatchAttributes ( void *usrdata )
{
  void         *sp;
  int          ncp, i;
  unsigned int *mk;

  sp = pkv_GetScratchMemTop ();
  ncp = (lastknot_u-degree_u)*(lastknot_v-degree_v);
  mk = (unsigned int*)pkv_GetScratchMemi ( ncp );
  if ( mk ) {
    for ( i = 0; i < ncp; i++ )
      mk[i] = (unsigned int)mkpoints[i];
    bsf_WritePointsMK ( ncp, mk );
  }
  pkv_SetScratchMemTop ( sp );
  return true;
} /*WritePatchAttributes*/

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
                             NULL, -1, WritePatchAttributes, NULL );
    bsf_CloseOutputFile ();
    return true;
  }
  else
    return false;
} /*SaveBSPatch*/

void BSPatchReader ( void *userData, const char *name, int ident,
                     int udeg, int lastknotu, const double *knotsu,
                     int vdeg, int lastknotv, const double *knotsv,
                     boolean closed_u, boolean closed_v,
                     int pitch, const point4d *_cpoints,
                     int spdimen, boolean rational )
{
  int ncp;

  degree_u = udeg;
  lastknot_u = lastknotu;
  memcpy ( knots_u, knotsu, (lastknotu+1)*sizeof(double) );
  degree_v = vdeg;
  lastknot_v = lastknotv;
  memcpy ( knots_v, knotsv, (lastknotv+1)*sizeof(double) );
  ncp = (lastknotu-udeg)*(lastknotv-vdeg);
  memcpy ( cpoints, _cpoints, ncp*sizeof(point4d) );
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
  sw_bind_blending = false;
  blending_mat_valid = false;
  n_blending_constraints = 0;
} /*BSPatchReader*/

void CPMarkReader ( void *userData, int ncp, unsigned int *mk )
{
  int i;

  if ( ncp > MAX_CPOINTS )
    ncp = MAX_CPOINTS;
  for ( i = 0; i < ncp; i++ )
    mkpoints[i] = (byte)mk[i];
} /*CPMarkReader*/

boolean ReadBSPatch ( char *filename )
{
  bsf_UserReaders readers;

  bsf_ClearReaders ( &readers );
  bsf_BSP4ReadFuncd ( &readers, BSPatchReader, MAX_DEGREE,
                      MAX_KNOTS*(MAX_DEGREE+1)-1 );
  bsf_CPMarkReadFunc ( &readers, CPMarkReader );
  readers.userData = &readers;
  if ( !bsf_ReadBSFiled ( filename, &readers ) ) {
    bsf_PrintErrorLocation ();
    ResetObject ();
    return false;
  }
  return true;
} /*ReadBSPatch*/

