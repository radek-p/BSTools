
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"
#include "convh.h"
#include "camerad.h"
#include "xgedit.h"
#include "xgeipc.h"

#include "render.h"
#include "spl3d.h"
#include "ed3ds.h"
#include "pomnijipc.h"

/* ////////////////////////////////////////////////////////////////////////// */
static void _PatchDataSize ( int *s0, int *s1, int *s2, int *s3 )
{
  *s0 = 5*sizeof(int);
  *s1 = (lastknot_u+1)*sizeof(double);
  *s2 = (lastknot_v+1)*sizeof(double);
  *s3 = ((lastknot_u-degree_u)*
         (lastknot_v-degree_v)*4)*sizeof(double);
} /*_PatchDataSize*/

int PatchDataSize ( void )
{
  int s0, s1, s2, s3;

  _PatchDataSize ( &s0, &s1, &s2, &s3 );
  return s0 + s1 + s2 + s3;
} /*PatchDataSize*/

void SendPatchToChild ( void )
{
  int s0, s1, s2, s3;
  int intpar[5];

  _PatchDataSize ( &s0, &s1, &s2, &s3 );
  intpar[0] = degree_u;
  intpar[1] = degree_v;
  intpar[2] = lastknot_u;
  intpar[3] = lastknot_v;
  intpar[4] = 4;
  write ( xge_pipe_in[1], intpar, s0 );
  write ( xge_pipe_in[1], knots_u, s1 );
  write ( xge_pipe_in[1], knots_v, s2 );
  write ( xge_pipe_in[1], cpoints, s3 );
} /*SendPatchToChild*/

int ConstraintDataSize ( void )
{
  return sizeof(int) + n_blending_constraints*sizeof(double) +
         n_blending_constraints*(lastknot_v-3)*sizeof(point4d);
} /*ConstraintDataSize*/

void SendConstraintsToChild ( void )
{
  write ( xge_pipe_in[1], &n_blending_constraints, sizeof(int) );
  write ( xge_pipe_in[1], &blending_constr_knots[1],
          n_blending_constraints*sizeof(double) );
  write ( xge_pipe_in[1], blending_constr_cp,
          n_blending_constraints*(lastknot_v-degree_v)*sizeof(point4d) );
} /*SendConstraintsToChild*/

void SendOptionsToChild ( void )
{
  ipc_options options;

  options.NLConst = xge_LogSlidebarValued ( NLBLENDING_CMIN, NLBLENDING_CMAX,
                                            blending_factor );
  options.maxiter = blending_lmt_iter;
  options.nkn1 = blending_quad1;
  options.nkn2 = blending_quad2;
  memcpy ( &options.opt_range[0], &blending_opt_part[0], 4*sizeof(int) );
  if ( sw_blending_constraints )
    options.nconstr = n_blending_constraints;
  else
    options.nconstr = 0;
  if ( sw_blending_g2 )
    options.gcont = 2;
  else
    options.gcont = 1;
  options.send_partial = sw_show_steps;
  options.closed = kwind.closed_u;
  options.dumpdata = sw_blending_opt_dump;
  options.trans = blending_opt_transform;
  write ( xge_pipe_in[1], &options, sizeof(options) );
} /*SendOptionsToChild*/

void GetPatchFromChild ( void )
{
  int i, ncp;
  int s1, s2, s3, dim;
  int intpar[5];

  read ( xge_pipe_out[0], intpar, 5*sizeof(int) );
  s1 = (intpar[2]+1)*sizeof(double);
  s2 = (intpar[3]+1)*sizeof(double);
  ncp = (intpar[2]-intpar[0])*(intpar[3]-intpar[1]);
  dim = intpar[4];
  s3 = ncp*dim*sizeof(double);
  degree_u = intpar[0];
  degree_v = intpar[1];
  lastknot_u = intpar[2];
  lastknot_v = intpar[3];
  read ( xge_pipe_out[0], knots_u, s1 );
  read ( xge_pipe_out[0], knots_v, s2 );
  read ( xge_pipe_out[0], cpoints, s3 );
  if ( dim == 3 ) {
    pkv_Rearranged ( ncp, 3, 3, 4, cpoints );
    for ( i = 0; i < ncp; i++ ) 
      cpoints[i].w = 1.0;
  }
  switch ( degree_u ) {
case 2:
    if ( knots_u[1] == knots_u[2] &&
         knots_v[1] == knots_v[2] ) {
      G1FreeBoundaryToClamped ();
      sw_clamped_blending = true;
    }
    else {
      G1ClampedBoundaryToFree ();
      sw_clamped_blending = false;
    }
    break;
case 3:
    if ( knots_u[1] == knots_u[3] &&
         knots_v[1] == knots_v[3] ) {
      G2FreeBoundaryToClamped ();
      sw_clamped_blending = true;
    }
    else {
      G2ClampedBoundaryToFree ();
      sw_clamped_blending = false;
    }
    break;
default:
    break;
  }
} /*GetPatchFromChild*/

void GetOptionsFromChild ( void )
{
} /*GetOptionsFromChild*/

void GetConstraintsFromChild ( void )
{
} /*GetConstraintsFromChild*/

