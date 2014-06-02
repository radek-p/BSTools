
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

#include "spl3d.h"

void ValidateG2ConstrKnots ( void )
{
  int i;

  for ( i = 1; i <= n_blending_constraints; i++ )
    blending_constr_poly_valid[i] =
      blending_constr_knots[i] > knots_u[blending_opt_part[0]] &&
      blending_constr_knots[i] < knots_u[blending_opt_part[1]+degree_u+1];
} /*ValidateG2ConstrKnots */

boolean InsertG2BlendingConstrKnot ( double newknot )
{
  int r, id;

  if ( n_blending_constraints < MAX_BLENDING_CONSTRAINTS &&
       newknot > knots_u[degree_u] &&
       newknot < knots_u[lastknot_u-degree_u] ) {
    if ( n_blending_constraints ) {
      r = mbs_KnotMultiplicityd ( n_blending_constraints,
                                  &blending_constr_knots[1], newknot );
      if ( r )
        return false;
      for ( r = n_blending_constraints;
            r > 0 && blending_constr_knots[r] > newknot;
            r -- ) {
        blending_constr_knots[r+1] = blending_constr_knots[r];
        memcpy ( &blending_constr_cp[r*(lastknot_v-degree_v)],
                 &blending_constr_cp[(r-1)*(lastknot_v-degree_v)],
                 (lastknot_v-degree_v)*sizeof(point4d) );
      }
      r ++;
      blending_constr_knots[r] = newknot;
      n_blending_constraints++;
      kwind.lastknot_u = n_blending_constraints + 1;
    }
    else {
      blending_constr_knots[1] = newknot;
      n_blending_constraints = r = 1;
      kwind.lastknot_u = 2;
    }
    FindG2BlendingConstrCP ( r );
    ValidateG2ConstrKnots ();
    for ( id = 0; id < 4; id ++ )
      ProjectG2BlendingConstrCPoly ( id );
    return true;
  }
  return false;
} /*InsertG2BlendingConstrKnot*/

boolean RemoveG2BlendingConstrKnot ( int knotnum )
{
  int i, j, id;

  if ( knotnum > 0 && knotnum <= n_blending_constraints ) {
    j = lastknot_v-degree_v;
    n_blending_constraints --;
    kwind.lastknot_u = n_blending_constraints + 1;
    for ( i = knotnum; i <= n_blending_constraints; i++ ) {
      blending_constr_knots[i] = blending_constr_knots[i+1];
      memcpy ( &blending_constr_cp[j*(i-1)], &blending_constr_cp[j*i],
               j*sizeof(point4d) );
    }
    ValidateG2ConstrKnots ();
    for ( id = 0; id < 4; id ++ )
      ProjectG2BlendingConstrCPoly ( id );
    return true;
  }
  return false;
} /*RemoveG2BlendingConstrKnot*/

void ChangeG2BlendingConstrKnot ( int oldnum, int knotnum )
{
  int id;

  if ( knotnum > oldnum ) {
    do {
      memcpy ( &blending_constr_cp[(oldnum-1)*(lastknot_v-degree_v)],
               &blending_constr_cp[oldnum*(lastknot_v-degree_v)],
               (lastknot_v-degree_v)*sizeof(point4d) );
      oldnum ++;
    } while ( knotnum > oldnum );
  }
  else if ( knotnum < oldnum ) {
    do {
      oldnum --;
      memcpy ( &blending_constr_cp[oldnum*(lastknot_v-degree_v)],
               &blending_constr_cp[(oldnum-1)*(lastknot_v-degree_v)],
               (lastknot_v-degree_v)*sizeof(point4d) );
    } while ( knotnum < oldnum );
  }
  FindG2BlendingConstrCP ( knotnum );
  ValidateG2ConstrKnots ();
  for ( id = 0; id < 4; id ++ )
    ProjectG2BlendingConstrCPoly ( id );
} /*ChangeG2BlendingConstrKnot*/

void ProjectG2BlendingConstrCP ( int id, int knotnum )
{
  int     i, j;

  id &= 0x03;
  if ( blending_constr_poly_valid[knotnum] ) {
    j = (knotnum-1)*(lastknot_v-degree_v );
    for ( i = 0; i < lastknot_v-degree_v; i++ )
      clblending_constr[id][j+i] = RzutujPunkt ( id,
                &blending_constr_cp[j+i], &blending_constr_rp[id][j+i] );
  }
} /*ProjectG2BlendingConstrCP*/

void ProjectG2BlendingConstrCPoly ( int id )
{
  int k;

  id &= 0x03;
  for ( k = 1; k <= n_blending_constraints; k++ )
    ProjectG2BlendingConstrCP ( id, k );
} /*ProjectG2BlendingConstrCPoly*/

void FindG2BlendingConstrCP ( int knotnum )
{
  if ( knotnum < 1 || knotnum > n_blending_constraints )
    return;
  blending_constr_poly_valid[knotnum] = false;
  if ( blending_constr_knots[knotnum] <= 3.0 ||
       blending_constr_knots[knotnum] >= lastknot_u-3.0 )
    return;
  mbs_multideBoord ( degree_u, lastknot_u, knots_u, 1,
                     4*(lastknot_v-degree_v), 0, (double*)cpoints,
                     blending_constr_knots[knotnum],
                     &blending_constr_cp[(knotnum-1)*(lastknot_v-degree_v)].x );
  blending_constr_poly_valid[knotnum] = true;
} /*FindG2BlendingConstrCP*/

void DisplayG2BlendingConstrCP ( int id )
{
  void   *sp;
  int    i, j, k;
  double *ccb;

  id &= 0x03;
  sp = pkv_GetScratchMemTop ();
  k = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  ccb = pkv_GetScratchMemd ( k*(degree_v+1)*4 );
  if ( ccb ) {
    xgeSetForeground ( xgec_Plum );
    for ( i = 1; i <= n_blending_constraints; i++ )
      if ( blending_constr_poly_valid[i] ) {
        mbs_BSToBezC4d ( degree_v, lastknot_v, knots_v,
                         &blending_constr_cp[(i-1)*(lastknot_v-degree_v)],
                         NULL, NULL, NULL, ccb );
        for ( j = 0; j < k; j++ )
          DisplayBezCurve4 ( &swind.CPos[id], degree_v,
                             (point4d*)&ccb[4*j*(degree_v+1)] );
      }
  }
  pkv_SetScratchMemTop ( sp );

  xgeSetForeground ( xgec_Green );
  for ( i = 1; i <= n_blending_constraints; i++ )
    if ( blending_constr_poly_valid[i] ) {
      k = (i-1)*(lastknot_v-degree_v);
      for ( j = 0; j < lastknot_v-degree_v-1; j++ )
        if ( clblending_constr[id][k+j] &&
             clblending_constr[id][k+j+1] )
          xgeDrawLine ( (int)(blending_constr_rp[id][k+j].x+0.5),
                        (int)(blending_constr_rp[id][k+j].y+0.5),
                        (int)(blending_constr_rp[id][k+j+1].x+0.5),
                        (int)(blending_constr_rp[id][k+j+1].y+0.5) );
        else
          DisplayLine4 ( id, &blending_constr_cp[k+j],
                         &blending_constr_cp[k+j+1] );
    }
  for ( i = 1; i <= n_blending_constraints; i++ )
    if ( blending_constr_poly_valid[i] ) {
      k = (i-1)*(lastknot_v-degree_v);
      for ( j = 0; j < lastknot_v-degree_v; j++ ) {
        if ( j < degree_v || j >= lastknot_v-2*degree_v )
          xgeSetForeground ( xgec_Green );
        else {
          if ( mkblending_cp[k+j] )
            xgeSetForeground ( xgec_OrangeRed );
          else
            xgeSetForeground ( xgec_Yellow );
        }
        if ( clblending_constr[id][k+j] )
          xgeFillRectangle ( 3, 3, (int)(blending_constr_rp[id][k+j].x+0.5),
                                   (int)(blending_constr_rp[id][k+j].y+0.5) );
      }
    }
} /*DisplayG2BlendingConstrCP*/

