
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <pthread.h>

#define CONST_

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivated.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define SOLUTION_ENTERED 1
#define SOLUTION_FINAL   2

static int EnterSolution ( double t0, double t1, double z,
                   boolean (*out)(void *usrptr, double t, boolean singular),
                   void *usrptr, boolean singular )
{
  if ( out ( usrptr, t0+z*(t1-t0), singular ) )
    return SOLUTION_ENTERED;
  else
    return SOLUTION_FINAL;
} /*EnterSolution*/

boolean rbez_FindRBezcHighlightPointsd ( int degree, point4d *cp,
                     double t0, double t1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, double t, boolean singular),
                     void *usrptr )
{
  typedef struct {
      double t0, t1;
      double *ff;
      int    level;
    } stack_el;

  void     *sp, *sp1;
  vector4d *scp, *dt;
  vector3d *rw, *ca, v;
  double   *f, *ff, z;
  int      n, n1, i, j;
  stack_el *stack;
  int      stp;

  sp = pkv_GetScratchMemTop ();
  n1 = 2*degree-1;
  n = n1+degree;
  f = pkv_GetScratchMemd ( n+1 );
  sp1 = pkv_GetScratchMemTop ();
  scp = pkv_GetScratchMem ( (2*degree+1)*sizeof(vector4d) );
  ca = pkv_GetScratchMem ( (3*degree+2)*sizeof(vector3d) );
  if ( !f || !scp || !ca )
    goto failure;
  dt = &scp[degree+1];
  rw = &ca[degree+1];
  memcpy ( scp, cp, (degree+1)*sizeof(vector4d) );
        /* translate the curve */
  for ( i = 0; i <= degree; i++ )
    AddVector3Md ( (vector3d*)&cp[i], a, -cp[i].w, &ca[i] );
        /* find the derivative of the curve */
  for ( i = 0; i < degree; i++ )
    SubtractPoints4d ( &cp[i+1], &cp[i], &dt[i] );
        /* construct the polynomial */
          /* find the second factor of the dot product in R^3 */
  mbs_multiBezScaled ( degree, 1, 1, 4, 0, &scp[0].x );
  mbs_multiBezScaled ( degree-1, 1, 1, 4, 0, &dt[0].x );
  memset ( rw, 0, (n1+1)*sizeof(vector3d) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j < degree; j++ ) {
      MultVector3d ( scp[i].w, (vector3d*)&dt[j], &v );
      AddVector3Md ( &v, (vector3d*)&scp[i], -dt[j].w, &v );
      AddVector3d ( &rw[i+j], &v, &rw[i+j] );
    }
          /* compute the dot product of the vector curves in R^3 */
  mbs_multiBezScaled ( degree, 1, 1, 3, 0, &ca[0].x );
  memset ( f, 0, (n+1)*sizeof(double) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j <= n1; j++ )
      f[i+j] += DotProduct3d ( &ca[i], &rw[j] );
  mbs_multiBezUnscaled ( n, 1, 1, 1, 0, f );
  pkv_SetScratchMemTop ( sp1 );
        /* prepare the stack */
  ff = pkv_GetScratchMemd ( maxlevel*(n+1) );
  stack = pkv_GetScratchMem ( (maxlevel+1)*sizeof(stack_el) );
  if ( !ff || !stack )
    goto failure;
  stack[0].t0 = t0;  stack[0].t1 = t1;
  stack[0].ff = f;
  stack[0].level = 0;
         /* recursive subdivision */
  stp = 1;
  do {
    stp --;
    t0 = stack[stp].t0;  t1 = stack[stp].t1;
    ff = stack[stp].ff;
    if ( _rbez_ConvexHullTest1d ( n, ff ) ) {
      if ( _rbez_UniquenessTest1d ( n, ff ) ) {
        if ( _rbez_NewtonMethod1d ( n, ff, &z ) ) {
          switch ( EnterSolution ( t0, t1, z, out, usrptr, false ) ) {
        case SOLUTION_ENTERED:
            break;
        case SOLUTION_FINAL:
            goto finish;
          }
        }
        else
          goto subdivide;
      }
      else {
subdivide:
        if ( stack[stp].level < maxlevel ) {
          stack[stp+1].ff = &ff[n+1];
          stack[stp+1].level = stack[stp].level = stack[stp].level+1;
          stack[stp+1].t0 = t0;
          stack[stp+1].t1 = stack[stp].t0 = 0.5*(t0+t1);
          mbs_BisectBC1d ( n, ff, stack[stp+1].ff );
          stp += 2;
        }
        else {
          switch ( EnterSolution ( t0, t1, 0.5, out, usrptr, true ) ) {
        case SOLUTION_ENTERED:
            break;
        case SOLUTION_FINAL:
            goto finish;
          }
        }
      }
    }
  } while ( stp > 0 );

finish:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*rbez_FindRBezcHighlightPointsd*/

