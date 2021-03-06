
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

#include "raybezprivatef.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define SOLUTION_ENTERED 1
#define SOLUTION_FINAL   2

static int EnterSolution ( float t0, float t1, float z,
                   boolean (*out)(void *usrptr, float t, boolean singular),
                   void *usrptr, boolean singular )
{
  if ( out ( usrptr, t0+z*(t1-t0), singular ) )
    return SOLUTION_ENTERED;
  else
    return SOLUTION_FINAL;
} /*EnterSolution*/

boolean rbez_FindBezcHighlightPointsf ( int degree, point3f *cp,
                     float t0, float t1,
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, float t, boolean singular),
                     void *usrptr )
{
  typedef struct {
      float t0, t1;
      float *ff;
      int   level;
    } stack_el;

  void     *sp, *sp1;
  vector3f *dt, *ca;
  float    *f, *ff, z;
  int      n, i, j;
  stack_el *stack;
  int      stp;

  sp = pkv_GetScratchMemTop ();
  n = degree+degree-1;
  f = pkv_GetScratchMemf ( n+1 );
  sp1 = pkv_GetScratchMemTop ();
  dt = pkv_GetScratchMem ( n+2*sizeof(vector3f) );
  if ( !f || !dt )
    goto failure;
  ca = &dt[degree];
        /* translate the curve */
  for ( i = 0; i <= degree; i++ )
    SubtractPoints3f ( &cp[i], a, &ca[i] );
        /* find the derivative of the curve */
  for ( i = 0; i < degree; i++ )
    SubtractPoints3f ( &cp[i+1], &cp[i], &dt[i] );
        /* construct the polynomial */
  mbs_multiBezScalef ( degree, 1, 1, 3, 0, &ca[0].x );
  mbs_multiBezScalef ( degree-1, 1, 1, 3, 0, &dt[0].x );
  memset ( f, 0, (n+1)*sizeof(float) );
  for ( i = 0; i <= degree; i++ )
    for ( j = 0; j < degree; j++ )
      f[i+j] += DotProduct3f ( &ca[i], &dt[j] );
  mbs_multiBezUnscalef ( n, 1, 1, 1, 0, f );
  pkv_SetScratchMemTop ( sp1 );
        /* prepare the stack */
  ff = pkv_GetScratchMemf ( maxlevel*(n+1) );
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
    if ( _rbez_ConvexHullTest1f ( n, ff ) ) {
      if ( _rbez_UniquenessTest1f ( n, ff ) ) {
        if ( _rbez_NewtonMethod1f ( n, ff, &z ) ) {
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
          mbs_BisectBC1f ( n, ff, stack[stp+1].ff );
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
} /*rbez_FindBezcHighlightPointsf*/

