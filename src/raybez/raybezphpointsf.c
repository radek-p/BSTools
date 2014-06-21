
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define CONST_

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

#include "raybezprivatef.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define SOLUTION_OUTSIDE 0
#define SOLUTION_ENTERED 1
#define SOLUTION_FINAL   2

static int EnterSolution ( float u0, float u1, float v0, float v1,
                   vector2f *z,
                   boolean (*out)(void *usrptr, point2f *uv, boolean singular ),
                   void *usrptr, boolean singular )
{
  point2f uv;

  if ( z->x >= 0.0 && z->x <= 1.0 && z->y >= 0.0 && z->y <= 1.0 ) {
    uv.x = u0 + z->x*(u1-u0);
    uv.y = v0 + z->y*(v1-v0);
    if ( out ( usrptr, &uv, singular ) )
      return SOLUTION_ENTERED;
    else
      return SOLUTION_FINAL;
  }
  else
    return SOLUTION_OUTSIDE;
} /*EnterSolution*/

boolean rbez_FindBezpHighlightPointsf ( int n, int m, point3f *cp,
                     float u0, float u1, float v0, float v1,
                     point3f *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2f *uv, boolean singular),
                     void *usrptr )
{
  typedef struct {
      float    u0, u1, v0, v1;
      vector2f *mcp;
      int      level;
    } stack_el;

  void     *sp, *sp1;
  int      ncp, ncpp;
  int      nn, mm;
  int      i, j, k, l;
  int      stp;
  stack_el *stack;
  vector3f *pa, *du, *dv;
  vector2f *cpp, *mapcp, p, pu, pv, z;
  float    *acc, K1, K2;

  sp = pkv_GetScratchMemTop ();
        /* degree of algebraic equations */
  nn = n+n;  mm = m+m;
  ncpp = (nn+1)*(mm+1);
  cpp = pkv_GetScratchMem ( ncpp*sizeof(vector2f) );
  sp1 = pkv_GetScratchMemTop ();
  ncp = (n+1)*(m+1);
  pa = pkv_GetScratchMem ( ncp*sizeof(vector3f) );
  du = pkv_GetScratchMem ( n*(m+1)*sizeof(vector3f) );
  dv = pkv_GetScratchMem ( (n+1)*m*sizeof(vector3f) );
  acc = pkv_GetScratchMemf ( 2*ncpp );
  if ( !cpp || !pa || !du || !dv || !acc )
    goto failure;
        /* construct cpp - control points of algebraic equations to solve */
          /* translate the patch and find its (scaled) derivatives */
  for ( i = 0; i <= ncp; i++ )
    SubtractPoints3f ( &cp[i], a, &pa[i] );
  pkn_SubtractMatrixf ( 1, 3*n*(m+1), 0, &cp[m+1].x, 0, &cp[0].x, 0, &du[0].x );
  pkn_SubtractMatrixf ( n+1, 3*m, 3*(m+1), &cp[1].x, 3*(m+1), &cp[0].x, 3*m, &dv[0].x );
          /* transform to scaled Bernstein bases */
  mbs_multiBezScalef ( n, 1, 1, 3*(m+1), 0, &pa[0].x );
  mbs_multiBezScalef ( m, n+1, 1, 3, 0, &pa[0].x );
  mbs_multiBezScalef ( n-1, 1, 1, 3*(m+1), 0, &du[0].x );
  mbs_multiBezScalef ( m, n, 1, 3, 0, &du[0].x );
  mbs_multiBezScalef ( n, 1, 1, 3*m, 0, &dv[0].x );
  mbs_multiBezScalef ( m-1, n+1, 1, 3, 0, &dv[0].x );
          /* multiply the polynomials pa and du in scaled bases */
  memset ( acc, 0, ncpp*sizeof(float) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k < n; k++ )
        for ( l = 0; l <= m; l++ )
          acc[(i+k)*(mm+1)+j+l] += DotProduct3f ( &pa[i*(m+1)+j], &du[k*(m+1)+l] );
          /* copy with degree elevation with respect to "u" */
  memcpy ( &acc[ncpp], acc, ncpp*sizeof(float) );
  pkn_AddMatrixf ( 1, nn*(mm+1), 0, acc, 0, &acc[ncpp+mm+1], 0, &acc[ncpp+mm+1] );
  pkv_Selectf ( ncpp, 1, 1, 2, &acc[ncpp], &cpp[0].x );
          /* multiply the polynomials pa and dv in scaled bases */
  memset ( acc, 0, ncpp*sizeof(float) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k <= n; k++ )
        for ( l = 0; l < m; l++ )
          acc[(i+k)*mm+j+l] += DotProduct3f ( &pa[i*(m+1)+j], &dv[k*m+l] );
          /* copy with degree elevation with respect to "v" */
  memset ( &acc[ncpp], 0, ncpp*sizeof(float) );
  pkv_Selectf ( nn+1, mm, mm, mm+1, acc, &acc[ncpp] );
  pkn_AddMatrixf ( nn+1, mm, mm, acc, mm+1, &acc[ncpp+1], mm+1, &acc[ncpp+1] );
  pkv_Selectf ( ncpp, 1, 1, 2, &acc[ncpp], &cpp[0].y );
          /* transform to the unscaled basis */
  mbs_multiBezUnscalef ( nn, 1, 1, 2*(mm+1), 0, &cpp[0].x );
  mbs_multiBezUnscalef ( mm, nn+1, 1, 2, 0, &cpp[0].x );
        /* prepare the stack */
  pkv_SetScratchMemTop ( sp1 );
  mapcp = pkv_GetScratchMem ( maxlevel*ncpp*sizeof(vector2f) );
  stack = pkv_GetScratchMem ( (maxlevel+1)*sizeof(stack_el) );
  if ( !mapcp || !stack )
    goto failure;
  stack[0].u0 = u0;  stack[0].u1 = u1;
  stack[0].v0 = v0;  stack[0].v1 = v1;
  stack[0].mcp = cpp;
  stack[0].level = 0;
        /* recursive subdivision and solving the system of equations */
  stp = 1;
  do {
    stp --;
    u0 = stack[stp].u0;  u1 = stack[stp].u1;
    v0 = stack[stp].v0;  v1 = stack[stp].v1;
    mapcp = stack[stp].mcp;
    if ( _rbez_ConvexHullTest2f ( ncpp, mapcp ) ) {
      if ( _rbez_UniquenessTest2f ( nn, mm, ncpp, mapcp,
                                    &p, &pu, &pv, &K1, &K2 ) ) {
        switch ( _rbez_NewtonMethod2f ( nn, mm, mapcp, &p, &pu, &pv, &z ) ) {
      case RBEZ_NEWTON_YES:
          /* regular solution found */
          switch ( EnterSolution ( u0, u1, v0, v1, &z, out, usrptr, false ) ) {
        case SOLUTION_OUTSIDE:
            if ( _rbez_SecondTest2f ( &z, nn, mm, K1, K2 ) )
              goto subdivide;
            break;
        case SOLUTION_ENTERED:
            break;
        case SOLUTION_FINAL:
            goto finish;
          }
          break;
      case RBEZ_NEWTON_NO:
          goto subdivide;
      case RBEZ_NEWTON_ERROR:
          goto failure;
        }
      }
      else {
subdivide:
        if ( stack[stp].level < maxlevel ) {
          stack[stp+1].mcp = &stack[stp].mcp[ncpp];
          stack[stp+1].level = stack[stp].level = stack[stp].level+1;
          stack[stp+1].u0 = u0;
          stack[stp+1].v0 = v0;
          if ( _rbez_SubdividePatch2f ( nn, mm, stack[stp].mcp,
                                        stack[stp+1].mcp ) ) {
            stack[stp].u0 = stack[stp+1].u1 = 0.5*(u0+u1);
            stack[stp+1].v1 = v1;
          }
          else {
            stack[stp+1].u1 = u1;
            stack[stp].v0 = stack[stp+1].v1 = 0.5*(v0+v1);
          }
          stp += 2;
        }
        else {
          /* solution at singularity found */
          z.x = z.y = 0.5;
          switch ( EnterSolution ( u0, u1, v0, v1, &z, out, usrptr, true ) ) {
case SOLUTION_OUTSIDE:  /* cannot be */
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
} /*rbez_FindBezpHighlightPointsf*/
