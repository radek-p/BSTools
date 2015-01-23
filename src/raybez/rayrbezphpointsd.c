
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014, 2015                            */
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

#include "raybezprivated.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define SOLUTION_OUTSIDE 0
#define SOLUTION_ENTERED 1
#define SOLUTION_FINAL   2

static int EnterSolution ( double u0, double u1, double v0, double v1,
                   vector2d *z,
                   boolean (*out)(void *usrptr, point2d *uv, boolean singular ),
                   void *usrptr, boolean singular )
{
  point2d uv;

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

boolean rbez_FindRBezpHighlightPointsd ( int n, int m, point4d *cp,
                     double u0, double u1, double v0, double v1,
                     point3d *a, int maxlevel,
                     boolean (*out)(void *usrptr, point2d *uv, boolean singular),
                     void *usrptr )
{
  typedef struct {
      double   u0, u1, v0, v1;
      vector2d *mcp;
      int      level;
    } stack_el;

  void     *sp;
  int      ncp, ncpp, size;
  int      nn, mm, nn1, mm1;
  int      i, j, k, l;
  int      stp;
  stack_el *stack;
  vector4d *scp, *du, *dv;
  vector3d *qw, *pa, v;
  vector2d *cpp, *mapcp, p, pu, pv, z;
  double   *acc, K1, K2;
  double   *workspace;

  ncp = (n+1)*(m+1);
  nn1 = n+n-1;  mm1 = m+m-1;
        /* degree of algebraic equations */
  nn = 3*n;    mm = 3*m;
  ncpp = (nn+1)*(mm+1);
  size = ncpp*sizeof(vector2d) +
         max ( ncp*sizeof(vector4d) + n*(m+1)*sizeof(vector4d) +
               (n+1)*m*sizeof(vector4d) + ncp*sizeof(vector3d) +
               2*ncpp*sizeof(double) +
               (nn1+1)*(mm1+1)*sizeof(vector3d),
               (maxlevel+1)*(ncpp*sizeof(vector2d)+
                             sizeof(stack_el)) +   
               (6+2*(mm+1))*2*sizeof(double) );
  sp = cpp = pkv_GetScratchMem ( size );
  if ( !cpp )
    goto failure;
  scp = (vector4d*)&cpp[ncpp];
  du = &scp[ncp];
  dv = &du[n*(m+1)];
  pa = (vector3d*)&dv[(n+1)*m];
  acc = (double*)&pa[ncp];
  qw = (vector3d*)&acc[2*ncpp];
  memcpy ( scp, cp, ncp*sizeof(vector4d) );
        /* construct cpp - control points of algebraic equations to solve */
          /* translate the patch and find its (scaled) derivatives */
  for ( i = 0; i <= ncp; i++ )
    AddVector3Md ( (vector3d*)&cp[i], a, -cp[i].w, &pa[i] );
  pkn_SubtractMatrixd ( 1, 4*n*(m+1), 0, &cp[m+1].x, 0, &cp[0].x, 0, &du[0].x );
  pkn_SubtractMatrixd ( n+1, 4*m, 4*(m+1), &cp[1].x, 4*(m+1), &cp[0].x, 4*m, &dv[0].x );
          /* transform to scaled Bernstein bases */
  mbs_multiBezScaled ( n, 1, 1, 4*(m+1), 0, &scp[0].x );
  mbs_multiBezScaled ( m, n+1, 1, 4, 0, &scp[0].x );
  mbs_multiBezScaled ( n, 1, 1, 3*(m+1), 0, &pa[0].x );
  mbs_multiBezScaled ( m, n+1, 1, 3, 0, &pa[0].x );
  mbs_multiBezScaled ( n-1, 1, 1, 4*(m+1), 0, &du[0].x );
  mbs_multiBezScaled ( m, n, 1, 4, 0, &du[0].x );
  mbs_multiBezScaled ( n, 1, 1, 4*m, 0, &dv[0].x );
  mbs_multiBezScaled ( m-1, n+1, 1, 4, 0, &dv[0].x );
          /* compute the second factor of the first dot product in R^3 */
  memset ( qw, 0, (nn1+1)*(mm1+2)*sizeof(vector3d) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k < n; k++ )
        for ( l = 0; l <= m; l++ ) {
          MultVector3d ( scp[i*(m+1)+j].w, (vector3d*)&du[k*(m+1)+l], &v );
          AddVector3Md ( &v, (vector3d*)&scp[i*(m+1)+j], -du[k*(m+1)+l].w, &v );
          AddVector3d ( &qw[(i+k)*(m+m+1)+j+l], &v, &qw[(i+k)*(m+m+1)+j+l] );
        }
          /* compute the first dot product of two patches in R^3 */
  memset ( acc, 0, ncpp*sizeof(double) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k < n+n; k++ )
        for ( l = 0; l <= m+m; l++ )
          acc[(i+k)*(mm+1)+j+l] += DotProduct3d ( &pa[i*(m+1)+j], &qw[k*(m+m+1)+l] );
          /* copy with degree elevation with respect to "u" */
  memcpy ( &acc[ncpp], acc, ncpp*sizeof(double) );
  pkn_AddMatrixd ( 1, nn*(mm+1), 0, acc, 0, &acc[ncpp+mm+1], 0, &acc[ncpp+mm+1] );
  pkv_Selectd ( ncpp, 1, 1, 2, &acc[ncpp], &cpp[0].x );
          /* compute the second factor of the second dot product in R^3 */
  memset ( qw, 0, (nn1+2)*(mm1+1)*sizeof(vector3d) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k <= n; k++ )
        for ( l = 0; l < m; l++ ) {
          MultVector3d ( scp[i*(m+1)+j].w, (vector3d*)&dv[k*m+l], &v );
          AddVector3Md ( &v, (vector3d*)&scp[i*(m+1)+j], -dv[k*m+l].w, &v );
          AddVector3d ( &qw[(i+k)*(m+m)+j+l], &v, &qw[(i+k)*(m+m)+j+l] );
        }
          /* multiply the polynomials pa and dv in scaled bases */
  memset ( acc, 0, ncpp*sizeof(double) );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j <= m; j++ )
      for ( k = 0; k <= n+n; k++ )
        for ( l = 0; l < m+m; l++ )
          acc[(i+k)*mm+j+l] += DotProduct3d ( &pa[i*(m+1)+j], &qw[k*(m+m)+l] );
          /* copy with degree elevation with respect to "v" */
  memset ( &acc[ncpp], 0, ncpp*sizeof(double) );
  pkv_Selectd ( nn+1, mm, mm, mm+1, acc, &acc[ncpp] );
  pkn_AddMatrixd ( nn+1, mm, mm, acc, mm+1, &acc[ncpp+1], mm+1, &acc[ncpp+1] );
  pkv_Selectd ( ncpp, 1, 1, 2, &acc[ncpp], &cpp[0].y );
          /* transform to the unscaled basis */
  mbs_multiBezUnscaled ( nn, 1, 1, 2*(mm+1), 0, &cpp[0].x );
  mbs_multiBezUnscaled ( mm, nn+1, 1, 2, 0, &cpp[0].x );
        /* prepare the stack */
  mapcp = &cpp[ncpp];
  stack = (stack_el*)&mapcp[maxlevel+1];
  workspace = (double*)&stack[maxlevel+1];
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
    if ( _rbez_ConvexHullTest2d ( ncpp, mapcp ) ) {
      if ( _rbez_UniquenessTest2d ( nn, mm, ncpp, mapcp,
                                    &p, &pu, &pv, &K1, &K2, workspace ) ) {
        switch ( _rbez_NewtonMethod2d ( nn, mm, mapcp,
                                        &p, &pu, &pv, &z, workspace ) ) {
      case RBEZ_NEWTON_YES:
          /* regular solution found */
          switch ( EnterSolution ( u0, u1, v0, v1, &z, out, usrptr, false ) ) {
        case SOLUTION_OUTSIDE:
            if ( _rbez_SecondTest2d ( &z, nn, mm, K1, K2 ) )
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
          if ( _rbez_SubdividePatch2d ( nn, mm, stack[stp].mcp,
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
} /*rbez_FindRBezpHighlightPointsd*/

