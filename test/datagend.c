
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "eg2holed.h"
#include "datagend.h"


/* built-in data generator */

static int hole_k;           /* number of hole sides, 3, 5, 6 or 8 */
int hole_np;                 /* number of B-spline control points */
point3d Bpt[MAX_BPTS];       /* B-spline control points */
point2d Dompt[MAX_BPTS];     /* B-spline domain control points */
double knots[11*GH_MAX_K];

/* ////////////////////////////////////////////////////////////////////////// */
void InitKnots ( int kk )
{
  int i;

  knots[0] = knots[1] = 0.0;
  for ( i = 2; i < 9; i++ )
    knots[i] = (double)(i-1)/8.0;
  knots[9] = knots[10] = 1.0;
  for ( i = 1; i < kk; i++ )
    memcpy ( &knots[11*i], &knots[0], 11*sizeof(double) );
} /*InitKnots*/

/* ////////////////////////////////////////////////////////////////////////// */
point3d hdp3[37] =
 {{0.0, 0.0, 0.0},
  {0.0, -1.0,0.0},{0.0, -2.0,0.0},{0.0, -3.0,0.0},
  {1.0, -1.0,0.0},{1.0, -2.0,0.0},{1.0, -3.0,0.0},
  {2.0, -1.0,0.0},{2.0, -2.0,0.0},{2.0, -3.0,0.0},
  {3.0, -1.0,0.0},{3.0, -2.0,0.0},{3.0, -3.0,0.0},
  {1.0, 0.0, 0.0},{2.0, 0.0, 0.0},{3.0,  0.0,0.0},
  {1.0, 0.0, 1.0},{2.0, 0.0, 1.0},{3.0,  0.0,1.0},
  {1.0, 0.0, 2.0},{2.0, 0.0, 2.0},{3.0,  0.0,2.0},
  {1.0, 0.0, 3.0},{2.0, 0.0, 3.0},{3.0,  0.0,3.0},
  {0.0, 0.0, 1.0},{0.0, 0.0, 2.0},{0.0,  0.0,3.0},
  {0.0, -1.0,1.0},{0.0, -1.0,2.0},{0.0, -1.0,3.0},
  {0.0, -2.0,1.0},{0.0, -2.0,2.0},{0.0, -2.0,3.0},
  {0.0, -3.0,1.0},{0.0, -3.0,2.0},{0.0, -3.0,3.0}};

point3d hdp5[61] =
 {{0.0, 0.0, 0.0},
  {0.0, -1.0,0.0},{0.0, -2.0,0.0},{0.0, -3.0,0.0},
  {1.0, -1.0,0.0},{1.0, -2.0,0.0},{1.0, -3.0,0.0},
  {2.0, -1.0,0.0},{2.0, -2.0,0.0},{2.0, -3.0,0.0},
  {3.0, -1.0,0.0},{3.0, -2.0,0.0},{3.0, -3.0,0.0},
  {1.0, 0.0, 0.0},{2.0,  0.0,0.0},{3.0,  0.0,0.0},
  {1.0, 1.0, 0.0},{2.0,  1.0,0.0},{3.0,  1.0,0.0},
  {1.0, 2.0, 0.0},{2.0,  2.0,0.0},{3.0,  2.0,0.0},
  {1.0, 3.0, 0.0},{2.0,  3.0,0.0},{3.0,  3.0,0.0},
  {0.0, 1.0, 0.0},{0.0,  2.0,0.0},{0.0,  3.0,0.0},
  {0.0, 1.0, 1.0},{0.0,  2.0,1.0},{0.0,  3.0,1.0},
  {0.0, 1.0, 2.0},{0.0,  2.0,2.0},{0.0,  3.0,2.0},
  {0.0, 1.0, 3.0},{0.0,  2.0,3.0},{0.0,  3.0,3.0},
  {0.0, 0.0, 1.0},{0.0,  0.0,2.0},{0.0,  0.0,3.0},
  {-1.0,0.0, 1.0},{-1.0, 0.0,2.0},{-1.0, 0.0,3.0},
  {-2.0,0.0, 1.0},{-2.0, 0.0,2.0},{-2.0, 0.0,3.0},
  {-3.0,0.0, 1.0},{-3.0, 0.0,2.0},{-3.0, 0.0,3.0},
  {-1.0,0.0, 0.0},{-2.0, 0.0,0.0},{-3.0, 0.0,0.0},
  {-1.0,-1.0,0.0},{-2.0,-1.0,0.0},{-3.0,-1.0,0.0},
  {-1.0,-2.0,0.0},{-2.0,-2.0,0.0},{-3.0,-2.0,0.0},
  {-1.0,-3.0,0.0},{-2.0,-3.0,0.0},{-3.0,-3.0,0.0}};

point3d hdp6[73] =
 {{0.0, 0.0,0.0},
  {0.0, 0.0,-1.0},{0.0, 0.0,-2.0},{0.0, 0.0,-3.0},
  {1.0, 0.0,-1.0},{1.0, 0.0,-2.0},{1.0, 0.0,-3.0},
  {2.0, 0.0,-1.0},{2.0, 0.0,-2.0},{2.0, 0.0,-3.0},
  {3.0, 0.0,-1.0},{3.0, 0.0,-2.0},{3.0, 0.0,-3.0},
  {1.0, 0.0, 0.0},{2.0,  0.0,0.0},{3.0,  0.0,0.0},
  {1.0, 1.0, 0.0},{2.0,  1.0,0.0},{3.0,  1.0,0.0},
  {1.0, 2.0, 0.0},{2.0,  2.0,0.0},{3.0,  2.0,0.0},
  {1.0, 3.0, 0.0},{2.0,  3.0,0.0},{3.0,  3.0,0.0},
  {0.0, 1.0, 0.0},{0.0,  2.0,0.0},{0.0,  3.0,0.0},
  {0.0, 1.0, 1.0},{0.0,  2.0,1.0},{0.0,  3.0,1.0},
  {0.0, 1.0, 2.0},{0.0,  2.0,2.0},{0.0,  3.0,2.0},
  {0.0, 1.0, 3.0},{0.0,  2.0,3.0},{0.0,  3.0,3.0},
  {0.0, 0.0, 1.0},{0.0,  0.0,2.0},{0.0,  0.0,3.0},
  {-1.0,0.0, 1.0},{-1.0, 0.0,2.0},{-1.0, 0.0,3.0},
  {-2.0,0.0, 1.0},{-2.0, 0.0,2.0},{-2.0, 0.0,3.0},
  {-3.0,0.0, 1.0},{-3.0, 0.0,2.0},{-3.0, 0.0,3.0},
  {-1.0,0.0, 0.0},{-2.0, 0.0,0.0},{-3.0, 0.0,0.0},
  {-1.0,-1.0,0.0},{-2.0,-1.0,0.0},{-3.0,-1.0,0.0},
  {-1.0,-2.0,0.0},{-2.0,-2.0,0.0},{-3.0,-2.0,0.0},
  {-1.0,-3.0,0.0},{-2.0,-3.0,0.0},{-3.0,-3.0,0.0},
  { 0.0,-1.0,0.0},{ 0.0,-2.0,0.0},{ 0.0,-3.0,0.0},
  {0.0,-1.0,-1.0},{0.0,-2.0,-1.0},{0.0,-3.0,-1.0},
  {0.0,-1.0,-2.0},{0.0,-2.0,-2.0},{0.0,-3.0,-2.0},
  {0.0,-1.0,-3.0},{0.0,-2.0,-3.0},{0.0,-3.0,-3.0}};

#define  s30  (2.0*0.5)
#define  js30 (2.0-s30)
#define  c30  (2.0*0.866025403)
#define  jc30 (2.0-c30)
#define  t30  (2.0*0.577350269)  /* 2 tan(30 deg) */
#define  jt30 (2.0-t30)

point3d hdp8[97] =
 {{0.0,0.0,0.0},
  {-s30,-s30,jc30},{-c30,-c30,js30},{-2,-2,2},
  {-s30,-3,jc30},{-c30,-3,js30},{-2,-3,2},
  {-s30,-4,jc30},{-c30,-4,js30},{-2,-4,2},
  {-s30,-5,jc30},{-c30,-5,js30},{-2,-5,2},
  {   0,-3,   0},{   0,-4,   0},{   0,-5,   0},
  { s30,-3,jc30},{ s30,-4,jc30},{ s30,-5,jc30},
  { c30,-3,js30},{ c30,-4,js30},{ c30,-5,js30},
  {   2,-3,   2},{   2,-4,   2},{   2,-5,   2},
  { s30,-s30,jc30},{ c30,-c30,js30},{ 2,-2,2},
  {   3,-s30,jc30},{   3,-c30,js30},{ 3,-2,2},
  {   4,-s30,jc30},{   4,-c30,js30},{ 4,-2,2},
  {   5,-s30,jc30},{   5,-c30,js30},{ 5,-2,2},
  {   3,   0,   0},{   4, 0,   0},{   5, 0,   0},
  {   3,s30,jc30},{   4,s30,jc30},{   5,s30,jc30},
  {   3,c30,js30},{   4,c30,js30},{   5,c30,js30},
  {   3, 2,   2},{   4, 2,   2},{   5, 2,   2},
  { s30, s30,jc30},{ c30, c30,js30},{ 2, 2,2},
  { s30, 3,jc30},{ c30, 3,js30},{ 2, 3,2},
  { s30, 4,jc30},{ c30, 4,js30},{ 2, 4,2},
  { s30, 5,jc30},{ c30, 5,js30},{ 2, 5,2},
  {   0, 3,   0},{   0, 4,   0},{   0, 5,   0},
  {-s30, 3,jc30},{-s30, 4,jc30},{-s30, 5,jc30},
  {-c30, 3,js30},{-c30, 4,js30},{-c30, 5,js30},
  {  -2, 3,   2},{  -2, 4,   2},{  -2, 5,   2},
  {-s30, s30,jc30},{-c30, c30,js30},{-2, 2,2},
  {  -3, s30,jc30},{  -3, c30,js30},{-3, 2,2},
  {  -4, s30,jc30},{  -4, c30,js30},{-4, 2,2},
  {  -5, s30,jc30},{  -5, c30,js30},{-5, 2,2},
  {  -3,   0,   0},{  -4, 0,   0},{  -5, 0,   0},
  {  -3,-s30,jc30},{  -4,-s30,jc30},{  -5,-s30,jc30},
  {  -3,-c30,js30},{  -4,-c30,js30},{  -5,-c30,js30},
  {  -3,-2,   2},{  -4,-2,   2},{  -5,-2,   2}};

/* ///////////////////////////////////////////////////////////////////////// */
static void CPt ( point3d *p00, point3d *p01, point3d *p10, point3d *p11,
                  point3d *p )
{
  double a, b, c, d, e, f, d0, d1;

  a = p01->x-p00->x;  b = p10->x-p11->x;  e = p10->x-p00->x;
  c = p01->y-p00->y;  d = p10->y-p11->y;  f = p10->y-p00->y;
  d0 = a*d-b*c;
  d1 = e*d-f*b;
  InterPoint3d ( p00, p01, d1/d0, p );
} /*CPt*/

void InitHole ( int kk, double *ip )
{
  int      i, j;
  double    ang, ca, sa;
  vector3d v1, v2;
  trans3d  t;
  point3d  *pp, *auxp1, *auxp2;
  double    zf;
  int      ind[16], ind1[16];
  void     *sp;

  sp = pkv_GetScratchMemTop ();
  auxp1 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point3d) );
  auxp2 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point3d) );
  if ( !auxp1 || !auxp2 )
    exit ( 1 );

  kk = min ( GH_MAX_K, kk );
  hole_k = kk = max ( 3, kk );
  SetPoint3d ( &auxp1[0], 0.0, 0.0, 0.0 );
  ang = 2.0*PI/hole_k;
  sa = sin ( ang );
  ca = cos ( ang );
  SetVector3d ( &v1, 0.0, -1.0, 0.0 );
  SetVector3d ( &v2, sa, -ca, 0.0 );
  for ( i = 1; i <= 3; i++ ) {
    AddVector3d ( &auxp1[i-1], &v1, &auxp1[i] );
    for ( j = 1; j <= 3; j++ )
      AddVector3d ( &auxp1[3*j-3+i], &v2, &auxp1[3*j+i] );
  }
  IdentTrans3d ( &t );
  RotZTrans3d ( &t, ang );
  for ( i = 1; i < hole_k; i++ ) {
    for ( j = 1; j <= 12; j++ )
      TransPoint3d ( &t, &auxp1[(i-1)*12+j], &auxp1[i*12+j] );
  }
  hole_np = 1+12*hole_k;
  memcpy ( auxp2, auxp1, hole_np*sizeof(point3d) );

  for ( i = 0; i < hole_k; i++ )
    MidPoint3d ( &auxp2[12*(i+1)], &auxp2[12*((i+1)%hole_k+1)],
                 &auxp2[12*((i+1)%hole_k)+3] );
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    InterPoint3d ( &auxp2[ind[12]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[13]] );
    InterPoint3d ( &auxp2[ind[12]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[14]] );
    InterPoint3d ( &auxp2[ind[3]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[7]] );
    InterPoint3d ( &auxp2[ind[3]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[11]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    MidPoint3d ( &auxp2[ind[4]], &auxp2[ind1[14]], &auxp2[ind[7]] );
    MidPoint3d ( &auxp2[ind[8]], &auxp2[ind1[13]], &auxp2[ind[11]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    CPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[7]], &auxp2[ind[4]],
          &auxp2[ind[6]] );
    CPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[11]], &auxp2[ind[8]],
          &auxp2[ind[10]] );
    CPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[7]], &auxp2[ind[4]],
          &auxp2[ind[5]] );
    CPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[11]], &auxp2[ind[8]],
          &auxp2[ind[9]] );
  }
  for ( i = 0; i < hole_np; i++ )
    InterPoint3d ( &auxp1[i], &auxp2[i], ip[1], &auxp1[i] );

  for ( i = 0; i < hole_np; i++ )
    auxp1[i].z = ip[2]*(sqrt(auxp1[i].x*auxp1[i].x+auxp1[i].y*auxp1[i].y)-2.0);

  zf = sqrt((2.0-ip[0])*ip[0]);
  switch ( hole_k ) {
case 3: pp = hdp3;  break;
case 5: pp = hdp5;  break;
case 6: pp = hdp6;  break;
case 8: pp = hdp8;  break;
default:
    exit ( 1 );
  }
  for ( i = 0; i < hole_np; i++ )
    SetPoint3d ( &Bpt[i], (1.0-ip[0])*auxp1[i].x+ip[0]*pp[i].x,
                          (1.0-ip[0])*auxp1[i].y+ip[0]*pp[i].y,
                          (1.0-ip[0])*auxp1[i].z+zf*pp[i].z );

  pkv_SetScratchMemTop ( sp );
} /*InitHole*/

/* ///////////////////////////////////////////////////////////////////////// */
static void DCPt ( point2d *p00, point2d *p01, point2d *p10, point2d *p11,
                   point2d *p )
{
  double a, b, c, d, e, f, d0, d1;

  a = p01->x-p00->x;  b = p10->x-p11->x;  e = p10->x-p00->x;
  c = p01->y-p00->y;  d = p10->y-p11->y;  f = p10->y-p00->y;
  d0 = a*d-b*c;
  d1 = e*d-f*b;
  InterPoint2d ( p00, p01, d1/d0, p );
} /*DCPt*/

void InitDomain ( int kk, double *ip )
{
  int      i, j;
  double    ang, ca, sa;
  vector2d v1, v2;
  trans2d  t;
  point2d  *auxp1, *auxp2;
  int      ind[16], ind1[16];
  void     *sp;

  sp = pkv_GetScratchMemTop ();
  auxp1 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point2d) );
  auxp2 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point2d) );
  if ( !auxp1 || !auxp2 )
    exit ( 1 );

  kk = min ( GH_MAX_K, kk );
  hole_k = kk = max ( 3, kk );
  SetPoint2d ( &auxp1[0], 0.0, 0.0 );
  ang = 2.0*PI/hole_k;
  sa = sin ( ang );
  ca = cos ( ang );
  SetVector2d ( &v1, 0.0, -1.0 );
  SetVector2d ( &v2, sa, -ca );
  for ( i = 1; i <= 3; i++ ) {
    AddVector2d ( &auxp1[i-1], &v1, &auxp1[i] );
    for ( j = 1; j <= 3; j++ )
      AddVector2d ( &auxp1[3*j-3+i], &v2, &auxp1[3*j+i] );
  }
  IdentTrans2d ( &t );
  RotTrans2d ( &t, ang );
  for ( i = 1; i < hole_k; i++ ) {
    for ( j = 1; j <= 12; j++ )
      TransPoint2d ( &t, &auxp1[(i-1)*12+j], &auxp1[i*12+j] );
  }
  hole_np = 1+12*hole_k;
  memcpy ( auxp2, auxp1, hole_np*sizeof(point2d) );

  for ( i = 0; i < hole_k; i++ )
    MidPoint2d ( &auxp2[12*(i+1)], &auxp2[12*((i+1)%hole_k+1)],
                 &auxp2[12*((i+1)%hole_k)+3] );
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    InterPoint2d ( &auxp2[ind[12]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[13]] );
    InterPoint2d ( &auxp2[ind[12]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[14]] );
    InterPoint2d ( &auxp2[ind[3]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[7]] );
    InterPoint2d ( &auxp2[ind[3]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[11]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    MidPoint2d ( &auxp2[ind[4]], &auxp2[ind1[14]], &auxp2[ind[7]] );
    MidPoint2d ( &auxp2[ind[8]], &auxp2[ind1[13]], &auxp2[ind[11]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    DCPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[7]], &auxp2[ind[4]],
           &auxp2[ind[6]] );
    DCPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[11]], &auxp2[ind[8]],
           &auxp2[ind[10]] );
    DCPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[7]], &auxp2[ind[4]],
           &auxp2[ind[5]] );
    DCPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[11]], &auxp2[ind[8]],
           &auxp2[ind[9]] );
  }
  for ( i = 0; i < hole_np; i++ )
    InterPoint2d ( &auxp1[i], &auxp2[i], ip[1], &Dompt[i] );

  pkv_SetScratchMemTop ( sp );
} /*InitDomain*/

