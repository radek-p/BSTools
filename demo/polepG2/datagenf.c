
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "eg2holef.h"
#include "datagenf.h"


/* built-in data generator */

static int hole_k;           /* number of hole sides, 3, 5, 6 or 8 */
static boolean fdomcp_valid = false;
static point2f fdomcp[MAX_BPTS];

/* ////////////////////////////////////////////////////////////////////////// */
void InitKnots ( int kk, float *knots )
{
  int i;

  knots[0] = knots[1] = 0.0;
  for ( i = 2; i < 9; i++ )
    knots[i] = (float)(i-1)/8.0;
  knots[9] = knots[10] = 1.0;
  for ( i = 1; i < kk; i++ )
    memcpy ( &knots[11*i], &knots[0], 11*sizeof(float) );
} /*InitKnots*/

/* ////////////////////////////////////////////////////////////////////////// */
point3f hdp3[37] =
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

point3f hdp5[61] =
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

point3f hdp6[73] =
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

point3f hdp6a[73] =
 {{0.0,0.0,0.0},
  {-s30,-s30,jc30},{-c30,-c30,js30},{-2,-2,2},        /*1-3*/
  {-s30,-3,jc30},{-c30,-3,js30},{-2,-3,2},            /*4-6*/
  {-s30,-4,jc30},{-c30,-4,js30},{-2,-4,2},            /*7-9*/
  {-s30,-5,jc30},{-c30,-5,js30},{-2,-5,2},            /*10-12*/
  {   0,-3,   0},{   0,-4,   0},{   0,-5,   0},       /*13-15*/
  { s30,-3,jc30},{ s30,-4,jc30},{ s30,-5,jc30},       /*16-18*/
  { c30,-3,js30},{ c30,-4,js30},{ c30,-5,js30},       /*19-21*/
  {   2,-3,   2},{   2,-4,   2},{   2,-5,   2},       /*22-24*/
  { s30,-s30,jc30},{ c30,-c30,js30},{ 2,-2,2},        /*25-27*/
  {   3,-s30,jc30},{   3,-c30,js30},{ 3,-2,2},        /*28-30*/
  {   4,-s30,jc30},{   4,-c30,js30},{ 4,-2,2},        /*31-33*/
  {   5,-s30,jc30},{   5,-c30,js30},{ 5,-2,2},        /*34-36*/
  {   3,   0,   0},{   4, 0,   0},{   5, 0,   0},     /*37-39*/
  {   3,s30,jc30},{   4,s30,jc30},{   5,s30,jc30},    /*40-42*/
  {   3,c30,js30},{   4,c30,js30},{   5,c30,js30},    /*43-45*/
  {   3, 2,   2},{   4, 2,   2},{   5, 2,   2},       /*46-48*/
  {   0, s30,jc30},{ 0, c30,js30},{ 0, 2,2},          /*49-51*/
  {  -3, s30,jc30},{  -3, c30,js30},{-3, 2,2},        /*52-54*/
  {  -4, s30,jc30},{  -4, c30,js30},{-4, 2,2},        /*55-57*/
  {  -5, s30,jc30},{  -5, c30,js30},{-5, 2,2},        /*58-60*/
  {  -3,   0,   0},{  -4, 0,   0},{  -5, 0,   0},     /*61-63*/
  {  -3,-s30,jc30},{  -4,-s30,jc30},{  -5,-s30,jc30}, /*64-66*/
  {  -3,-c30,js30},{  -4,-c30,js30},{  -5,-c30,js30}, /*67-69*/
  {  -3,-2,   2},{  -4,-2,   2},{  -5,-2,   2}};      /*70-72*/

point3f hdp8[97] =
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
static void CPt ( point3f *p00, point3f *p01, point3f *p10, point3f *p11,
                  point3f *p )
{
  float a, b, c, d, e, f, d0, d1;

  a = p01->x-p00->x;  b = p10->x-p11->x;  e = p10->x-p00->x;
  c = p01->y-p00->y;  d = p10->y-p11->y;  f = p10->y-p00->y;
  d0 = a*d-b*c;
  d1 = e*d-f*b;
  InterPoint3f ( p00, p01, d1/d0, p );
} /*CPt*/

int InitHole ( int kk, float *ip, point3f *surfcp )
{
  int      i, j, hole_np;
  float    ang, ca, sa;
  vector3f v1, v2;
  trans3f  t;
  point3f  *pp, *auxp1, *auxp2;
  float    zf;
  int      ind[16], ind1[16];
  void     *sp;

  sp = pkv_GetScratchMemTop ();
  auxp1 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point3f) );
  auxp2 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point3f) );
  if ( !auxp1 || !auxp2 )
    exit ( 1 );

  kk = min ( GH_MAX_K, kk );
  hole_k = kk = max ( 3, kk );
  SetPoint3f ( &auxp1[0], 0.0, 0.0, 0.0 );
  ang = 2.0*PI/hole_k;
  sa = sin ( ang );
  ca = cos ( ang );
  SetVector3f ( &v1, 0.0, -1.0, 0.0 );
  SetVector3f ( &v2, sa, -ca, 0.0 );
  for ( i = 1; i <= 3; i++ ) {
    AddVector3f ( &auxp1[i-1], &v1, &auxp1[i] );
    for ( j = 1; j <= 3; j++ )
      AddVector3f ( &auxp1[3*j-3+i], &v2, &auxp1[3*j+i] );
  }
  IdentTrans3f ( &t );
  RotZTrans3f ( &t, ang );
  for ( i = 1; i < hole_k; i++ ) {
    for ( j = 1; j <= 12; j++ )
      TransPoint3f ( &t, &auxp1[(i-1)*12+j], &auxp1[i*12+j] );
  }
  hole_np = 1+12*hole_k;
  memcpy ( auxp2, auxp1, hole_np*sizeof(point3f) );

  for ( i = 0; i < hole_k; i++ )
    MidPoint3f ( &auxp2[12*(i+1)], &auxp2[12*((i+1)%hole_k+1)],
                 &auxp2[12*((i+1)%hole_k)+3] );
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    InterPoint3f ( &auxp2[ind[12]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[13]] );
    InterPoint3f ( &auxp2[ind[12]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[14]] );
    InterPoint3f ( &auxp2[ind[0]], &auxp2[ind[12]], 1/3.0, &auxp2[ind[4]] );
    InterPoint3f ( &auxp2[ind[0]], &auxp2[ind[12]], 2/3.0, &auxp2[ind[8]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    MidPoint3f ( &auxp2[ind[4]], &auxp2[ind1[14]], &auxp2[ind[7]] );
    MidPoint3f ( &auxp2[ind[8]], &auxp2[ind1[13]], &auxp2[ind[11]] );
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
    InterPoint3f ( &auxp1[i], &auxp2[i], ip[1], &auxp1[i] );

  for ( i = 0; i < hole_np; i++ )
    auxp1[i].z = ip[2]*(sqrt(auxp1[i].x*auxp1[i].x+auxp1[i].y*auxp1[i].y)-2.0);

  if ( hole_k == 6 ) {
    for ( i = 0; i < hole_np; i++ )
      InterPoint3f ( &auxp1[i], &hdp6a[i], ip[3], &auxp1[i] );
  }

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
    SetPoint3f ( &surfcp[i], (1.0-ip[0])*auxp1[i].x+ip[0]*pp[i].x,
                             (1.0-ip[0])*auxp1[i].y+ip[0]*pp[i].y,
                             (1.0-ip[0])*auxp1[i].z+zf*pp[i].z );

  pkv_SetScratchMemTop ( sp );
  return hole_np;
} /*InitHole*/

/* ///////////////////////////////////////////////////////////////////////// */
static void DCPt ( point2f *p00, point2f *p01, point2f *p10, point2f *p11,
                   point2f *p )
{
  float a, b, c, d, e, f, d0, d1;

  a = p01->x-p00->x;  b = p10->x-p11->x;  e = p10->x-p00->x;
  c = p01->y-p00->y;  d = p10->y-p11->y;  f = p10->y-p00->y;
  d0 = a*d-b*c;
  d1 = e*d-f*b;
  InterPoint2f ( p00, p01, d1/d0, p );
} /*DCPt*/

int InitDomain ( int kk, vector2f *part, float *ip, point2f *domcp )
{
  int      i, j, k, l, hole_np;
  float    a;
  point2f  *auxp1, *auxp2;
  int      ind[16], ind1[16];
  void     *sp;

  if ( fdomcp_valid && kk != hole_k )
    fdomcp_valid = false;

  sp = pkv_GetScratchMemTop ();
  auxp1 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point2f) );
  auxp2 = pkv_GetScratchMem ( MAX_BPTS*sizeof(point2f) );
  if ( !auxp1 || !auxp2 )
    exit ( 1 );

  kk = min ( GH_MAX_K, kk );
  hole_k = kk = max ( 3, kk );
  SetPoint2f ( &auxp1[0], 0.0, 0.0 );
  for ( k = 0, l = hole_k-1;  k < hole_k;  l = k, k++ )
    for ( i = 1; i <= 3; i++ ) {
      MultVector2f ( (float)i/3.0, &part[l], &auxp1[12*k+i] );
      for ( j = 1; j <= 3; j++ )
        AddVector2Mf ( &auxp1[12*k+3*(j-1)+i], &part[k], 1.0/3.0, &auxp1[12*k+3*j+i] );
    }
  hole_np = 1+12*hole_k;
  memcpy ( auxp2, auxp1, hole_np*sizeof(point2f) );

  for ( i = 0; i < hole_k; i++ )
    MidPoint2f ( &auxp2[12*(i+1)], &auxp2[12*((i+1) % hole_k+1)],
                 &auxp2[12*((i+1) % hole_k)+3] );

  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    InterPoint2f ( &auxp2[ind[12]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[13]] );
    InterPoint2f ( &auxp2[ind[12]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[14]] );
    InterPoint2f ( &auxp2[ind[0]], &auxp2[ind[12]], 1/3.0, &auxp2[ind[4]] );
    InterPoint2f ( &auxp2[ind[0]], &auxp2[ind[12]], 2/3.0, &auxp2[ind[8]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    MidPoint2f ( &auxp2[ind[4]], &auxp2[ind1[14]], &auxp2[ind[7]] );
    MidPoint2f ( &auxp2[ind[8]], &auxp2[ind1[13]], &auxp2[ind[11]] );
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

  domcp[0] = auxp1[0];
  for ( k = 0; k < hole_k; k++ )
    for ( i = 1; i <= 12; i++ ) {
      if ( i == 1 || i == 4 )
        a = ip[1];
      else
        a = ip[0];
      InterPoint2f ( &auxp1[12*k+i], &auxp2[12*k+i], a, &domcp[12*k+i] );
    }
  if ( fdomcp_valid ) {
    for ( i = 0; i <= 12*hole_k; i++ )
      InterPoint2f ( &domcp[i], &fdomcp[i], ip[2], &domcp[i] );
  }
  pkv_SetScratchMemTop ( sp );
  return hole_np;
} /*InitDomain*/

void StoreDomain ( int k, point2f *domcp )
{
  hole_k = k;
  memcpy ( fdomcp, domcp, (12*hole_k+1)*sizeof(point2f) );
  fdomcp_valid = true;
} /*StoreDomain*/

