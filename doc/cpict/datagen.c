
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"

#include "datagen.h"


/* built-in data generator */

int hole_k = 5;              /* number of hole sides, 3, 5, 6 or 8 */
int hole_np;                 /* number of B-spline control points */
point3f Bpt[MAX_BPTS];       /* B-spline control points */

float dataparam[3] =         /* data initialization parameters */
  { 0.0, 0.0, 0.0 };


/* displaying switches */
boolean netk3 = false;
boolean netk5 = true;
boolean netk6 = false;
boolean netk8 = false;


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

void GetBspInd ( int i, int j, int *ind )
{
  int  l;
  int  ind1[16];

  switch ( j ) {
case 0:
    i = i % hole_k;
    ind[0] = 0;
    for ( j = 1; j <= 3; j++ ) {
      for ( l = 0; l <= 3; l++ )
        ind[j+4*l] = 12*i+3*l+j;
    }
    i = (i+1) % hole_k;
    for ( j = 1; j <= 3; j++ )
      ind[4*j] = 12*i+j;
    break;

case 1:
    GetBspInd ( i, 0, ind );
    GetBspInd ( i+1, 0, ind1 );
    for ( i = 0; i <= 3; i++ ) {
      memmove ( &ind[4*i+1], &ind[4*i], 3*sizeof(int) );
      ind[4*i] = ind1[4+i];
    }
    break;

case 2:
    GetBspInd ( i, 0, ind );
    GetBspInd ( i+1, 0, ind1 );
    for ( i = 0; i <= 3; i++ ) {
      memmove ( &ind[4*i+2], &ind[4*i], 2*sizeof(int) );
      ind[4*i]   = ind1[8+i];
      ind[4*i+1] = ind1[4+i];
    }
    break;

default:
    exit ( 1 );
  }
} /*GetBspInd*/

void InitHole ( int kk, float ip1, float ip2, float ip3 )
{
  int      i, j;
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

  kk = min ( MAX_K, kk );
  hole_k = kk = max ( 3, kk );
  SetPoint3f ( &auxp1[0], 0.0, 0.0, 0.0 );
  ang = (float)(2.0*PI/hole_k);
  sa = (float)sin ( ang );
  ca = (float)cos ( ang );
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
    GetBspInd ( i, 0, ind );
    InterPoint3f ( &auxp2[ind[12]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[13]] );
    InterPoint3f ( &auxp2[ind[12]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[14]] );
    InterPoint3f ( &auxp2[ind[3]], &auxp2[ind[15]], 1/3.0, &auxp2[ind[7]] );
    InterPoint3f ( &auxp2[ind[3]], &auxp2[ind[15]], 2/3.0, &auxp2[ind[11]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    GetBspInd ( i, 0, ind );
    GetBspInd ( i+1, 0, ind1 );
    MidPoint3f ( &auxp2[ind[7]], &auxp2[ind1[13]], &auxp2[ind[4]] );
    MidPoint3f ( &auxp2[ind[11]], &auxp2[ind1[14]], &auxp2[ind[8]] );
  }
  for ( i = 0; i < hole_k; i++ ) {
    GetBspInd ( i, 0, ind );
    CPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[4]], &auxp2[ind[7]],
          &auxp2[ind[5]] );
    CPt ( &auxp2[ind[1]], &auxp2[ind[13]], &auxp2[ind[8]], &auxp2[ind[11]],
          &auxp2[ind[9]] );
    CPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[4]], &auxp2[ind[7]],
          &auxp2[ind[6]] );
    CPt ( &auxp2[ind[2]], &auxp2[ind[14]], &auxp2[ind[8]], &auxp2[ind[11]],
          &auxp2[ind[10]] );
  }
  for ( i = 0; i < hole_np; i++ )
    InterPoint3f ( &auxp1[i], &auxp2[i], ip2, &auxp1[i] );

  for ( i = 0; i < hole_np; i++ )
    auxp1[i].z = (float)(ip3*(sqrt(auxp1[i].x*auxp1[i].x+auxp1[i].y*auxp1[i].y)-2.0));

  zf = (float)sqrt((2.0-ip1)*ip1);
  switch ( hole_k ) {
case 3: pp = hdp3;  break;
case 5: pp = hdp5;  break;
case 6: pp = hdp6;  break;
case 8: pp = hdp8;  break;
default:
    exit ( 1 );
  }
  for ( i = 0; i < hole_np; i++ )
    SetPoint3f ( &Bpt[i], (float)((1.0-ip1)*auxp1[i].x+ip1*pp[i].x),
                          (float)((1.0-ip1)*auxp1[i].y+ip1*pp[i].y),
                          (float)((1.0-ip1)*auxp1[i].z+zf*pp[i].z) );

  pkv_SetScratchMemTop ( sp );
} /*InitHole*/

