
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"

#include "datagenf.new.h"

void InitGHKnotsf ( int hole_k, float *knots )
{
  int i;

  knots[0] = knots[1] = 0.0;
  for ( i = 2; i < 9; i++ )
    knots[i] = (i-1)/8.0;
  knots[9] = knots[10] = 1.0;
  for ( i = 1; i < hole_k; i++ )
    memcpy ( &knots[11*i], &knots[0], 11*sizeof(float) );
} /*InitGHKnotsf*/

void InitGHVectors2f ( int hole_k, vector2f *v )
{
  int i;
  trans2f tr;
  float  a;

  a = 2.0*PI/(float)hole_k;
  IdentTrans2f ( &tr );
  RotTrans2f ( &tr, a );
  SetVector2f ( &v[0], 3.0*cos(a), -3.0*sin(a) );
  for ( i = 1; i < hole_k; i++ )
    TransVector2f ( &tr, &v[i-1], &v[i] );
} /*InitGHVectors2f*/

void InitGHVectors3f ( int hole_k, vector3f *v )
{
  int     i;
  trans3f tr;
  float  a;

  a = 2.0*PI/(float)hole_k;
  IdentTrans3f ( &tr );
  RotZTrans3f ( &tr, a );
  SetVector3f ( &v[0], 3.0*cos(a), -3.0*sin(a), 0.0 );
  for ( i = 1; i < hole_k; i++ )
    TransVector3f ( &tr, &v[i-1], &v[i] );
} /*InitGHVectors3f*/

void LineInters2f ( const point2f *p0, const point2f *p1,
                    const point2f *q0, const point2f *q1,
                    point2f *p )
{
  float a[6];

  a[0] = p1->x-p0->x;  a[1] = q0->x-q1->x;
  a[2] = p1->y-p0->y;  a[3] = q0->y-q1->y;
  a[4] = q0->x-p0->x;  a[5] = q0->y-p0->y;
  pkn_multiGaussSolveLinEqf ( 2, a, 1, 1, &a[4] );
  InterPoint2f ( p0, p1, a[4], p );
} /*LineInters2f*/

void LineInters3f ( const point3f *p0, const point3f *p1,
                    const point3f *q0, const point3f *q1,
                    point3f *p )
{
#define TOL 1.0e-8
  float   a[12];
  vector3f u, v, w;

  SubtractPoints3f ( p1, p0, &u );
  SubtractPoints3f ( q1, q0, &v );
  CrossProduct3f ( &u, &v, &w );
  if ( fabs(w.x)+fabs(w.y)+fabs(w.z) >= TOL ) {
    a[0] = u.x;  a[1] = -v.x;  a[2] = -w.x;
    a[3] = u.y;  a[4] = -v.y;  a[5] = -w.y;
    a[6] = u.z;  a[7] = -v.z;  a[8] = -w.z;
    a[9] = q0->x-p0->x;  a[10] = q0->y-p0->y;  a[11] = q0->z-p0->z;
    pkn_multiGaussSolveLinEqf ( 3, a, 1, 1, &a[9] );
    InterPoint3f ( p0, p1, a[9], &u );
    InterPoint3f ( q0, q1, a[10], &v );
  }
  else {
    MidPoint3f ( p0, p1, &u );
    MidPoint3f ( q0, q1, &v );
  }
  MidPoint3f ( &u, &v, p );
#undef TOL
} /*LineInters3f*/

static void SetAuxCP2 ( point2f *p0, point2f *p1, point2f *p2, point2f *p3,
                        point2f *p4, point2f *p5, point2f *p6,
                        float phi, int j, point2f *cp )
{
  vector3f vxi, veta, vzeta;
  float a[9], ai[9], x[3], r;
  int    P[2], Q[2];

  r = 1.0/sin(phi);
  SetVector3f ( &vxi, p1->x-p0->x, p1->y-p0->y, 0.0 );
  NormalizeVector3f ( &vxi );
  SetVector3f ( &veta, p3->x-p2->x, p3->y-p2->y, 0.0 );
  CrossProduct3f ( &veta, &vxi, &vzeta );
  MultVector3f ( 0.5, &vzeta, &vzeta );
  SetVector3f ( &veta, p5->x-p4->x, p5->y-p4->y, 0.0 );
  a[0] = vxi.x;  a[1] = veta.x;  a[2] = vzeta.x;
  a[3] = vxi.y;  a[4] = veta.y;  a[5] = vzeta.y;
  a[6] = vxi.z;  a[7] = veta.z;  a[8] = vzeta.z;
  memcpy ( ai, a, 9*sizeof(float) );
  pkn_GaussDecomposePLUQf ( 3, ai, P, Q );
  memcpy ( x, p6, 3*sizeof(float) );
  pkn_multiSolvePLUQf ( 3, ai, P, Q, 1, 1, x );
  AddVector2Mf ( p4, (vector2f*)&veta, r*sin((float)j/3.0*phi), cp );
  AddVector2Mf ( cp, (vector2f*)&vzeta, r*(cos((float)j/3.0*phi)-cos(phi)), cp );
} /*SetAuxCP2*/

int InitGHDomainNetf ( int hole_k, const vector2f *v,
                       int nparams, const float *param,
                       point2f *domcp )
{
  void     *sp;
  int      i, j, k, l, m, npoints;
  point2f  *acp1, *acp2, *acp3;
  vector2f av;
  float   phi;

  sp = pkv_GetScratchMemTop ();
  npoints = 1+12*hole_k;

        /* rhomboidal */
  SetPoint2f ( &domcp[0], 0.0, 0.0 );
  for ( i = 0, k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;
    for ( l = 1; l < 4; l++ ) {
      AddVector2Mf ( &domcp[0], &v[i], ((float)l)/3.0, &domcp[k+l] );
      for ( m = 1; m < 4; m++ )
        AddVector2Mf ( &domcp[k+l], &v[j], ((float)m)/3.0, &domcp[k+l+3*m] );
    }    
  }
  if ( !(acp1 = (point2f*)pkv_GetScratchMem ( 3*npoints*sizeof(point2f) )) )
    goto way_out;
  acp2 = &acp1[npoints];
  acp3 = &acp2[npoints];

        /* 1st modification */
  memcpy ( acp1, domcp, npoints*sizeof(point2f) );
  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    MidPoint2f ( &acp1[12*i+12], &acp1[12*j+12], &acp1[12*j+3] );
    InterPoint2f ( &acp1[0], &acp1[12*j+3], 1.0/3.0, &acp1[12*j+1] );
    InterPoint2f ( &acp1[0], &acp1[12*j+3], 2.0/3.0, &acp1[12*j+2] );
    InterPoint2f ( &acp1[12*i+12], &acp1[12*j+3], 1.0/3.0, &acp1[12*i+11] );
    InterPoint2f ( &acp1[12*i+12], &acp1[12*j+3], 2.0/3.0, &acp1[12*i+10] );
    InterPoint2f ( &acp1[12*j+3], &acp1[12*j+12], 1.0/3.0, &acp1[12*j+6] );
    InterPoint2f ( &acp1[12*j+3], &acp1[12*j+12], 2.0/3.0, &acp1[12*j+9] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2f ( &acp1[k+2], &acp1[k+11], &acp1[m+2], &acp1[k+9], &acp1[k+8] );
    LineInters2f ( &acp1[k+2], &acp1[k+11], &acp1[m+1], &acp1[k+6], &acp1[k+5] );
    LineInters2f ( &acp1[k+1], &acp1[k+10], &acp1[m+2], &acp1[k+9], &acp1[k+7] );
    LineInters2f ( &acp1[k+1], &acp1[k+10], &acp1[m+1], &acp1[k+6], &acp1[k+4] );
  }
  pkn_SubtractMatrixf ( 1, 2*npoints, 0, (float*)acp1, 0, (float*)domcp,
                        0, (float*)acp1 );
  pkn_MultMatrixNumf ( 1, 2*npoints, 0, (float*)acp1, param[0],
                       0, (float*)acp1 );
        /* 2nd modification */
  memcpy ( acp2, domcp, npoints*sizeof(point2f) );
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    MidPoint2f ( &acp2[k+12], &acp2[m+12], &acp2[m+3] );
    InterPoint2f ( &acp2[k+12], &acp2[m+3], 1.0/3.0, &acp2[k+11] );
    InterPoint2f ( &acp2[k+12], &acp2[m+3], 2.0/3.0, &acp2[k+10] );
    InterPoint2f ( &acp2[m+3], &acp2[m+12], 1.0/3.0, &acp2[m+6] );
    InterPoint2f ( &acp2[m+3], &acp2[m+12], 2.0/3.0, &acp2[m+9] );
  }
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2f ( &acp2[0], &acp2[m+3], &acp2[k+6], &acp2[m+10], &acp2[m+1] );
    LineInters2f ( &acp2[0], &acp2[m+3], &acp2[k+9], &acp2[m+11], &acp2[m+2] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2f ( &acp2[k+2], &acp2[k+11], &acp2[m+2], &acp2[k+9], &acp2[k+8] );
    LineInters2f ( &acp2[k+2], &acp2[k+11], &acp2[m+1], &acp2[k+6], &acp2[k+5] );
    LineInters2f ( &acp2[k+1], &acp2[k+10], &acp2[m+2], &acp2[k+9], &acp2[k+7] );
    LineInters2f ( &acp2[k+1], &acp2[k+10], &acp2[m+1], &acp2[k+6], &acp2[k+4] );
  }
  pkn_SubtractMatrixf ( 1, 2*npoints, 0, (float*)acp2, 0, (float*)domcp,
                        0, (float*)acp2 );
  pkn_AddMatrixMf ( 1, 2*npoints, 0, (float*)acp1, 0, (float*)acp2,
                    param[1], 0, (float*)acp1 );
  pkn_AddMatrixf ( 1, 2*npoints, 0, (float*)domcp, 0, (float*)acp1,
                   0, (float*)acp3 );

        /* 3rd modification */
  phi = 0.5*param[2]*PI;
  if ( phi > 0.5*PI ) phi = 0.5*PI;
  if ( phi > 0.0 ) {
    memset ( acp1, 0, npoints*sizeof(point2f) );
    memcpy ( acp2, acp3, npoints*sizeof(point2f) );
    SetVector2f ( &acp2[0], 0.0, 0.0 );
    for ( i = k = 0;  i < hole_k;  i += 2, k += 24 ) {
      j = (i+1) % hole_k;  m = 12*j;
      SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[0], phi, 0, &av );
      AddVector2f ( &acp2[0], &av, &acp2[0] );
      SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[k+1], phi, 1, &acp2[k+1] );
      SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[k+2], phi, 2, &acp2[k+2] );
      for ( l = 1; l < 4; l++ ) {
        SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l],
                   phi, 0, &acp2[m+l] );
        SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l+3],
                   phi, 1, &acp2[m+l+3] );
        SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l+6],
                   phi, 2, &acp2[m+l+6] );
        SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[k+3*l+3], &acp3[k+3*l+1],
                   phi, 1, &acp2[k+3*l+1] );
        SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[k+3*l+3], &acp3[k+3*l+2],
                   phi, 2, &acp2[k+3*l+2] );
      }
    }
    if ( hole_k & 0x01 ) {
      MultVector2f ( 2.0/(float)(hole_k+1), &acp2[0], &av ) ;
      AddVector2f ( &domcp[0], &av, &acp2[0] );
      AddVector2f ( &domcp[12], &av, &acp2[12] );
      InterPoint2f ( &acp2[3], &acp2[12], 1.0/3.0, &acp2[6] );
      InterPoint2f ( &acp2[3], &acp2[12], 2.0/3.0, &acp2[9] );
      InterPoint2f ( &acp2[15], &acp2[12], 1.0/3.0, &acp2[10] );
      InterPoint2f ( &acp2[15], &acp2[12], 2.0/3.0, &acp2[11] );
      LineInters2f ( &acp2[1], &acp2[10], &acp2[6], &acp2[13], &acp2[4] );
      LineInters2f ( &acp2[1], &acp2[10], &acp2[9], &acp2[14], &acp2[7] );
      LineInters2f ( &acp2[2], &acp2[11], &acp2[6], &acp2[13], &acp2[5] );
      LineInters2f ( &acp2[2], &acp2[11], &acp2[9], &acp2[14], &acp2[8] );
    }
    else {
      AddVector2Mf ( &domcp[0], &acp2[0], 2.0/(float)hole_k, &acp2[0] );
    }
    pkn_SubtractMatrixf ( 1, 2*npoints, 0, (float*)acp2, 0, (float*)domcp,
                          0, (float*)acp2 );
    pkn_AddMatrixf ( 1, 2*npoints, 0, (float*)acp1, 0, (float*)acp2,
                     0, (float*)acp1 );
  }

  pkn_AddMatrixf ( 1, 2*npoints, 0, (float*)domcp, 0, (float*)acp1,
                   0, (float*)domcp );

way_out:
  pkv_SetScratchMemTop ( sp );
  return npoints;
} /*InitGHDomainNet*/

/* ///////////////////////////////////////////////////////////////////////// */
static void SetAuxCP3 ( point3f *p0, point3f *p1, point3f *p2, point3f *p3,
                        point3f *p4, point3f *p5, point3f *p6,
                        float phi, int j, point3f *cp )
{
  vector3f vxi, veta, vzeta;
  float a[9], ai[9], x[3], r;
  int    P[2], Q[2];

  r = 1.0/sin(phi);
  SubtractPoints3f ( p1, p0, &vxi );
  NormalizeVector3f ( &vxi );
  SubtractPoints3f ( p3, p2, &veta );
  CrossProduct3f ( &veta, &vxi, &vzeta );
  MultVector3f ( 0.5, &vzeta, &vzeta );
  SubtractPoints3f ( p5, p4, &veta );
  a[0] = vxi.x;  a[1] = veta.x;  a[2] = vzeta.x;
  a[3] = vxi.y;  a[4] = veta.y;  a[5] = vzeta.y;
  a[6] = vxi.z;  a[7] = veta.z;  a[8] = vzeta.z;
  memcpy ( ai, a, 9*sizeof(float) );
  pkn_GaussDecomposePLUQf ( 3, ai, P, Q );
  memcpy ( x, p6, 3*sizeof(float) );
  pkn_multiSolvePLUQf ( 3, ai, P, Q, 1, 1, x );
  AddVector3Mf ( p4, &veta, r*sin((float)j/3.0*phi), cp );
  AddVector3Mf ( cp, &vzeta, r*(cos((float)j/3.0*phi)-cos(phi)), cp );
} /*SetAuxCP3*/

int InitGHSurfNetf ( int hole_k, const vector3f *v,
                     int nparams, const float *param,
                     point3f *surfcp )
{
  void     *sp;
  int      i, j, k, l, m, npoints;
  point3f  *acp1, *acp2, *acp3;
  vector3f *cv, av;
  float   phi;

  sp = pkv_GetScratchMemTop ();
  npoints = 1+12*hole_k;

  cv = (vector3f*)pkv_GetScratchMem ( hole_k*sizeof(vector3f) );
  if ( !cv ) {
    memset ( surfcp, 0, npoints*sizeof(point3f) );
    goto way_out;
  }
  memcpy ( cv, v, hole_k*sizeof(vector3f) );
        /* modify the construction vectors */
  for ( i = 0; i < hole_k-1; i += 2 ) {
    cv[i].z += param[3]*sqrt(cv[i].x*cv[i].x + cv[i].y*cv[i].y);
    cv[i+1].z -= param[3]*sqrt(cv[i+1].x*cv[i+1].x + cv[i+1].y*cv[i+1].y);
  }
        /* rhomboidal */
  SetPoint3f ( &surfcp[0], 0.0, 0.0, 0.0 );
  for ( i = 0, k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;
    for ( l = 1; l < 4; l++ ) {
      AddVector3Mf ( &surfcp[0], &cv[i], ((float)l)/3.0, &surfcp[k+l] );
      for ( m = 1; m < 4; m++ )
        AddVector3Mf ( &surfcp[k+l], &cv[j], ((float)m)/3.0, &surfcp[k+l+3*m] );
    }    
  }
  if ( !(acp1 = (point3f*)pkv_GetScratchMem ( 3*npoints*sizeof(point3f) )) )
    goto way_out;
  acp2 = &acp1[npoints];
  acp3 = &acp2[npoints];

        /* 1st modification */
  memcpy ( acp1, surfcp, npoints*sizeof(point3f) );
  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    MidPoint3f ( &acp1[12*i+12], &acp1[12*j+12], &acp1[12*j+3] );
    InterPoint3f ( &acp1[0], &acp1[12*j+3], 1.0/3.0, &acp1[12*j+1] );
    InterPoint3f ( &acp1[0], &acp1[12*j+3], 2.0/3.0, &acp1[12*j+2] );
    InterPoint3f ( &acp1[12*i+12], &acp1[12*j+3], 1.0/3.0, &acp1[12*i+11] );
    InterPoint3f ( &acp1[12*i+12], &acp1[12*j+3], 2.0/3.0, &acp1[12*i+10] );
    InterPoint3f ( &acp1[12*j+3], &acp1[12*j+12], 1.0/3.0, &acp1[12*j+6] );
    InterPoint3f ( &acp1[12*j+3], &acp1[12*j+12], 2.0/3.0, &acp1[12*j+9] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3f ( &acp1[k+2], &acp1[k+11], &acp1[m+2], &acp1[k+9], &acp1[k+8] );
    LineInters3f ( &acp1[k+2], &acp1[k+11], &acp1[m+1], &acp1[k+6], &acp1[k+5] );
    LineInters3f ( &acp1[k+1], &acp1[k+10], &acp1[m+2], &acp1[k+9], &acp1[k+7] );
    LineInters3f ( &acp1[k+1], &acp1[k+10], &acp1[m+1], &acp1[k+6], &acp1[k+4] );
  }
  pkn_SubtractMatrixf ( 1, 3*npoints, 0, (float*)acp1, 0, (float*)surfcp,
                        0, (float*)acp1 );
  pkn_MultMatrixNumf ( 1, 3*npoints, 0, (float*)acp1, param[0],
                       0, (float*)acp1 );
        /* 2nd modification */
  memcpy ( acp2, surfcp, npoints*sizeof(point3f) );
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    MidPoint3f ( &acp2[k+12], &acp2[m+12], &acp2[m+3] );
    InterPoint3f ( &acp2[k+12], &acp2[m+3], 1.0/3.0, &acp2[k+11] );
    InterPoint3f ( &acp2[k+12], &acp2[m+3], 2.0/3.0, &acp2[k+10] );
    InterPoint3f ( &acp2[m+3], &acp2[m+12], 1.0/3.0, &acp2[m+6] );
    InterPoint3f ( &acp2[m+3], &acp2[m+12], 2.0/3.0, &acp2[m+9] );
  }
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3f ( &acp2[0], &acp2[m+3], &acp2[k+6], &acp2[m+10], &acp2[m+1] );
    LineInters3f ( &acp2[0], &acp2[m+3], &acp2[k+9], &acp2[m+11], &acp2[m+2] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3f ( &acp2[k+2], &acp2[k+11], &acp2[m+2], &acp2[k+9], &acp2[k+8] );
    LineInters3f ( &acp2[k+2], &acp2[k+11], &acp2[m+1], &acp2[k+6], &acp2[k+5] );
    LineInters3f ( &acp2[k+1], &acp2[k+10], &acp2[m+2], &acp2[k+9], &acp2[k+7] );
    LineInters3f ( &acp2[k+1], &acp2[k+10], &acp2[m+1], &acp2[k+6], &acp2[k+4] );
  }
  pkn_SubtractMatrixf ( 1, 3*npoints, 0, (float*)acp2, 0, (float*)surfcp,
                        0, (float*)acp2 );
  pkn_AddMatrixMf ( 1, 3*npoints, 0, (float*)acp1, 0, (float*)acp2,
                    param[1], 0, (float*)acp1 );
  pkn_AddMatrixf ( 1, 3*npoints, 0, (float*)surfcp, 0, (float*)acp1,
                   0, (float*)acp3 );
        /* 3rd modification */
  phi = 0.5*param[4]*PI;
  if ( phi > 0.5*PI ) phi = 0.5*PI;
  if ( phi > 0.0 ) {
    memset ( acp1, 0, npoints*sizeof(point3f) );
    memcpy ( acp2, acp3, npoints*sizeof(point3f) );
    SetVector3f ( &acp2[0], 0.0, 0.0, 0.0 );
    for ( i = k = 0;  i < hole_k;  i += 2, k += 24 ) {
      j = (i+1) % hole_k;  m = 12*j;
      SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[0], phi, 0, &av );
      AddVector3f ( &acp2[0], &av, &acp2[0] );
      SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[k+1], phi, 1, &acp2[k+1] );
      SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[k+2], phi, 2, &acp2[k+2] );
      for ( l = 1; l < 4; l++ ) {
        SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l],
                   phi, 0, &acp2[m+l] );
        SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l+3],
                   phi, 1, &acp2[m+l+3] );
        SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[m+l+9], &acp3[m+l+6],
                   phi, 2, &acp2[m+l+6] );
        SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[k+3*l+3], &acp3[k+3*l+1],
                   phi, 1, &acp2[k+3*l+1] );
        SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                   &acp3[m+l], &acp3[k+3*l+3], &acp3[k+3*l+2],
                   phi, 2, &acp2[k+3*l+2] );
      }
    }
    if ( hole_k & 0x01 ) {
      MultVector3f ( 2.0/(float)(hole_k+1), &acp2[0], &av ) ;
      AddVector3f ( &surfcp[0], &av, &acp2[0] );
      AddVector3f ( &surfcp[12], &av, &acp2[12] );
      InterPoint3f ( &acp2[3], &acp2[12], 1.0/3.0, &acp2[6] );
      InterPoint3f ( &acp2[3], &acp2[12], 2.0/3.0, &acp2[9] );
      InterPoint3f ( &acp2[15], &acp2[12], 1.0/3.0, &acp2[10] );
      InterPoint3f ( &acp2[15], &acp2[12], 2.0/3.0, &acp2[11] );
      LineInters3f ( &acp2[1], &acp2[10], &acp2[6], &acp2[13], &acp2[4] );
      LineInters3f ( &acp2[1], &acp2[10], &acp2[9], &acp2[14], &acp2[7] );
      LineInters3f ( &acp2[2], &acp2[11], &acp2[6], &acp2[13], &acp2[5] );
      LineInters3f ( &acp2[2], &acp2[11], &acp2[9], &acp2[14], &acp2[8] );
    }
    else {
      AddVector3Mf ( &surfcp[0], &acp2[0], 2.0/(float)hole_k, &acp2[0] );
    }
    pkn_SubtractMatrixf ( 1, 3*npoints, 0, (float*)acp2, 0, (float*)surfcp,
                          0, (float*)acp2 );
    pkn_AddMatrixf ( 1, 3*npoints, 0, (float*)acp1, 0, (float*)acp2,
                     0, (float*)acp1 );
  }
        /* 4th modification */
  memcpy ( acp2, surfcp, npoints*sizeof(point3f) );
  for ( i = 0; i < npoints; i++ ) {
    acp2[i].z -= sqrt(acp2[i].x*acp2[i].x + acp2[i].y*acp2[i].y);
  }
  pkn_SubtractMatrixf ( 1, 3*npoints, 0, (float*)acp2, 0, (float*)surfcp,
                        0, (float*)acp2 );
  pkn_AddMatrixMf ( 1, 3*npoints, 0, (float*)acp1, 0, (float*)acp2,
                    param[2], 0, (float*)acp1 );
  pkn_AddMatrixf ( 1, 3*npoints, 0, (float*)surfcp, 0, (float*)acp1,
                   0, (float*)surfcp );

way_out:
  pkv_SetScratchMemTop ( sp );
  return npoints;
} /*InitGHSurfNet*/

