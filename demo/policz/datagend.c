
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
#include "egholed.h"

#include "datagend.h"

void InitGHKnotsd ( int hole_k, double *knots )
{
  int i;

  knots[0] = knots[1] = 0.0;
  for ( i = 2; i < 9; i++ )
    knots[i] = (i-1)/8.0;
  knots[9] = knots[10] = 1.0;
  for ( i = 1; i < hole_k; i++ )
    memcpy ( &knots[11*i], &knots[0], 11*sizeof(double) );
} /*InitGHKnotsd*/

void InitGHVectors2d ( int hole_k, vector2d *v )
{
  int i;
  trans2d tr;
  double  a;

  a = 2.0*PI/(double)hole_k;
  IdentTrans2d ( &tr );
  RotTrans2d ( &tr, a );
  SetVector2d ( &v[0], 3.0*cos(a), -3.0*sin(a) );
  for ( i = 1; i < hole_k; i++ )
    TransVector2d ( &tr, &v[i-1], &v[i] );
} /*InitGHVectors2d*/

void InitGHVectors3d ( int hole_k, vector3d *v )
{
  int     i;
  trans3d tr;
  double  a;

  a = 2.0*PI/(double)hole_k;
  IdentTrans3d ( &tr );
  RotZTrans3d ( &tr, a );
  SetVector3d ( &v[0], 3.0*cos(a), -3.0*sin(a), 0.0 );
  for ( i = 1; i < hole_k; i++ )
    TransVector3d ( &tr, &v[i-1], &v[i] );
} /*InitGHVectors3d*/

void LineInters2d ( const point2d *p0, const point2d *p1,
                    const point2d *q0, const point2d *q1,
                    point2d *p )
{
  double a[6];

  a[0] = p1->x-p0->x;  a[1] = q0->x-q1->x;
  a[2] = p1->y-p0->y;  a[3] = q0->y-q1->y;
  a[4] = q0->x-p0->x;  a[5] = q0->y-p0->y;
  pkn_multiGaussSolveLinEqd ( 2, a, 1, 1, &a[4] );
  InterPoint2d ( p0, p1, a[4], p );
} /*LineInters2d*/

void LineInters3d ( const point3d *p0, const point3d *p1,
                    const point3d *q0, const point3d *q1,
                    point3d *p )
{
#define TOL 1.0e-8
  double   a[12];
  vector3d u, v, w;

  SubtractPoints3d ( p1, p0, &u );
  SubtractPoints3d ( q1, q0, &v );
  CrossProduct3d ( &u, &v, &w );
  if ( fabs(w.x)+fabs(w.y)+fabs(w.z) >= TOL ) {
    a[0] = u.x;  a[1] = -v.x;  a[2] = -w.x;
    a[3] = u.y;  a[4] = -v.y;  a[5] = -w.y;
    a[6] = u.z;  a[7] = -v.z;  a[8] = -w.z;
    a[9] = q0->x-p0->x;  a[10] = q0->y-p0->y;  a[11] = q0->z-p0->z;
    pkn_multiGaussSolveLinEqd ( 3, a, 1, 1, &a[9] );
    InterPoint3d ( p0, p1, a[9], &u );
    InterPoint3d ( q0, q1, a[10], &v );
  }
  else {
    MidPoint3d ( p0, p1, &u );
    MidPoint3d ( q0, q1, &v );
  }
  MidPoint3d ( &u, &v, p );
#undef TOL
} /*LineInters3d*/

static void SetAuxCP2 ( point2d *p0, point2d *p1, point2d *p2, point2d *p3,
                        point2d *p4, point2d *p5, point2d *p6,
                        double phi, int j, point2d *cp )
{
  vector3d vxi, veta, vzeta;
  double a[9], ai[9], x[3], r;
  int    P[2], Q[2];

  r = 1.0/sin(phi);
  SetVector3d ( &vxi, p1->x-p0->x, p1->y-p0->y, 0.0 );
  NormalizeVector3d ( &vxi );
  SetVector3d ( &veta, p3->x-p2->x, p3->y-p2->y, 0.0 );
  CrossProduct3d ( &veta, &vxi, &vzeta );
  MultVector3d ( 0.5, &vzeta, &vzeta );
  SetVector3d ( &veta, p5->x-p4->x, p5->y-p4->y, 0.0 );
  a[0] = vxi.x;  a[1] = veta.x;  a[2] = vzeta.x;
  a[3] = vxi.y;  a[4] = veta.y;  a[5] = vzeta.y;
  a[6] = vxi.z;  a[7] = veta.z;  a[8] = vzeta.z;
  memcpy ( ai, a, 9*sizeof(double) );
  pkn_GaussDecomposePLUQd ( 3, ai, P, Q );
  memcpy ( x, p6, 3*sizeof(double) );
  pkn_multiSolvePLUQd ( 3, ai, P, Q, 1, 1, x );
  AddVector2Md ( p4, (vector2d*)(void*)&veta, r*sin((double)j/3.0*phi), cp );
  AddVector2Md ( cp, (vector2d*)(void*)&vzeta, r*(cos((double)j/3.0*phi)-cos(phi)), cp );
} /*SetAuxCP2*/

int InitGHDomainNetd ( int hole_k, const vector2d *v,
                       int nparams, const double *param,
                       point2d *domcp )
{
  void     *sp;
  int      i, j, k, l, m, npoints;
  point2d  *acp1, *acp2, *acp3;
  vector2d av;
  double   phi;

  sp = pkv_GetScratchMemTop ();
  npoints = 1+12*hole_k;

        /* rhomboidal */
  SetPoint2d ( &domcp[0], 0.0, 0.0 );
  for ( i = 0, k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;
    for ( l = 1; l < 4; l++ ) {
      AddVector2Md ( &domcp[0], &v[i], ((double)l)/3.0, &domcp[k+l] );
      for ( m = 1; m < 4; m++ )
        AddVector2Md ( &domcp[k+l], &v[j], ((double)m)/3.0, &domcp[k+l+3*m] );
    }    
  }
  if ( !(acp1 = (point2d*)pkv_GetScratchMem ( 3*npoints*sizeof(point2d) )) )
    goto way_out;
  acp2 = &acp1[npoints];
  acp3 = &acp2[npoints];

        /* 1st modification */
  memcpy ( acp1, domcp, npoints*sizeof(point2d) );
  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    MidPoint2d ( &acp1[12*i+12], &acp1[12*j+12], &acp1[12*j+3] );
    InterPoint2d ( &acp1[0], &acp1[12*j+3], 1.0/3.0, &acp1[12*j+1] );
    InterPoint2d ( &acp1[0], &acp1[12*j+3], 2.0/3.0, &acp1[12*j+2] );
    InterPoint2d ( &acp1[12*i+12], &acp1[12*j+3], 1.0/3.0, &acp1[12*i+11] );
    InterPoint2d ( &acp1[12*i+12], &acp1[12*j+3], 2.0/3.0, &acp1[12*i+10] );
    InterPoint2d ( &acp1[12*j+3], &acp1[12*j+12], 1.0/3.0, &acp1[12*j+6] );
    InterPoint2d ( &acp1[12*j+3], &acp1[12*j+12], 2.0/3.0, &acp1[12*j+9] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2d ( &acp1[k+2], &acp1[k+11], &acp1[m+2], &acp1[k+9], &acp1[k+8] );
    LineInters2d ( &acp1[k+2], &acp1[k+11], &acp1[m+1], &acp1[k+6], &acp1[k+5] );
    LineInters2d ( &acp1[k+1], &acp1[k+10], &acp1[m+2], &acp1[k+9], &acp1[k+7] );
    LineInters2d ( &acp1[k+1], &acp1[k+10], &acp1[m+1], &acp1[k+6], &acp1[k+4] );
  }
  pkn_SubtractMatrixd ( 1, 2*npoints, 0, (double*)acp1, 0, (double*)domcp,
                        0, (double*)acp1 );
  pkn_MultMatrixNumd ( 1, 2*npoints, 0, (double*)acp1, param[0],
                       0, (double*)acp1 );
        /* 2nd modification */
  memcpy ( acp2, domcp, npoints*sizeof(point2d) );
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    MidPoint2d ( &acp2[k+12], &acp2[m+12], &acp2[m+3] );
    InterPoint2d ( &acp2[k+12], &acp2[m+3], 1.0/3.0, &acp2[k+11] );
    InterPoint2d ( &acp2[k+12], &acp2[m+3], 2.0/3.0, &acp2[k+10] );
    InterPoint2d ( &acp2[m+3], &acp2[m+12], 1.0/3.0, &acp2[m+6] );
    InterPoint2d ( &acp2[m+3], &acp2[m+12], 2.0/3.0, &acp2[m+9] );
  }
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2d ( &acp2[0], &acp2[m+3], &acp2[k+6], &acp2[m+10], &acp2[m+1] );
    LineInters2d ( &acp2[0], &acp2[m+3], &acp2[k+9], &acp2[m+11], &acp2[m+2] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters2d ( &acp2[k+2], &acp2[k+11], &acp2[m+2], &acp2[k+9], &acp2[k+8] );
    LineInters2d ( &acp2[k+2], &acp2[k+11], &acp2[m+1], &acp2[k+6], &acp2[k+5] );
    LineInters2d ( &acp2[k+1], &acp2[k+10], &acp2[m+2], &acp2[k+9], &acp2[k+7] );
    LineInters2d ( &acp2[k+1], &acp2[k+10], &acp2[m+1], &acp2[k+6], &acp2[k+4] );
  }
  pkn_SubtractMatrixd ( 1, 2*npoints, 0, (double*)acp2, 0, (double*)domcp,
                        0, (double*)acp2 );
  pkn_AddMatrixMd ( 1, 2*npoints, 0, (double*)acp1, 0, (double*)acp2,
                    param[1], 0, (double*)acp1 );
  pkn_AddMatrixd ( 1, 2*npoints, 0, (double*)domcp, 0, (double*)acp1,
                   0, (double*)acp3 );

        /* 3rd modification */
  phi = 0.5*param[2]*PI;
  if ( phi > 0.5*PI ) phi = 0.5*PI;
  if ( phi > 0.0 ) {
    memset ( acp1, 0, npoints*sizeof(point2d) );
    memcpy ( acp2, acp3, npoints*sizeof(point2d) );
    SetVector2d ( &acp2[0], 0.0, 0.0 );
    for ( i = k = 0;  i < hole_k;  i += 2, k += 24 ) {
      j = (i+1) % hole_k;  m = 12*j;
      SetAuxCP2 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[0], phi, 0, &av );
      AddVector2d ( &acp2[0], &av, &acp2[0] );
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
      MultVector2d ( 2.0/(double)(hole_k+1), &acp2[0], &av ) ;
      AddVector2d ( &domcp[0], &av, &acp2[0] );
      AddVector2d ( &domcp[12], &av, &acp2[12] );
      InterPoint2d ( &acp2[3], &acp2[12], 1.0/3.0, &acp2[6] );
      InterPoint2d ( &acp2[3], &acp2[12], 2.0/3.0, &acp2[9] );
      InterPoint2d ( &acp2[15], &acp2[12], 1.0/3.0, &acp2[10] );
      InterPoint2d ( &acp2[15], &acp2[12], 2.0/3.0, &acp2[11] );
      LineInters2d ( &acp2[1], &acp2[10], &acp2[6], &acp2[13], &acp2[4] );
      LineInters2d ( &acp2[1], &acp2[10], &acp2[9], &acp2[14], &acp2[7] );
      LineInters2d ( &acp2[2], &acp2[11], &acp2[6], &acp2[13], &acp2[5] );
      LineInters2d ( &acp2[2], &acp2[11], &acp2[9], &acp2[14], &acp2[8] );
    }
    else {
      AddVector2Md ( &domcp[0], &acp2[0], 2.0/(double)hole_k, &acp2[0] );
    }
    pkn_SubtractMatrixd ( 1, 2*npoints, 0, (double*)acp2, 0, (double*)domcp,
                          0, (double*)acp2 );
    pkn_AddMatrixd ( 1, 2*npoints, 0, (double*)acp1, 0, (double*)acp2,
                     0, (double*)acp1 );
  }

  pkn_AddMatrixd ( 1, 2*npoints, 0, (double*)domcp, 0, (double*)acp1,
                   0, (double*)domcp );

        /* 4th modification */
  switch ( hole_k ) {
case 3:  acp1 = &egh_eigendom3d[0];   break;
case 5:  acp1 = &egh_eigendom5d[0];   break;
case 6:  acp1 = &egh_eigendom6d[0];   break;
case 7:  acp1 = &egh_eigendom7d[0];   break;
case 8:  acp1 = &egh_eigendom8d[0];   break;
case 9:  acp1 = &egh_eigendom9d[0];   break;
case 10: acp1 = &egh_eigendom10d[0];  break;
case 11: acp1 = &egh_eigendom11d[0];  break;
case 12: acp1 = &egh_eigendom12d[0];  break;
case 13: acp1 = &egh_eigendom13d[0];  break;
case 14: acp1 = &egh_eigendom14d[0];  break;
case 15: acp1 = &egh_eigendom15d[0];  break;
case 16: acp1 = &egh_eigendom16d[0];  break;
default: goto way_out;
  }
  pkn_MatrixLinCombd ( 1, 2*(12*hole_k+1), 0, (double*)domcp, 1.0-param[3],
                       0, (double*)acp1, param[3], 0, (double*)domcp );

way_out:
  pkv_SetScratchMemTop ( sp );
  return npoints;
} /*InitGHDomainNet*/

/* ///////////////////////////////////////////////////////////////////////// */
static void SetAuxCP3 ( point3d *p0, point3d *p1, point3d *p2, point3d *p3,
                        point3d *p4, point3d *p5, point3d *p6,
                        double phi, int j, point3d *cp )
{
  vector3d vxi, veta, vzeta;
  double a[9], ai[9], x[3], r;
  int    P[2], Q[2];

  r = 1.0/sin(phi);
  SubtractPoints3d ( p1, p0, &vxi );
  NormalizeVector3d ( &vxi );
  SubtractPoints3d ( p3, p2, &veta );
  CrossProduct3d ( &veta, &vxi, &vzeta );
  MultVector3d ( 0.5, &vzeta, &vzeta );
  SubtractPoints3d ( p5, p4, &veta );
  a[0] = vxi.x;  a[1] = veta.x;  a[2] = vzeta.x;
  a[3] = vxi.y;  a[4] = veta.y;  a[5] = vzeta.y;
  a[6] = vxi.z;  a[7] = veta.z;  a[8] = vzeta.z;
  memcpy ( ai, a, 9*sizeof(double) );
  pkn_GaussDecomposePLUQd ( 3, ai, P, Q );
  memcpy ( x, p6, 3*sizeof(double) );
  pkn_multiSolvePLUQd ( 3, ai, P, Q, 1, 1, x );
  AddVector3Md ( p4, &veta, r*sin((double)j/3.0*phi), cp );
  AddVector3Md ( cp, &vzeta, r*(cos((double)j/3.0*phi)-cos(phi)), cp );
} /*SetAuxCP3*/

int InitGHSurfNetd ( int hole_k, const vector3d *v,
                     int nparams, const double *param,
                     point3d *surfcp )
{
  void     *sp;
  int      i, j, k, l, m, npoints;
  point3d  *acp1, *acp2, *acp3;
  vector3d *cv, av;
  double   phi;

  sp = pkv_GetScratchMemTop ();
  npoints = 1+12*hole_k;

  cv = (vector3d*)pkv_GetScratchMem ( hole_k*sizeof(vector3d) );
  if ( !cv ) {
    memset ( surfcp, 0, npoints*sizeof(point3d) );
    goto way_out;
  }
  memcpy ( cv, v, hole_k*sizeof(vector3d) );
        /* modify the construction vectors */
  for ( i = 0; i < hole_k-1; i += 2 ) {
    cv[i].z += param[3]*sqrt(cv[i].x*cv[i].x + cv[i].y*cv[i].y);
    cv[i+1].z -= param[3]*sqrt(cv[i+1].x*cv[i+1].x + cv[i+1].y*cv[i+1].y);
  }
        /* rhomboidal */
  SetPoint3d ( &surfcp[0], 0.0, 0.0, 0.0 );
  for ( i = 0, k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;
    for ( l = 1; l < 4; l++ ) {
      AddVector3Md ( &surfcp[0], &cv[i], ((double)l)/3.0, &surfcp[k+l] );
      for ( m = 1; m < 4; m++ )
        AddVector3Md ( &surfcp[k+l], &cv[j], ((double)m)/3.0, &surfcp[k+l+3*m] );
    }    
  }
  if ( !(acp1 = (point3d*)pkv_GetScratchMem ( 3*npoints*sizeof(point3d) )) )
    goto way_out;
  acp2 = &acp1[npoints];
  acp3 = &acp2[npoints];

        /* 1st modification */
  memcpy ( acp1, surfcp, npoints*sizeof(point3d) );
  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    MidPoint3d ( &acp1[12*i+12], &acp1[12*j+12], &acp1[12*j+3] );
    InterPoint3d ( &acp1[0], &acp1[12*j+3], 1.0/3.0, &acp1[12*j+1] );
    InterPoint3d ( &acp1[0], &acp1[12*j+3], 2.0/3.0, &acp1[12*j+2] );
    InterPoint3d ( &acp1[12*i+12], &acp1[12*j+3], 1.0/3.0, &acp1[12*i+11] );
    InterPoint3d ( &acp1[12*i+12], &acp1[12*j+3], 2.0/3.0, &acp1[12*i+10] );
    InterPoint3d ( &acp1[12*j+3], &acp1[12*j+12], 1.0/3.0, &acp1[12*j+6] );
    InterPoint3d ( &acp1[12*j+3], &acp1[12*j+12], 2.0/3.0, &acp1[12*j+9] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3d ( &acp1[k+2], &acp1[k+11], &acp1[m+2], &acp1[k+9], &acp1[k+8] );
    LineInters3d ( &acp1[k+2], &acp1[k+11], &acp1[m+1], &acp1[k+6], &acp1[k+5] );
    LineInters3d ( &acp1[k+1], &acp1[k+10], &acp1[m+2], &acp1[k+9], &acp1[k+7] );
    LineInters3d ( &acp1[k+1], &acp1[k+10], &acp1[m+1], &acp1[k+6], &acp1[k+4] );
  }
  pkn_SubtractMatrixd ( 1, 3*npoints, 0, (double*)acp1, 0, (double*)surfcp,
                        0, (double*)acp1 );
  pkn_MultMatrixNumd ( 1, 3*npoints, 0, (double*)acp1, param[0],
                       0, (double*)acp1 );
        /* 2nd modification */
  memcpy ( acp2, surfcp, npoints*sizeof(point3d) );
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    MidPoint3d ( &acp2[k+12], &acp2[m+12], &acp2[m+3] );
    InterPoint3d ( &acp2[k+12], &acp2[m+3], 1.0/3.0, &acp2[k+11] );
    InterPoint3d ( &acp2[k+12], &acp2[m+3], 2.0/3.0, &acp2[k+10] );
    InterPoint3d ( &acp2[m+3], &acp2[m+12], 1.0/3.0, &acp2[m+6] );
    InterPoint3d ( &acp2[m+3], &acp2[m+12], 2.0/3.0, &acp2[m+9] );
  }
  for ( i = k = 0;  i < hole_k;  i+= 2, k += 24 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3d ( &acp2[0], &acp2[m+3], &acp2[k+6], &acp2[m+10], &acp2[m+1] );
    LineInters3d ( &acp2[0], &acp2[m+3], &acp2[k+9], &acp2[m+11], &acp2[m+2] );
  }
  for ( i = k = 0;  i < hole_k;  i++, k += 12 ) {
    j = (i+1) % hole_k;  m = 12*j;
    LineInters3d ( &acp2[k+2], &acp2[k+11], &acp2[m+2], &acp2[k+9], &acp2[k+8] );
    LineInters3d ( &acp2[k+2], &acp2[k+11], &acp2[m+1], &acp2[k+6], &acp2[k+5] );
    LineInters3d ( &acp2[k+1], &acp2[k+10], &acp2[m+2], &acp2[k+9], &acp2[k+7] );
    LineInters3d ( &acp2[k+1], &acp2[k+10], &acp2[m+1], &acp2[k+6], &acp2[k+4] );
  }
  pkn_SubtractMatrixd ( 1, 3*npoints, 0, (double*)acp2, 0, (double*)surfcp,
                        0, (double*)acp2 );
  pkn_AddMatrixMd ( 1, 3*npoints, 0, (double*)acp1, 0, (double*)acp2,
                    param[1], 0, (double*)acp1 );
  pkn_AddMatrixd ( 1, 3*npoints, 0, (double*)surfcp, 0, (double*)acp1,
                   0, (double*)acp3 );
        /* 3rd modification */
  phi = 0.5*param[4]*PI;
  if ( phi > 0.5*PI ) phi = 0.5*PI;
  if ( phi > 0.0 ) {
    memset ( acp1, 0, npoints*sizeof(point3d) );
    memcpy ( acp2, acp3, npoints*sizeof(point3d) );
    SetVector3d ( &acp2[0], 0.0, 0.0, 0.0 );
    for ( i = k = 0;  i < hole_k;  i += 2, k += 24 ) {
      j = (i+1) % hole_k;  m = 12*j;
      SetAuxCP3 ( &acp3[0], &acp3[m+3], &acp3[m+12], &acp3[k+12],
                 &acp3[0], &acp3[k+3], &acp3[0], phi, 0, &av );
      AddVector3d ( &acp2[0], &av, &acp2[0] );
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
      MultVector3d ( 2.0/(double)(hole_k+1), &acp2[0], &av ) ;
      AddVector3d ( &surfcp[0], &av, &acp2[0] );
      AddVector3d ( &surfcp[12], &av, &acp2[12] );
      InterPoint3d ( &acp2[3], &acp2[12], 1.0/3.0, &acp2[6] );
      InterPoint3d ( &acp2[3], &acp2[12], 2.0/3.0, &acp2[9] );
      InterPoint3d ( &acp2[15], &acp2[12], 1.0/3.0, &acp2[10] );
      InterPoint3d ( &acp2[15], &acp2[12], 2.0/3.0, &acp2[11] );
      LineInters3d ( &acp2[1], &acp2[10], &acp2[6], &acp2[13], &acp2[4] );
      LineInters3d ( &acp2[1], &acp2[10], &acp2[9], &acp2[14], &acp2[7] );
      LineInters3d ( &acp2[2], &acp2[11], &acp2[6], &acp2[13], &acp2[5] );
      LineInters3d ( &acp2[2], &acp2[11], &acp2[9], &acp2[14], &acp2[8] );
    }
    else {
      AddVector3Md ( &surfcp[0], &acp2[0], 2.0/(double)hole_k, &acp2[0] );
    }
    pkn_SubtractMatrixd ( 1, 3*npoints, 0, (double*)acp2, 0, (double*)surfcp,
                          0, (double*)acp2 );
    pkn_AddMatrixd ( 1, 3*npoints, 0, (double*)acp1, 0, (double*)acp2,
                     0, (double*)acp1 );
  }
        /* 4th modification */
  memcpy ( acp2, surfcp, npoints*sizeof(point3d) );
  for ( i = 0; i < npoints; i++ ) {
    acp2[i].z -= sqrt(acp2[i].x*acp2[i].x + acp2[i].y*acp2[i].y);
  }
  pkn_SubtractMatrixd ( 1, 3*npoints, 0, (double*)acp2, 0, (double*)surfcp,
                        0, (double*)acp2 );
  pkn_AddMatrixMd ( 1, 3*npoints, 0, (double*)acp1, 0, (double*)acp2,
                    param[2], 0, (double*)acp1 );
  pkn_AddMatrixd ( 1, 3*npoints, 0, (double*)surfcp, 0, (double*)acp1,
                   0, (double*)surfcp );

way_out:
  pkv_SetScratchMemTop ( sp );
  return npoints;
} /*InitGHSurfNet*/

