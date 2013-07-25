
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "glutki.h"


void glutSolidCube ( double a )
{
  a *= 0.5;
  glBegin ( GL_POLYGON );
  glNormal3d ( 0.0, 0.0, -1.0 );
  glVertex3d ( -a, -a, -a );
  glVertex3d ( a, -a, -a );
  glVertex3d ( a, a, -a );
  glVertex3d ( -a, a, -a );
  glEnd ();
  glBegin ( GL_POLYGON );
  glNormal3d ( 0.0, 0.0, 1.0 );
  glVertex3d ( -a, -a, a );
  glVertex3d ( a, -a, a );
  glVertex3d ( a, a, a );
  glVertex3d ( -a, a, a );
  glEnd ();
  glBegin ( GL_POLYGON );
  glNormal3d ( 0.0, -1.0, 0.0 );
  glVertex3d ( -a, -a, -a );
  glVertex3d ( a, -a, -a );
  glVertex3d ( a, -a, a );
  glVertex3d ( -a, -a, a );
  glEnd ();
  glBegin ( GL_POLYGON );
  glNormal3d ( 0.0, 1.0, 0.0 );
  glVertex3d ( -a, a, -a );
  glVertex3d ( a, a, -a );
  glVertex3d ( a, a, a );
  glVertex3d ( -a, a, a );
  glEnd ();
  glBegin ( GL_POLYGON );
  glNormal3d ( -1.0, 0.0, 0.0 );
  glVertex3d ( -a, -a, -a );
  glVertex3d ( -a, a, -a );
  glVertex3d ( -a, a, a );
  glVertex3d ( -a, -a, a );
  glEnd ();
  glBegin ( GL_POLYGON );
  glNormal3d ( 1.0, 0.0, 0.0 );
  glVertex3d ( a, -a, -a );
  glVertex3d ( a, a, -a );
  glVertex3d ( a, a, a );
  glVertex3d ( a, -a, a );
  glEnd ();
} /*glutSolidCube*/

void glutSolidSphere ( double r, int a, int b )
{
  int      i, j;
  double   alpha, beta;
  point3d  p, q, s;
  vector3d np, nq, ns;
  trans3d  tr1, tr2;

  r = fabs ( r );
  alpha = PI/a;
  beta = 2.0*PI/b;
  SetVector3d ( &np, 0.0, 0.0, 1.0 );
  IdentTrans3d ( &tr1 );
  RotXTrans3d ( &tr1, alpha );
  IdentTrans3d ( &tr2 );
  RotZTrans3d ( &tr2, beta );
  glBegin ( GL_TRIANGLE_FAN );
    MultVector3d ( r, &np, &p );
    glNormal3dv ( &np.x );
    glVertex3dv ( &p.x );
    TransVector3d ( &tr1, &np, &nq );
    ns = nq;
    for ( j = 0; j < b; j++ ) {
      glNormal3dv ( &nq.x );
      MultVector3d ( r, &nq, &q );
      glVertex3dv ( &q.x );
      TransVector3d ( &tr2, &nq, &nq );
    }
    glNormal3dv ( &ns.x );
    MultVector3d ( r, &ns, &s );
    glVertex3dv ( &s.x );
  glEnd ();

  for ( i = 1; i < a-1; i++ ) {
    np = ns;
    TransVector3d ( &tr1, &ns, &nq );
    ns = nq;
    glBegin ( GL_QUAD_STRIP );
      for ( j = 0; j <= b; j++ ) {
        glNormal3dv ( &np.x );
        MultVector3d ( r, &np, &p );
        glVertex3dv ( &p.x );
        glNormal3dv ( &nq.x );
        MultVector3d ( r, &nq, &q );
        glVertex3dv ( &q.x );
        TransVector3d ( &tr2, &np, &np );
        TransVector3d ( &tr2, &nq, &nq );
      }
    glEnd ();
  }

  IdentTrans3d ( &tr2 );
  RotZTrans3d ( &tr2, -beta );
  SetVector3d ( &np, 0.0, 0.0, -1.0 );
  glBegin ( GL_TRIANGLE_FAN );
    MultVector3d ( r, &np, &p );
    glNormal3dv ( &np.x );
    glVertex3dv ( &p.x );
    nq = ns;
    for ( j = 0; j < b; j++ ) {
      glNormal3dv ( &nq.x );
      MultVector3d ( r, &nq, &q );
      glVertex3dv ( &q.x );
      TransVector3d ( &tr2, &nq, &nq );
    }
    glNormal3dv ( &ns.x );
    MultVector3d ( r, &ns, &s );
    glVertex3dv ( &s.x );
  glEnd ();
} /*glutSolidSphere*/

void glutSolidTeapot ( double r )
{
  static const point3d TeapotPoints[306] =
   {{ 1.40000,  0.00000,  2.40000}, { 1.40000, -0.78400,  2.40000},
    { 0.78400, -1.40000,  2.40000}, { 0.00000, -1.40000,  2.40000},
    { 1.33750,  0.00000,  2.53125}, { 1.33750, -0.74900,  2.53125},
    { 0.74900, -1.33750,  2.53125}, { 0.00000, -1.33750,  2.53125},
    { 1.43750,  0.00000,  2.53125}, { 1.43750, -0.80500,  2.53125},
    { 0.80500, -1.43750,  2.53125}, { 0.00000, -1.43750,  2.53125},
    { 1.50000,  0.00000,  2.40000}, { 1.50000, -0.84000,  2.40000},
    { 0.84000, -1.50000,  2.40000}, { 0.00000, -1.50000,  2.40000},
    {-0.78400, -1.40000,  2.40000}, {-1.40000, -0.78400,  2.40000},
    {-1.40000,  0.00000,  2.40000}, {-0.74900, -1.33750,  2.53125},
    {-1.33750, -0.74900,  2.53125}, {-1.33750,  0.00000,  2.53125},
    {-0.80500, -1.43750,  2.53125}, {-1.43750, -0.80500,  2.53125},
    {-1.43750,  0.00000,  2.53125}, {-0.84000, -1.50000,  2.40000},
    {-1.50000, -0.84000,  2.40000}, {-1.50000,  0.00000,  2.40000},
    {-1.40000,  0.78400,  2.40000}, {-0.78400,  1.40000,  2.40000},
    { 0.00000,  1.40000,  2.40000}, {-1.33750,  0.74900,  2.53125},
    {-0.74900,  1.33750,  2.53125}, { 0.00000,  1.33750,  2.53125},
    {-1.43750,  0.80500,  2.53125}, {-0.80500,  1.43750,  2.53125},
    { 0.00000,  1.43750,  2.53125}, {-1.50000,  0.84000,  2.40000},
    {-0.84000,  1.50000,  2.40000}, { 0.00000,  1.50000,  2.40000},
    { 0.78400,  1.40000,  2.40000}, { 1.40000,  0.78400,  2.40000},
    { 0.74900,  1.33750,  2.53125}, { 1.33750,  0.74900,  2.53125},
    { 0.80500,  1.43750,  2.53125}, { 1.43750,  0.80500,  2.53125},
    { 0.84000,  1.50000,  2.40000}, { 1.50000,  0.84000,  2.40000},
    { 1.75000,  0.00000,  1.87500}, { 1.75000, -0.98000,  1.87500},
    { 0.98000, -1.75000,  1.87500}, { 0.00000, -1.75000,  1.87500},
    { 2.00000,  0.00000,  1.35000}, { 2.00000, -1.12000,  1.35000},
    { 1.12000, -2.00000,  1.35000}, { 0.00000, -2.00000,  1.35000},
    { 2.00000,  0.00000,  0.90000}, { 2.00000, -1.12000,  0.90000},
    { 1.12000, -2.00000,  0.90000}, { 0.00000, -2.00000,  0.90000},
    {-0.98000, -1.75000,  1.87500}, {-1.75000, -0.98000,  1.87500},
    {-1.75000,  0.00000,  1.87500}, {-1.12000, -2.00000,  1.35000},
    {-2.00000, -1.12000,  1.35000}, {-2.00000,  0.00000,  1.35000},
    {-1.12000, -2.00000,  0.90000}, {-2.00000, -1.12000,  0.90000},
    {-2.00000,  0.00000,  0.90000}, {-1.75000,  0.98000,  1.87500},
    {-0.98000,  1.75000,  1.87500}, { 0.00000,  1.75000,  1.87500},
    {-2.00000,  1.12000,  1.35000}, {-1.12000,  2.00000,  1.35000},
    { 0.00000,  2.00000,  1.35000}, {-2.00000,  1.12000,  0.90000},
    {-1.12000,  2.00000,  0.90000}, { 0.00000,  2.00000,  0.90000},
    { 0.98000,  1.75000,  1.87500}, { 1.75000,  0.98000,  1.87500},
    { 1.12000,  2.00000,  1.35000}, { 2.00000,  1.12000,  1.35000},
    { 1.12000,  2.00000,  0.90000}, { 2.00000,  1.12000,  0.90000},
    { 2.00000,  0.00000,  0.45000}, { 2.00000, -1.12000,  0.45000},
    { 1.12000, -2.00000,  0.45000}, { 0.00000, -2.00000,  0.45000},
    { 1.50000,  0.00000,  0.22500}, { 1.50000, -0.84000,  0.22500},
    { 0.84000, -1.50000,  0.22500}, { 0.00000, -1.50000,  0.22500},
    { 1.50000,  0.00000,  0.15000}, { 1.50000, -0.84000,  0.15000},
    { 0.84000, -1.50000,  0.15000}, { 0.00000, -1.50000,  0.15000},
    {-1.12000, -2.00000,  0.45000}, {-2.00000, -1.12000,  0.45000},
    {-2.00000,  0.00000,  0.45000}, {-0.84000, -1.50000,  0.22500},
    {-1.50000, -0.84000,  0.22500}, {-1.50000,  0.00000,  0.22500},
    {-0.84000, -1.50000,  0.15000}, {-1.50000, -0.84000,  0.15000},
    {-1.50000,  0.00000,  0.15000}, {-2.00000,  1.12000,  0.45000},
    {-1.12000,  2.00000,  0.45000}, { 0.00000,  2.00000,  0.45000},
    {-1.50000,  0.84000,  0.22500}, {-0.84000,  1.50000,  0.22500},
    { 0.00000,  1.50000,  0.22500}, {-1.50000,  0.84000,  0.15000},
    {-0.84000,  1.50000,  0.15000}, { 0.00000,  1.50000,  0.15000},
    { 1.12000,  2.00000,  0.45000}, { 2.00000,  1.12000,  0.45000},
    { 0.84000,  1.50000,  0.22500}, { 1.50000,  0.84000,  0.22500},
    { 0.84000,  1.50000,  0.15000}, { 1.50000,  0.84000,  0.15000},
    {-1.60000,  0.00000,  2.02500}, {-1.60000, -0.30000,  2.02500},
    {-1.50000, -0.30000,  2.25000}, {-1.50000,  0.00000,  2.25000},
    {-2.30000,  0.00000,  2.02500}, {-2.30000, -0.30000,  2.02500},
    {-2.50000, -0.30000,  2.25000}, {-2.50000,  0.00000,  2.25000},
    {-2.70000,  0.00000,  2.02500}, {-2.70000, -0.30000,  2.02500},
    {-3.00000, -0.30000,  2.25000}, {-3.00000,  0.00000,  2.25000},
    {-2.70000,  0.00000,  1.80000}, {-2.70000, -0.30000,  1.80000},
    {-3.00000, -0.30000,  1.80000}, {-3.00000,  0.00000,  1.80000},
    {-1.50000,  0.30000,  2.25000}, {-1.60000,  0.30000,  2.02500},
    {-2.50000,  0.30000,  2.25000}, {-2.30000,  0.30000,  2.02500},
    {-3.00000,  0.30000,  2.25000}, {-2.70000,  0.30000,  2.02500},
    {-3.00000,  0.30000,  1.80000}, {-2.70000,  0.30000,  1.80000},
    {-2.70000,  0.00000,  1.57500}, {-2.70000, -0.30000,  1.57500},
    {-3.00000, -0.30000,  1.35000}, {-3.00000,  0.00000,  1.35000},
    {-2.50000,  0.00000,  1.12500}, {-2.50000, -0.30000,  1.12500},
    {-2.65000, -0.30000,  0.93750}, {-2.65000,  0.00000,  0.93750},
    {-2.00000, -0.30000,  0.90000}, {-1.90000, -0.30000,  0.60000},
    {-1.90000,  0.00000,  0.60000}, {-3.00000,  0.30000,  1.35000},
    {-2.70000,  0.30000,  1.57500}, {-2.65000,  0.30000,  0.93750},
    {-2.50000,  0.30000,  1.12500}, {-1.90000,  0.30000,  0.60000},
    {-2.00000,  0.30000,  0.90000}, { 1.70000,  0.00000,  1.42500},
    { 1.70000, -0.66000,  1.42500}, { 1.70000, -0.66000,  0.60000},
    { 1.70000,  0.00000,  0.60000}, { 2.60000,  0.00000,  1.42500},
    { 2.60000, -0.66000,  1.42500}, { 3.10000, -0.66000,  0.82500},
    { 3.10000,  0.00000,  0.82500}, { 2.30000,  0.00000,  2.10000},
    { 2.30000, -0.25000,  2.10000}, { 2.40000, -0.25000,  2.02500},
    { 2.40000,  0.00000,  2.02500}, { 2.70000,  0.00000,  2.40000},
    { 2.70000, -0.25000,  2.40000}, { 3.30000, -0.25000,  2.40000},
    { 3.30000,  0.00000,  2.40000}, { 1.70000,  0.66000,  0.60000},
    { 1.70000,  0.66000,  1.42500}, { 3.10000,  0.66000,  0.82500},
    { 2.60000,  0.66000,  1.42600}, { 2.40000,  0.25000,  2.02500},
    { 2.30000,  0.25000,  2.10000}, { 3.30000,  0.25000,  2.40000},
    { 2.70000,  0.25000,  2.40000}, { 2.80000,  0.00000,  2.47500},
    { 2.80000, -0.25000,  2.47500}, { 3.52500, -0.25000,  2.49375},
    { 3.52500,  0.00000,  2.49375}, { 2.90000,  0.00000,  2.47500},
    { 2.90000, -0.15000,  2.47500}, { 3.45000, -0.15000,  2.51250},
    { 3.45000,  0.00000,  2.51250}, { 2.80000,  0.00000,  2.40000},
    { 2.80000, -0.15000,  2.40000}, { 3.20000, -0.15000,  2.40000},
    { 3.20000,  0.00000,  2.40000}, { 3.52500,  0.25000,  2.49375},
    { 2.80000,  0.25000,  2.47500}, { 3.45000,  0.15000,  2.51250},
    { 2.90000,  0.15000,  2.47500}, { 3.20000,  0.15000,  2.40000},
    { 2.80000,  0.15000,  2.40000}, { 0.00000,  0.00000,  3.15000},
    { 0.00000, -0.00200,  3.15000}, { 0.00200,  0.00000,  3.15000},
    { 0.80000,  0.00000,  3.15000}, { 0.80000, -0.45000,  3.15000},
    { 0.45000, -0.80000,  3.15000}, { 0.00000, -0.80000,  3.15000},
    { 0.00000,  0.00000,  2.85000}, { 0.20000,  0.00000,  2.70000},
    { 0.20000, -0.11200,  2.70000}, { 0.11200, -0.20000,  2.70000},
    { 0.00000, -0.20000,  2.70000}, {-0.00200,  0.00000,  3.15000},
    {-0.45000, -0.80000,  3.15000}, {-0.80000, -0.45000,  3.15000},
    {-0.80000,  0.00000,  3.15000}, {-0.11200, -0.20000,  2.70000},
    {-0.20000, -0.11200,  2.70000}, {-0.20000,  0.00000,  2.70000},
    { 0.00000,  0.00200,  3.15000}, {-0.80000,  0.45000,  3.15000},
    {-0.45000,  0.80000,  3.15000}, { 0.00000,  0.80000,  3.15000},
    {-0.20000,  0.11200,  2.70000}, {-0.11200,  0.20000,  2.70000},
    { 0.00000,  0.20000,  2.70000}, { 0.45000,  0.80000,  3.15000},
    { 0.80000,  0.45000,  3.15000}, { 0.11200,  0.20000,  2.70000},
    { 0.20000,  0.11200,  2.70000}, { 0.40000,  0.00000,  2.55000},
    { 0.40000, -0.22400,  2.55000}, { 0.22400, -0.40000,  2.55000},
    { 0.00000, -0.40000,  2.55000}, { 1.30000,  0.00000,  2.55000},
    { 1.30000, -0.72800,  2.55000}, { 0.72800, -1.30000,  2.55000},
    { 0.00000, -1.30000,  2.55000}, { 1.30000,  0.00000,  2.40000},
    { 1.30000, -0.72800,  2.40000}, { 0.72800, -1.30000,  2.40000},
    { 0.00000, -1.30000,  2.40000}, {-0.22400, -0.40000,  2.55000},
    {-0.40000, -0.22400,  2.55000}, {-0.40000,  0.00000,  2.55000},
    {-0.72800, -1.30000,  2.55000}, {-1.30000, -0.72800,  2.55000},
    {-1.30000,  0.00000,  2.55000}, {-0.72800, -1.30000,  2.40000},
    {-1.30000, -0.72800,  2.40000}, {-1.30000,  0.00000,  2.40000},
    {-0.40000,  0.22400,  2.55000}, {-0.22400,  0.40000,  2.55000},
    { 0.00000,  0.40000,  2.55000}, {-1.30000,  0.72800,  2.55000},
    {-0.72800,  1.30000,  2.55000}, { 0.00000,  1.30000,  2.55000},
    {-1.30000,  0.72800,  2.40000}, {-0.72800,  1.30000,  2.40000},
    { 0.00000,  1.30000,  2.40000}, { 0.22400,  0.40000,  2.55000},
    { 0.40000,  0.22400,  2.55000}, { 0.72800,  1.30000,  2.55000},
    { 1.30000,  0.72800,  2.55000}, { 0.72800,  1.30000,  2.40000},
    { 1.30000,  0.72800,  2.40000}, { 0.00000,  0.00000,  0.00000},
    { 1.50000,  0.00000,  0.15000}, { 1.50000,  0.84000,  0.15000},
    { 0.84000,  1.50000,  0.15000}, { 0.00000,  1.50000,  0.15000},
    { 1.50000,  0.00000,  0.07500}, { 1.50000,  0.84000,  0.07500},
    { 0.84000,  1.50000,  0.07500}, { 0.00000,  1.50000,  0.07500},
    { 1.42500,  0.00000,  0.00000}, { 1.42500,  0.79800,  0.00000},
    { 0.79800,  1.42500,  0.00000}, { 0.00000,  1.42500,  0.00000},
    {-0.84000,  1.50000,  0.15000}, {-1.50000,  0.84000,  0.15000},
    {-1.50000,  0.00000,  0.15000}, {-0.84000,  1.50000,  0.07500},
    {-1.50000,  0.84000,  0.07500}, {-1.50000,  0.00000,  0.07500},
    {-0.79800,  1.42500,  0.00000}, {-1.42500,  0.79800,  0.00000},
    {-1.42500,  0.00000,  0.00000}, {-1.50000, -0.84000,  0.15000},
    {-0.84000, -1.50000,  0.15000}, { 0.00000, -1.50000,  0.15000},
    {-1.50000, -0.84000,  0.07500}, {-0.84000, -1.50000,  0.07500},
    { 0.00000, -1.50000,  0.07500}, {-1.42500, -0.79800,  0.00000},
    {-0.79800, -1.42500,  0.00000}, { 0.00000, -1.42500,  0.00000},
    { 0.84000, -1.50000,  0.15000}, { 1.50000, -0.84000,  0.15000},
    { 0.84000, -1.50000,  0.07500}, { 1.50000, -0.84000,  0.07500},
    { 0.79800, -1.42500,  0.00000}, { 1.42500, -0.79800,  0.00000}};
  static const int TeapotBPatches[32][16] =
      /* body */
   {{  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16},
    {  4, 17, 18, 19,  8, 20, 21, 22, 12, 23, 24, 25, 16, 26, 27, 28},
    { 19, 29, 30, 31, 22, 32, 33, 34, 25, 35, 36, 37, 28, 38, 39, 40},
    { 31, 41, 42,  1, 34, 43, 44,  5, 37, 45, 46,  9, 40, 47, 48, 13},
    { 13, 14, 15, 16, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60},
    { 16, 26, 27, 28, 52, 61, 62, 63, 56, 64, 65, 66, 60, 67, 68, 69},
    { 28, 38, 39, 40, 63, 70, 71, 72, 66, 73, 74, 75, 69, 76, 77, 78},
    { 40, 47, 48, 13, 72, 79, 80, 49, 75, 81, 82, 53, 78, 83, 84, 57},
    { 57, 58, 59, 60, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96},
    { 60, 67, 68, 69, 88, 97, 98, 99, 92,100,101,102, 96,103,104,105},
    { 69, 76, 77, 78, 99,106,107,108,102,109,110,111,105,112,113,114},
    { 78, 83, 84, 57,108,115,116, 85,111,117,118, 89,114,119,120, 93},
      /* handle */
    {121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136},
    {124,137,138,121,128,139,140,125,132,141,142,129,136,143,144,133},
    {133,134,135,136,145,146,147,148,149,150,151,152, 69,153,154,155},
    {136,143,144,133,148,156,157,145,152,158,159,149,155,160,161, 69},
      /* spout */
    {162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177},
    {165,178,179,162,169,180,181,166,173,182,183,170,177,184,185,174},
    {174,175,176,177,186,187,188,189,190,191,192,193,194,195,196,197},
    {177,184,185,174,189,198,199,186,193,200,201,190,197,202,203,194},
      /* lid */
    {204,204,204,204,207,208,209,210,211,211,211,211,212,213,214,215},
    {204,204,204,204,210,217,218,219,211,211,211,211,215,220,221,222},
    {204,204,204,204,219,224,225,226,211,211,211,211,222,227,228,229},
    {204,204,204,204,226,230,231,207,211,211,211,211,229,232,233,212},
    {212,213,214,215,234,235,236,237,238,239,240,241,242,243,244,245},
    {215,220,221,222,237,246,247,248,241,249,250,251,245,252,253,254},
    {222,227,228,229,248,255,256,257,251,258,259,260,254,261,262,263},
    {229,232,233,212,257,264,265,234,260,266,267,238,263,268,269,242},
      /* bottom */
    {270,270,270,270,279,280,281,282,275,276,277,278, 93,120,119,114},
    {270,270,270,270,282,289,290,291,278,286,287,288,114,113,112,105},
    {270,270,270,270,291,298,299,300,288,295,296,297,105,104,103, 96},
    {270,270,270,270,300,305,306,279,297,303,304,275, 96, 95, 94, 93}};

  int i, j;
  point3d cnet[16], p;

  for ( i = 0; i < 32; i++ ) {
    for ( j = 0; j < 16; j++ ) {
      MultVector3d ( 0.5*r, &TeapotPoints[TeapotBPatches[i][j]-1], &p );
      SetPoint3d ( &cnet[j], p.x, p.z-0.75, -p.y );
    }
    glMap2d ( GL_MAP2_VERTEX_3, 0.0, 1.0, 3, 4, 0.0, 1.0, 12, 4,
              &cnet[0].x );
    glEnable ( GL_MAP2_VERTEX_3 );
    glEnable ( GL_AUTO_NORMAL );
    glMapGrid2d ( 10, 0.0, 1.0, 10, 0.0, 1.0 );
    glEvalMesh2 ( GL_FILL, 0, 10, 0, 10 );
    glDisable ( GL_AUTO_NORMAL );
    glDisable ( GL_MAP2_VERTEX_3 );
  }
} /*glutSolidTeapot*/

