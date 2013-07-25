
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

#include "testgraphd.h"

/*
#define DUMP
*/

/* these procedures obtain a Bezier patch in R^3; it describes a */
/* graph of a function of two variables: z = f(x,y). They are intended */
/* to compute a number of points of the graph of the laplacian of f */

double ComputeLaplacian ( int n, int m, const point2d *xycp, const double *zcp,
                         double u, double v )
{
  vector2d xy, xyu, xyv, xyuu, xyuv, xyvv;
  vector2d gx, gy, gxx, gxy, gyy;
  double    xu, xv, xuu, xuv, xvv, yu, yv, yuu, yuv, yvv,
            z, zu, zv, zuu, zuv, zvv;
  double    lap, jac, a[9], b[6], c[6], trd[5];

  mbs_BCHornerDer2P1d ( n, m, zcp, u, v, &z, &zu, &zv, &zuu, &zuv, &zvv );
  mbs_BCHornerDer2P2d ( n, m, xycp, u, v,
                        &xy, &xyu, &xyv, &xyuu, &xyuv, &xyvv );

  xu = xyu.x;    yu = xyu.y;
  xv = xyv.x;    yv = xyv.y;
  xuu = xyuu.x;  yuu = xyuu.y;
  xuv = xyuv.x;  yuv = xyuv.y;
  xvv = xyvv.x;  yvv = xyvv.y;

  jac = fabs ( xu*yv-yu*xv );
  a[0] = xu;  a[1] = yu;  b[0] = 1.0;  b[1] = 0.0;
  a[2] = xv;  a[3] = yv;  b[2] = 0.0;  b[3] = 1.0;
  pkn_multiSolveRLSQd ( 2, 2, a, 2, 2, b, 2, c );
  gx.x = c[0];  gx.y = c[1];
  gy.x = c[2];  gy.y = c[3];

/* here c[0]*zu + c[1]*zv = zx;  c[2]*zu + c[3]*zv = zy */

  a[0] = xu*xu;  a[1] = 2.0*xu*yu;    a[2] = yu*yu;
  a[3] = xu*xv;  a[4] = xu*yv+xv*yu;  a[5] = yu*yv;
  a[6] = xv*xv;  a[7] = 2.0*xv*yv;    a[8] = yv*yv;
  b[0] = -(xuu*gx.x+yuu*gy.x);  b[1] = -(xuu*gx.y+yuu*gy.y);
  b[2] = -(xuv*gx.x+yuv*gy.x);  b[3] = -(xuv*gx.y+yuv*gy.y);
  b[4] = -(xvv*gx.x+yvv*gy.x);  b[5] = -(xvv*gx.y+yvv*gy.y);
  pkn_multiSolveRLSQd ( 3, 3, a, 2, 2, b, 2, c );
  gxx.x = c[0];  gxx.y = c[1];
  gxy.x = c[2];  gxy.y = c[3];
  gyy.x = c[4];  gyy.y = c[5];

  trd[0]  = gxx.x;      trd[1]  = gxx.y;
  trd[2]  = gx.x*gx.x;  trd[3]  = 2.0*gx.x*gx.y;  trd[4] = gx.y*gx.y;
  trd[0] += gyy.x;      trd[1] += gyy.y;
  trd[2] += gy.x*gy.x;  trd[3] += 2.0*gy.x*gy.y;  trd[4] += gy.y*gy.y;

  lap = trd[0]*zu + trd[1]*zv + trd[2]*zuu + trd[3]*zuv + trd[4]*zvv;
  return lap;
} /*ComputeLaplacian*/

void ComputeLaplacianGrad ( int n, int m,
                            const point2d *xycp, const double *zcp,
                            double u, double v, double *jac, vector2d *lgr )
{
#ifdef DUMP
  FILE *f;
#endif
  vector2d xy, xyu, xyv, xyuu, xyuv, xyvv, xyuuu, xyuuv, xyuvv, xyvvv;
  vector2d gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy;
  double   xu, xv, xuu, xuv, xvv, xuuu, xuuv, xuvv, xvvv,
           yu, yv, yuu, yuv, yvv, yuuu, yuuv, yuvv, yvvv,
           z, zu, zv, zuu, zuv, zvv, zuuu, zuuv, zuvv, zvvv;
  double   a[16], b[8], c[8], trd[18];
#ifdef DUMP
  f = fopen ( "lapgr.dat", "w+" );
  fprintf ( f, "n=%d, m=%d, u=%f, v=%f\n", n, m, u, v );
#endif
  mbs_BCHornerDer3P1d ( n, m, zcp, u, v, &z, &zu, &zv, &zuu, &zuv, &zvv,
                       &zuuu, &zuuv, &zuvv, &zvvv );
  mbs_BCHornerDer3P2d ( n, m, (double*)xycp, u, v,
                       (double*)&xy, (double*)&xyu, (double*)&xyv,
                       (double*)&xyuu, (double*)&xyuv, (double*)&xyvv,
                       (double*)&xyuuu, (double*)&xyuuv, (double*)&xyuvv,
                       (double*)&xyvvv );
#ifdef DUMP
  fprintf ( f, "z=%f, zu=%f, zv=%f,\nzuu=%f, zuv=%f, zvv=%f\n",
            z, zu, zv, zuu, zuv, zvv );
  fprintf ( f, "zuuu=%f, zuuv=%f,\nzuvv=%f, zvvv=%f\n",
            zuuu, zuuv, zuvv, zvvv );
#endif
  xu = xyu.x;    yu = xyu.y;
  xv = xyv.x;    yv = xyv.y;
  xuu = xyuu.x;  yuu = xyuu.y;
  xuv = xyuv.x;  yuv = xyuv.y;
  xvv = xyvv.x;  yvv = xyvv.y;
  xuuu = xyuuu.x;  yuuu = xyuuu.y;
  xuuv = xyuuv.x;  yuuv = xyuuv.y;
  xuvv = xyuvv.x;  yuvv = xyuvv.y;
  xvvv = xyvvv.x;  yvvv = xyvvv.y;
#ifdef DUMP
  fprintf ( f, "x=%f, xu=%f, xv=%f,\nxuu=%f, xuv=%f, xvv=%f\n",
            xy.x, xu, xv, xuu, xuv, xvv );
  fprintf ( f, "y=%f, yu=%f, yv=%f,\nyuu=%f, yuv=%f, yvv=%f\n",
            xy.y, yu, yv, yuu, yuv, yvv );
  fprintf ( f, "xuuu=%f, xuuv=%f, xuvv=%f, xvvv=%f\n",
            xuuu, xuuv, xuvv, xvvv );
  fprintf ( f, "yuuu=%f, yuuv=%f, yuvv=%f, yvvv=%f\n",
            yuuu, yuuv, yuvv, yvvv );
#endif
  *jac = fabs ( xu*yv-yu*xv );
#ifdef DUMP
  fprintf ( f, "jac=%f\n", *jac );
#endif
  a[0] = xu;  a[1] = yu;  b[0] = 1.0;  b[1] = 0.0;
  a[2] = xv;  a[3] = yv;  b[2] = 0.0;  b[3] = 1.0;
  pkn_multiSolveRLSQd ( 2, 2, a, 2, 2, b, 2, c );
  gx.x = c[0];  gx.y = c[1];
  gy.x = c[2];  gy.y = c[3];

/* here c[0]*zu + c[1]*zv = zx;  c[2]*zu + c[3]*zv = zy */

  a[0] = xu*xu;  a[1] = 2.0*xu*yu;    a[2] = yu*yu;
  a[3] = xu*xv;  a[4] = xu*yv+xv*yu;  a[5] = yu*yv;
  a[6] = xv*xv;  a[7] = 2.0*xv*yv;    a[8] = yv*yv;
  b[0] = -(xuu*gx.x+yuu*gy.x);  b[1] = -(xuu*gx.y+yuu*gy.y);
  b[2] = -(xuv*gx.x+yuv*gy.x);  b[3] = -(xuv*gx.y+yuv*gy.y);
  b[4] = -(xvv*gx.x+yvv*gy.x);  b[5] = -(xvv*gx.y+yvv*gy.y);
  pkn_multiSolveRLSQd ( 3, 3, a, 2, 2, b, 2, c );
  gxx.x = c[0];  gxx.y = c[1];
  gxy.x = c[2];  gxy.y = c[3];
  gyy.x = c[4];  gyy.y = c[5];

/* here is the code to compute the laplacian */
/*
  trd[0]  = gxx.x;      trd[1]  = gxx.y;
  trd[2]  = gx.x*gx.x;  trd[3]  = 2.0*gx.x*gx.y;  trd[4] = gx.y*gx.y;
  trd[0] += gyy.x;      trd[1] += gyy.y;
  trd[2] += gy.x*gy.x;  trd[3] += 2.0*gy.x*gy.y;  trd[4] += gy.y*gy.y;

  lap = trd[0]*zu + trd[1]*zv + trd[2]*zuu + trd[3]*zuv + trd[4]*zvv;
*/

        /* third order derivatives */
  a[0]  = xu*xu*xu;  a[1]  = 3.0*xu*xu*yu;           a[2]  = 3.0*xu*yu*yu;           a[3]  = yu*yu*yu;
  a[4]  = xu*xu*xv;  a[5]  = 2.0*xu*xv*yu+xu*xu*yv;  a[6]  = 2.0*xu*yu*yv+xv*yu*yu;  a[7]  = yu*yu*yv;
  a[8]  = xu*xv*xv;  a[9]  = 2.0*xu*xv*yv+xv*xv*yu;  a[10] = 2.0*xv*yu*yv+xu*yv*yv;  a[11] = yu*yv*yv;
  a[12] = xv*xv*xv;  a[13] = 3.0*xv*xv*yv;           a[14] = 3.0*xv*yv*yv;           a[15] = yv*yv*yv;
  b[0] = -(3.0*xu*xuu*gxx.x + 3.0*(xuu*yu+xu*yuu)*gxy.x + 3.0*yu*yuu*gyy.x + xuuu*gx.x + yuuu*gy.x);
  b[1] = -(3.0*xu*xuu*gxx.y + 3.0*(xuu*yu+xu*yuu)*gxy.y + 3.0*yu*yuu*gyy.y + xuuu*gx.y + yuuu*gy.y);
  b[2] = -((xuu*xv+2.0*xu*xuv)*gxx.x + (xuu*yv+2.0*(xu*yuv+xuv*yu)+xv*yuu)*gxy.x + (yuu*yv+2.0*yu*yuv)*gyy.x + xuuv*gx.x + yuuv*gy.x);
  b[3] = -((xuu*xv+2.0*xu*xuv)*gxx.y + (xuu*yv+2.0*(xu*yuv+xuv*yu)+xv*yuu)*gxy.y + (yuu*yv+2.0*yu*yuv)*gyy.y + xuuv*gx.y + yuuv*gy.y);
  b[4] = -((xvv*xu+2.0*xv*xuv)*gxx.x + (xvv*yu+2.0*(xv*yuv+xuv*yv)+xu*yvv)*gxy.x + (yu*yvv+2.0*yv*yuv)*gyy.x + xuvv*gx.x + yuvv*gy.x);
  b[5] = -((xvv*xu+2.0*xv*xuv)*gxx.y + (xvv*yu+2.0*(xv*yuv+xuv*yv)+xu*yvv)*gxy.y + (yu*yvv+2.0*yv*yuv)*gyy.y + xuvv*gx.y + yuvv*gy.y);
  b[6] = -(3.0*xv*xvv*gxx.x + 3.0*(xvv*yv+xv*yvv)*gxy.x + 3.0*yv*yvv*gyy.x + xvvv*gx.x + yvvv*gy.x);
  b[7] = -(3.0*xv*xvv*gxx.y + 3.0*(xvv*yv+xv*yvv)*gxy.y + 3.0*yv*yvv*gyy.y + xvvv*gx.y + yvvv*gy.y);
  pkn_multiSolveRLSQd ( 4, 4, a, 2, 2, b, 2, c );
  gxxx.x = c[0];  gxxx.y = c[1];
  gxxy.x = c[2];  gxxy.y = c[3];
  gxyy.x = c[4];  gxyy.y = c[5];
  gyyy.x = c[6];  gyyy.y = c[7];

  trd[ 0] = gxxx.x + gxyy.x;
  trd[ 1] = gxxx.y + gxyy.y;
  trd[ 2] = 3.0*gx.x*gxx.x + (gyy.x*gx.x+2.0*gy.x*gxy.x);
  trd[ 3] = 3.0*(gxx.x*gx.y+gx.x*gxx.y) +
            (gyy.x*gx.y+2.0*(gy.x*gxy.y+gxy.x*gy.y)+gx.x*gyy.y);
  trd[ 4] = 3.0*gx.y*gxx.y + (gyy.y*gx.y+2.0*gy.y*gxy.y);
  trd[ 5] = gx.x*gx.x*gx.x + gx.x*gy.x*gy.x;
  trd[ 6] = 3.0*gx.x*gx.x*gx.y + (2.0*gx.x*gy.x*gy.y+gy.x*gy.x*gx.y);
  trd[ 7] = 3.0*gx.x*gx.y*gx.y + (2.0*gy.x*gx.y*gy.y+gx.x*gy.y*gy.y);
  trd[ 8] = gx.y*gx.y*gx.y + gx.y*gy.y*gy.y;
  trd[ 9] = gxxy.x + gyyy.x;
  trd[10] = gxxy.y + gyyy.y;
  trd[11] = (gxx.x*gy.x+2.0*gx.x*gxy.x) + 3.0*gy.x*gyy.x;
  trd[12] = (gxx.x*gy.y+2.0*(gx.x*gxy.y+gxy.x*gx.y)+gy.x*gxx.y) +
            3.0*(gyy.x*gy.y+gy.x*gyy.y);
  trd[13] = (gxx.y*gy.y+2.0*gx.y*gxy.y) + 3.0*gy.y*gyy.y;
  trd[14] = gx.x*gx.x*gy.x + gy.x*gy.x*gy.x;
  trd[15] = (2.0*gx.x*gy.x*gx.y+gx.x*gx.x*gy.y) + 3.0*gy.x*gy.x*gy.y;
  trd[16] = (2.0*gx.x*gx.y*gy.y+gy.x*gx.y*gx.y) + 3.0*gy.x*gy.y*gy.y;
  trd[17] = gx.y*gx.y*gy.y + gy.y*gy.y*gy.y;

  lgr->x = trd[ 0]*zu + trd[ 1]*zv + trd[ 2]*zuu + trd[ 3]*zuv + trd[ 4]*zvv +
           trd[ 5]*zuuu + trd[ 6]*zuuv + trd[ 7]*zuvv + trd[ 8]*zvvv;
  lgr->y = trd[ 9]*zu + trd[10]*zv + trd[11]*zuu + trd[12]*zuv + trd[13]*zvv +
           trd[14]*zuuu + trd[15]*zuuv + trd[16]*zuvv + trd[17]*zvvv;
#ifdef DUMP
  fprintf ( f, "lgr=(%f,%f)\n", lgr->x, lgr->y );
  fclose ( f );
  exit ( 0 );
#endif
} /*ComputeLaplacianGrad*/

void TestTabLaplacianU ( int n, int m, const point3d *cp, double u, int dd,
                         point3d *gp )
{
  void     *sp;
  point2d  *xycp;
  double    *zcp;
  double    v;
  int      i;

  sp   = pkv_GetScratchMemTop ();
  xycp = (point2d*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  zcp  = pkv_GetScratchMemd ( (n+1)*(m+1) );
  if ( !xycp || !zcp )
    exit ( 1 );

  pkv_Selectd ( (n+1)*(m+1), 2, 3, 2, cp, xycp );
  pkv_Selectd ( (n+1)*(m+1), 1, 3, 1, &cp[0].z, zcp );

  for ( i = 0; i <= dd; i++ ) {
    v = (double)i/dd;
    mbs_BCHornerP2d ( n, m, (double*)xycp, u, v, &gp[i].x );
    gp[i].z = ComputeLaplacian ( n, m, xycp, zcp, u, v );
  }

  pkv_SetScratchMemTop ( sp );
} /*TestTabLaplacianU*/

void TestTabLaplacianV ( int n, int m, const point3d *cp, double v, int dd,
                         point3d *gp )
{
  void     *sp;
  point2d  *xycp;
  double    *zcp;
  double    u;
  int      i;

  sp   = pkv_GetScratchMemTop ();
  xycp = (point2d*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  zcp  = pkv_GetScratchMemd ( (n+1)*(m+1) );
  if ( !xycp || !zcp )
    exit ( 1 );

  pkv_Selectd ( (n+1)*(m+1), 2, 3, 2, cp, xycp );
  pkv_Selectd ( (n+1)*(m+1), 1, 3, 1, &cp[0].z, zcp );

  for ( i = 0; i <= dd; i++ ) {
    u = (double)i/dd;
    mbs_BCHornerP2d ( n, m, (double*)xycp, u, v, &gp[i].x );
    gp[i].z = ComputeLaplacian ( n, m, xycp, zcp, u, v );
  }

  pkv_SetScratchMemTop ( sp );
} /*TestTabLaplacianV*/

void TestTabLapGradU ( int n, int m, const point3d *cp, double u, int dd,
                       char sxy, point3d *gp )
{
  void     *sp;
  point2d  *xycp;
  double    *zcp;
  double    v, jac;
  int      i;
  vector2d lgr;

  sp   = pkv_GetScratchMemTop ();
  xycp = (point2d*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  zcp  = pkv_GetScratchMemd ( (n+1)*(m+1) );
  if ( !xycp || !zcp )
    exit ( 1 );

  pkv_Selectd ( (n+1)*(m+1), 2, 3, 2, cp, xycp );
  pkv_Selectd ( (n+1)*(m+1), 1, 3, 1, &cp[0].z, zcp );

  for ( i = 0; i <= dd; i++ ) {
    v = (double)i/dd;
    mbs_BCHornerP2d ( n, m, (double*)xycp, u, v, &gp[i].x );
    ComputeLaplacianGrad ( n, m, xycp, zcp, u, v, &jac, &lgr );
    if ( sxy == 0 ) gp[i].z = lgr.x;  else gp[i].z = lgr.y;
  }

  pkv_SetScratchMemTop ( sp );
} /*TestTabLapGradU*/

void TestTabLapGradV ( int n, int m, const point3d *cp, double v, int dd,
                       char sxy, point3d *gp )
{
  void     *sp;
  point2d  *xycp;
  double    *zcp;
  double    u, jac;
  int      i;
  vector2d lgr;

  sp   = pkv_GetScratchMemTop ();
  xycp = (point2d*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  zcp  = pkv_GetScratchMemd ( (n+1)*(m+1) );
  if ( !xycp || !zcp )
    exit ( 1 );

  pkv_Selectd ( (n+1)*(m+1), 2, 3, 2, cp, xycp );
  pkv_Selectd ( (n+1)*(m+1), 1, 3, 1, &cp[0].z, zcp );

  for ( i = 0; i <= dd; i++ ) {
    u = (double)i/dd;
    mbs_BCHornerP2d ( n, m, (double*)xycp, u, v, &gp[i].x );
    ComputeLaplacianGrad ( n, m, xycp, zcp, u, v, &jac, &lgr );
    if ( sxy == 0 ) gp[i].z = lgr.x;  else gp[i].z = lgr.y;
  }

  pkv_SetScratchMemTop ( sp );
} /*TestTabLapGradV*/

