
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

#include "msgpool.h"

/* //////////////////////////////////////////////////////////////////// */
/* the procedures in this file set up the matrices, which appear in the */
/* formulae expressing the partial derivatives up to the order 4        */
/* of a composition of a function R^2 -> R^2, with another function.    */

void pkn_Setup2DerA11Matrixd ( double xu, double yu, double xv, double yv,
                               double *A11 )
{
  A11[0] = xu;  A11[1] = yu;
  A11[2] = xv;  A11[3] = yv;
} /*pkn_Setup2DerA11Matrixd*/

void pkn_Setup2DerA21Matrixd ( double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A21 )
{
  A21[0] = xuu;  A21[1] = yuu;
  A21[2] = xuv;  A21[3] = yuv;
  A21[4] = xvv;  A21[5] = yvv;
} /*pkn_Setup2DerA21Matrixd*/

void pkn_Setup2DerA22Matrixd ( double xu, double yu, double xv, double yv,
                               double *A22 )
{
  A22[0] = xu*xu;  A22[1] = 2.0*xu*yu;    A22[2] = yu*yu;
  A22[3] = xu*xv;  A22[4] = xu*yv+yu*xv;  A22[5] = yu*yv;
  A22[6] = xv*xv;  A22[7] = 2.0*xv*yv;    A22[8] = yv*yv;
} /*pkn_Setup2DerA22Matrixd*/

void pkn_Setup2DerA31Matrixd ( double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double *A31 )
{
  A31[0] = xuuu;  A31[1] = yuuu;
  A31[2] = xuuv;  A31[3] = yuuv;
  A31[4] = xuvv;  A31[5] = yuvv;
  A31[6] = xvvv;  A31[7] = yvvv;
} /*pkn_Setup2DerA31Matrixd*/

void pkn_Setup2DerA32Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A32 )
{
  A32[0] = 3.0*xu*xuu;  A32[1] = 3.0*(xuu*yu+xu*yuu);  A32[2] = 3.0*yu*yuu;

  A32[3] = xuu*xv+2.0*xu*xuv;  A32[4] = xuu*yv+2.0*(xu*yuv+xuv*yu)+xv*yuu;
  A32[5] = yuu*yv+2.0*yu*yuv;

  A32[6] = xvv*xu+2.0*xv*xuv;  A32[7] = xvv*yu+2.0*(xv*yuv+xuv*yv)+xu*yvv;
  A32[8] = yu*yvv+2.0*yv*yuv;

  A32[9] = 3.0*xv*xvv;  A32[10] = 3.0*(xvv*yv+xv*yvv);  A32[11] = 3.0*yv*yvv;
} /*pkn_Setup2DerA32Matrixd*/

void pkn_Setup2DerA33Matrixd ( double xu, double yu, double xv, double yv,
                               double *A33 )
{
  A33[0] = xu*xu*xu;  A33[1] = 3.0*xu*xu*yu;  A33[2] = 3.0*xu*yu*yu;
  A33[3] = yu*yu*yu;

  A33[4] = xu*xu*xv;  A33[5] = 2.0*xu*xv*yu+xu*xu*yv;
  A33[6] = 2.0*xu*yu*yv+xv*yu*yu;  A33[7] = yu*yu*yv;

  A33[8] = xu*xv*xv;  A33[9] = 2.0*xu*xv*yv+xv*xv*yu;
  A33[10] = 2.0*xv*yu*yv+xu*yv*yv;  A33[11] = yu*yv*yv;

  A33[12] = xv*xv*xv;  A33[13] = 3.0*xv*xv*yv;  A33[14] = 3.0*xv*yv*yv;
  A33[15] = yv*yv*yv;
} /*pkn_Setup2DerA33Matrixd*/

void pkn_Setup2DerA41Matrixd ( double xuuuu, double yuuuu, double xuuuv, double yuuuv,
                               double xuuvv, double yuuvv, double xuvvv, double yuvvv,
                               double xvvvv, double yvvvv,
                               double *A41 )
{
  A41[0] = xuuuu;  A41[1] = yuuuu;
  A41[2] = xuuuv;  A41[3] = yuuuv;
  A41[4] = xuuvv;  A41[5] = yuuvv;
  A41[6] = xuvvv;  A41[7] = yuvvv;
  A41[8] = xvvvv;  A41[9] = yvvvv;
} /*pkn_Setup2DerA41Matrixd*/

void pkn_Setup2DerA42Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double *A42 )
{
  A42[0] = 3.0*xuu*xuu+4.0*xu*xuuu;  A42[1] = 4.0*(xuuu*yu+xu*yuuu)+6.0*xuu*yuu;
  A42[2] = 3.0*yuu*yuu+4.0*yu*yuuu;

  A42[3] = xuuu*xv+3.0*(xuu*xuv+xu*xuuv);
  A42[4] = xuuu*yv+3.0*(xuu*yuv+xu*yuuv+xuuv*yu+xuv*yuu)+xv*yuuu;
  A42[5] = yuuu*yv+3.0*(yuu*yuv+yu*yuuv);

  A42[6] = 2.0*(xv*xuuv+xuv*xuv+xu*xuvv)+xuu*xvv;
  A42[7] = 2.0*(xuuv*yv+xu*yuvv+xuvv*yu+xv*yuuv)+4.0*xuv*yuv+xuu*yvv+xvv*yuu;
  A42[8] = 2.0*(yv*yuuv+yuv*yuv+yu*yuvv)+yuu*yvv;

  A42[9] = xvvv*xu+3.0*(xvv*xuv+xv*xuvv);
  A42[10] = xvvv*yu+3.0*(xvv*yuv+xv*yuvv+xuvv*yv+xuv*yvv)+xu*yvvv;
  A42[11] = yvvv*yu+3.0*(yvv*yuv+yv*yuvv);

  A42[12] = 3.0*xvv*xvv+4.0*xv*xvvv;  A42[13] = 4.0*(xvvv*yv+xv*yvvv)+6.0*xvv*yvv;
  A42[14] = 3.0*yvv*yvv+4.0*yv*yvvv;
} /*pkn_Setup2DerA42Matrixd*/

void pkn_Setup2DerA43Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A43 )
{
  A43[0] = 6.0*xu*xu*xuu;  A43[1] = 12.0*xu*xuu*yu+6.0*xu*xu*yuu;
  A43[2] = 12.0*xu*yu*yuu+6.0*xuu*yu*yu;  A43[3] = 6.0*yu*yu*yuu;

  A43[4] = 3.0*(xu*xuu*xv+xu*xu*xuv);
  A43[5] = 3.0*(xuu*xv*yu+xu*xv*yuu+xu*xuu*yv+xu*xu*yuv)+6.0*xu*xuv*yu;
  A43[6] = 3.0*(xuu*yu*yv+xu*yuu*yv+xuv*yu*yu+xv*yu*yuu)+6.0*xu*yu*yuv;
  A43[7] = 3.0*(yu*yuu*yv+yu*yu*yuv);

  A43[8] = xu*xu*xvv+4.0*xu*xv*xuv+xuu*xv*xv;
  A43[9] = 4.0*(xv*xuv*yu+xu*xuv*yv+xu*xv*yuv)+
           2.0*(xuu*xv*yv+xu*xvv*yu)+xu*xu*yvv+xv*xv*yuu;
  A43[10] = 4.0*(xuv*yu*yv+xu*yuv*yv+xv*yu*yuv)+
            2.0*(xu*yu*yvv+xv*yuu*yv)+xvv*yu*yu+xuu*yv*yv;
  A43[11] = yu*yu*yvv+4.0*yu*yuv*yv+yuu*yv*yv;

  A43[12] = 3.0*(xv*xvv*xu+xv*xv*xuv);
  A43[13] = 3.0*(xvv*xu*yv+xv*xu*yvv+xv*xvv*yu+xv*xv*yuv)+6.0*xv*xuv*yv;
  A43[14] = 3.0*(xvv*yv*yu+xv*yvv*yu+xuv*yv*yv+xu*yv*yvv)+6.0*xv*yv*yuv;
  A43[15] = 3.0*(yv*yvv*yu+yv*yv*yuv);

  A43[16] = 6.0*xv*xv*xvv;  A43[17] = 12.0*xv*xvv*yv+6.0*xv*xv*yvv;
  A43[18] = 12.0*xv*yv*yvv+6.0*xvv*yv*yv;  A43[19] = 6.0*yv*yv*yvv;
} /*pkn_Setup2DerA43Matrixd*/

void pkn_Setup2DerA44Matrixd ( double xu, double yu, double xv, double yv,
                               double *A44 )
{
  A44[0] = xu*xu*xu*xu;  A44[1] = 4.0*xu*xu*xu*yu;
  A44[2] = 6.0*xu*xu*yu*yu;
  A44[3] = 4.0*xu*yu*yu*yu;  A44[4] = yu*yu*yu*yu;

  A44[5] = xu*xu*xu*xv;  A44[6] = 3.0*xu*xu*xv*yu+xu*xu*xu*yv;
  A44[7] = 3.0*(xu*xv*yu*yu+xu*xu*yu*yv);
  A44[8] = 3.0*xu*yu*yu*yv+xv*yu*yu*yu;  A44[9] = yu*yu*yu*yv;

  A44[10] = xu*xu*xv*xv;  A44[11] = 2.0*(xu*xu*xv*yv+xu*xv*xv*yu);
  A44[12] = xu*xu*yv*yv+4.0*xu*xv*yu*yv+xv*xv*yu*yu;
  A44[13] = 2.0*(xu*yu*yv*yv+xv*yu*yu*yv);  A44[14] = yu*yu*yv*yv;

  A44[15] = xv*xv*xv*xu;  A44[16] = 3.0*xv*xv*xu*yv+xv*xv*xv*yu;
  A44[17] = 3.0*(xv*xu*yv*yv+xv*xv*yv*yu);
  A44[18] = 3.0*xv*yv*yv*yu+xu*yv*yv*yv;  A44[19] = yv*yv*yv*yu;

  A44[20] = xv*xv*xv*xv;  A44[21] = 4.0*xv*xv*xv*yv;
  A44[22] = 6.0*xv*xv*yv*yv;
  A44[23] = 4.0*xv*yv*yv*yv;  A44[24] = yv*yv*yv*yv;
} /*pkn_Setup2DerA44Matrixd*/

/* //////////////////////////////////////////////////////////////////// */
/* the procedures below compute the partial derivatives up to the order */
/* 4 of a composition h of a function f: R^2 -> R^2 (), described with  */
/* the scalar functions x, y with a function g: R^2 -> R^d.             */
/* The parameter spdimen specifies d.                                   */

void pkn_Comp2Derivatives1d ( double xu, double yu, double xv, double yv,
                              int spdimen, const double *gx, const double *gy,
                              double *hu, double *hv )
{
  void   *sp;
  double *A11;
  int    i;

  sp = pkv_GetScratchMemTop ();
  A11 = pkv_GetScratchMemd ( 4 );
  if ( !A11 )
    pkv_SignalError ( LIB_PKNUM, 3, ERRMSG_0 );

  pkn_Setup2DerA11Matrixd ( xu, yu, xv, yv, A11 );

      /* manual matrix multiplication instead of a call to a special */
      /* procedure because we know nothing about the values of       */
      /* the pointers gx, gy, hx and hy. */
  for ( i = 0; i < spdimen; i++ ) {
    hu[i] = A11[0]*gx[i] + A11[1]*gy[i];
    hv[i] = A11[2]*gx[i] + A11[3]*gy[i];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2Derivatives1d*/

void pkn_Comp2Derivatives2d ( double xu, double yu, double xv, double yv,
                              double xuu, double yuu, double xuv,
                              double yuv, double xvv, double yvv,
                              int spdimen,
                              const double *gx, const double *gy,
                              const double *gxx, const double *gxy,
                              const double *gyy,
                              double *huu, double *huv, double *hvv )
{
  void   *sp;
  double *A21, *A22;
  int    i;

  sp = pkv_GetScratchMemTop ();
  A21 = pkv_GetScratchMemd ( 15 );
  if ( !A21 )
    pkv_SignalError ( LIB_PKNUM, 4, ERRMSG_0 );
  A22 = &A21[6];

  pkn_Setup2DerA21Matrixd ( xuu, yuu, xuv, yuv, xvv, yvv, A21 );
  pkn_Setup2DerA22Matrixd ( xu, yu, xv, yv, A22 );

  for ( i = 0; i < spdimen; i++ ) {
    huu[i] = A21[0]*gx[i] + A21[1]*gy[i] +
             A22[0]*gxx[i] + A22[1]*gxy[i] + A22[2]*gyy[i];
    huv[i] = A21[2]*gx[i] + A21[3]*gy[i] +
             A22[3]*gxx[i] + A22[4]*gxy[i] + A22[5]*gyy[i];
    hvv[i] = A21[4]*gx[i] + A21[5]*gy[i] +
             A22[6]*gxx[i] + A22[7]*gxy[i] + A22[8]*gyy[i];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2Derivatives2d*/

void pkn_Comp2Derivatives3d ( double xu, double yu, double xv, double yv,
                              double xuu, double yuu, double xuv,
                              double yuv, double xvv, double yvv,
                              double xuuu, double yuuu, double xuuv, double yuuv,
                              double xuvv, double yuvv, double xvvv, double yvvv,
                              int spdimen,
                              const double *gx, const double *gy,
                              const double *gxx, const double *gxy,
                              const double *gyy,
                              const double *gxxx, const double *gxxy,
                              const double *gxyy, const double *gyyy,
                              double *huuu, double *huuv,
                              double *huvv, double *hvvv )
{
  void   *sp;
  double *A31, *A32, *A33;
  int    i;

  sp = pkv_GetScratchMemTop ();
  A31 = pkv_GetScratchMemd ( 36 );
  if ( !A31 )
    pkv_SignalError ( LIB_PKNUM, 5, ERRMSG_0 );
  A32 = &A31[8];  A33 = &A32[12];

  pkn_Setup2DerA31Matrixd ( xuuu, yuuu, xuuv, yuuv,
                            xuvv, yuvv, xvvv, yvvv, A31 );
  pkn_Setup2DerA32Matrixd ( xu, yu, xv, yv,
                            xuu, yuu, xuv, yuv, xvv, yvv, A32 );
  pkn_Setup2DerA33Matrixd ( xu, yu, xv, yv, A33 );

  for ( i = 0; i < spdimen; i++ ) {
    huuu[i] = A31[0]*gx[i] + A31[1]*gy[i] +
              A32[0]*gxx[i] + A32[1]*gxy[i] + A32[2]*gyy[i] +
              A33[0]*gxxx[i] + A33[1]*gxxy[i] + A33[2]*gxyy[i] + A33[3]*gyyy[i];
    huuv[i] = A31[2]*gx[i] + A31[3]*gy[i] +
              A32[3]*gxx[i] + A32[4]*gxy[i] + A32[5]*gyy[i] +
              A33[4]*gxxx[i] + A33[5]*gxxy[i] + A33[6]*gxyy[i] + A33[7]*gyyy[i];
    huvv[i] = A31[4]*gx[i] + A31[5]*gy[i] +
              A32[6]*gxx[i] + A32[7]*gxy[i] + A32[8]*gyy[i] +
              A33[8]*gxxx[i] + A33[9]*gxxy[i] + A33[10]*gxyy[i] + A33[11]*gyyy[i];
    hvvv[i] = A31[6]*gx[i] + A31[7]*gy[i] +
              A32[9]*gxx[i] + A32[10]*gxy[i] + A32[11]*gyy[i] +
              A33[12]*gxxx[i] + A33[13]*gxxy[i] + A33[14]*gxyy[i] + A33[15]*gyyy[i];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2Derivatives3d*/

void pkn_Comp2Derivatives4d ( double xu, double yu, double xv, double yv,
                              double xuu, double yuu, double xuv,
                              double yuv, double xvv, double yvv,
                              double xuuu, double yuuu, double xuuv, double yuuv,
                              double xuvv, double yuvv, double xvvv, double yvvv,
                              double xuuuu, double yuuuu, double xuuuv,
                              double yuuuv, double xuuvv, double yuuvv,
                              double xuvvv, double yuvvv,
                              double xvvvv, double yvvvv,
                              int spdimen,
                              const double *gx, const double *gy,
                              const double *gxx, const double *gxy,
                              const double *gyy,
                              const double *gxxx, const double *gxxy,
                              const double *gxyy, const double *gyyy,
                              const double *gxxxx, const double *gxxxy,
                              const double *gxxyy, const double *gxyyy,
                              const double *gyyyy,
                              double *huuuu, double *huuuv, double *huuvv,
                              double *huvvv, double *hvvvv )
{
  void   *sp;
  double *A41, *A42, *A43, *A44;
  int    i;

  sp = pkv_GetScratchMemTop ();
  A41 = pkv_GetScratchMemd ( 70 );
  if ( !A41 )
    pkv_SignalError ( LIB_PKNUM, 6, ERRMSG_0 );
  A42 = &A41[10];  A43 = &A42[15];  A44 = &A43[20];

  pkn_Setup2DerA41Matrixd ( xuuuu, yuuuu, xuuuv, yuuuv, xuuvv, yuuvv,
                            xuvvv, yuvvv, xvvvv, yvvvv, A41);
  pkn_Setup2DerA42Matrixd ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            xuuu, yuuu, xuuv, yuuv, xuvv, yuvv, xvvv, yvvv,
                            A42 );
  pkn_Setup2DerA43Matrixd ( xu, yu, xv, yv,
                            xuu, yuu, xuv, yuv, xvv, yvv, A43 );
  pkn_Setup2DerA44Matrixd ( xu, yu, xv, yv, A44 );

  for ( i = 0; i < spdimen; i++ ) {
    huuuu[i] = A41[0]*gx[i] + A41[1]*gy[i] +
               A42[0]*gxx[i] + A42[1]*gxy[i] + A42[2]*gyy[i] +
               A43[0]*gxxx[i] + A43[1]*gxxy[i] + A43[2]*gxyy[i] + A43[3]*gyyy[i] +
               A44[0]*gxxxx[i] + A44[1]*gxxxy[i] + A44[2]*gxxyy[i] +
               A44[3]*gxyyy[i] + A44[4]*gyyyy[i];
    huuuv[i] = A41[2]*gx[i] + A41[3]*gy[i] +
               A42[3]*gxx[i] + A42[4]*gxy[i] + A42[5]*gyy[i] +
               A43[4]*gxxx[i] + A43[5]*gxxy[i] + A43[6]*gxyy[i] + A43[7]*gyyy[i] +
               A44[5]*gxxxx[i] + A44[6]*gxxxy[i] + A44[7]*gxxyy[i] +
               A44[8]*gxyyy[i] + A44[9]*gyyyy[i];
    huuvv[i] = A41[4]*gx[i] + A41[5]*gy[i] +
               A42[6]*gxx[i] + A42[7]*gxy[i] + A42[8]*gyy[i] +
               A43[8]*gxxx[i] + A43[9]*gxxy[i] + A43[10]*gxyy[i] + A43[11]*gyyy[i] +
               A44[10]*gxxxx[i] + A44[11]*gxxxy[i] + A44[12]*gxxyy[i] +
               A44[13]*gxyyy[i] + A44[14]*gyyyy[i];
    huvvv[i] = A41[6]*gx[i] + A41[7]*gy[i] +
               A42[9]*gxx[i] + A42[10]*gxy[i] + A42[11]*gyy[i] +
               A43[12]*gxxx[i] + A43[13]*gxxy[i] + A43[14]*gxyy[i] + A43[15]*gyyy[i] +
               A44[15]*gxxxx[i] + A44[16]*gxxxy[i] + A44[17]*gxxyy[i] +
               A44[18]*gxyyy[i] + A44[19]*gyyyy[i];
    hvvvv[i] = A41[8]*gx[i] + A41[9]*gy[i] +
               A42[12]*gxx[i] + A42[13]*gxy[i] + A42[14]*gyy[i] +
               A43[16]*gxxx[i] + A43[17]*gxxy[i] + A43[18]*gxyy[i] + A43[19]*gyyy[i] +
               A44[20]*gxxxx[i] + A44[21]*gxxxy[i] + A44[22]*gxxyy[i] +
               A44[23]*gxyyy[i] + A44[24]*gyyyy[i];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2Derivatives4d*/

/* //////////////////////////////////////////////////////////////////// */
/* the procedures below compute the partial derivatives up to the order */
/* 4 of a composition g of a function f^-1: R^2 -> R^2 (), where f is   */
/* described with the scalar functions x, y with a function             */
/* h: R^2 -> R^d. The parameter spdimen specifies d.                    */

void pkn_Comp2iDerivatives1d ( double xu, double yu, double xv, double yv,
                               int spdimen, const double *hu, const double *hv,
                               double *gx, double *gy )
{
  void   *sp;
  double *A11, *hc;

  sp = pkv_GetScratchMemTop ();
  A11 = pkv_GetScratchMemd ( 4+2*spdimen );
  if ( !A11 )
    pkv_SignalError ( LIB_PKNUM, 7, ERRMSG_0 );
  hc = &A11[4];

  pkn_Setup2DerA11Matrixd ( xu, yu, xv, yv, A11 );
  memcpy ( hc, hu, spdimen*sizeof(double) );
  memcpy ( &hc[spdimen], hv, spdimen*sizeof(double) );

  pkn_multiGaussSolveLinEqd ( 2, A11, spdimen, spdimen, hc );
  memcpy ( gx, hc, spdimen*sizeof(double) );
  memcpy ( gy, &hc[spdimen], spdimen*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2iDerivatives1d*/

void pkn_Comp2iDerivatives2d ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               int spdimen,
                               const double *hu, const double *hv,
                               const double *huu, const double *huv,
                               const double *hvv,
                               double *gx, double *gy,
                               double *gxx, double *gxy, double *gyy )
{
  void   *sp;
  double *A21, *A22, *hc;
  int    i;

  sp = pkv_GetScratchMemTop ();
  pkn_Comp2iDerivatives1d ( xu, yu, xv, yv, spdimen, hu, hv, gx, gy );

  A21 = pkv_GetScratchMemd ( 15+3*spdimen );
  if ( !A21 )
    pkv_SignalError ( LIB_PKNUM, 8, ERRMSG_0 );
  A22 = &A21[6];  hc = &A22[9];

  pkn_Setup2DerA21Matrixd ( xuu, yuu, xuv, yuv, xvv, yvv, A21 );
  pkn_Setup2DerA22Matrixd ( xu, yu, xv, yv, A22 );

  memcpy ( hc, huu, spdimen*sizeof(double) );
  memcpy ( &hc[spdimen], huv, spdimen*sizeof(double) );
  memcpy ( &hc[2*spdimen], hvv, spdimen*sizeof(double) );
  for ( i = 0; i < spdimen; i++ ) {
    hc[i]           -= A21[0]*gx[i] + A21[1]*gy[i];
    hc[spdimen+i]   -= A21[2]*gx[i] + A21[3]*gy[i];
    hc[2*spdimen+i] -= A21[4]*gx[i] + A21[5]*gy[i];
  }

  pkn_multiGaussSolveLinEqd ( 3, A22, spdimen, spdimen, hc );
  memcpy ( gxx, hc, spdimen*sizeof(double) );
  memcpy ( gxy, &hc[spdimen], spdimen*sizeof(double) );
  memcpy ( gyy, &hc[2*spdimen], spdimen*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2iDerivatives2d*/

void pkn_Comp2iDerivatives3d ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               int spdimen,
                               const double *hu, const double *hv,
                               const double *huu, const double *huv,
                               const double *hvv,
                               const double *huuu, const double *huuv,
                               const double *huvv, const double *hvvv,
                               double *gx, double *gy,
                               double *gxx, double *gxy, double *gyy,
                               double *gxxx, double *gxxy,
                               double *gxyy, double *gyyy )
{
  void   *sp;
  double *A31, *A32, *A33, *hc;
  int    i;

  sp = pkv_GetScratchMemTop ();
  pkn_Comp2iDerivatives2d ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            spdimen, hu, hv, huu, huv, hvv,
                            gx, gy, gxx, gxy, gyy );

  A31 = pkv_GetScratchMemd ( 36+4*spdimen );
  if ( !A31 )
    pkv_SignalError ( LIB_PKNUM, 9, ERRMSG_0 );
  A32 = &A31[8];  A33 = &A32[12];  hc = &A33[16];

  pkn_Setup2DerA31Matrixd ( xuuu, yuuu, xuuv, yuuv,
                            xuvv, yuvv, xvvv, yvvv, A31 );
  pkn_Setup2DerA32Matrixd ( xu, yu, xv, yv,
                            xuu, yuu, xuv, yuv, xvv, yvv, A32 );
  pkn_Setup2DerA33Matrixd ( xu, yu, xv, yv, A33 );

  memcpy ( hc, huuu, spdimen*sizeof(double) );
  memcpy ( &hc[spdimen], huuv, spdimen*sizeof(double) );
  memcpy ( &hc[2*spdimen], huvv, spdimen*sizeof(double) );
  memcpy ( &hc[3*spdimen], hvvv, spdimen*sizeof(double) );
  for ( i = 0; i < spdimen; i++ ) {
    hc[i]           -= A31[0]*gx[i] + A31[1]*gy[i] +
                       A32[0]*gxx[i] + A32[1]*gxy[i] + A32[2]*gyy[i];
    hc[spdimen+i]   -= A31[2]*gx[i] + A31[3]*gy[i] +
                       A32[3]*gxx[i] + A32[4]*gxy[i] + A32[5]*gyy[i];
    hc[2*spdimen+i] -= A31[4]*gx[i] + A31[5]*gy[i] +
                       A32[6]*gxx[i] + A32[7]*gxy[i] + A32[8]*gyy[i];
    hc[3*spdimen+i] -= A31[6]*gx[i] + A31[7]*gy[i] +
                       A32[9]*gxx[i] + A32[10]*gxy[i] + A32[11]*gyy[i];
  }

  pkn_multiGaussSolveLinEqd ( 4, A33, spdimen, spdimen, hc );
  memcpy ( gxxx, hc, spdimen*sizeof(double) );
  memcpy ( gxxy, &hc[spdimen], spdimen*sizeof(double) );
  memcpy ( gxyy, &hc[2*spdimen], spdimen*sizeof(double) );
  memcpy ( gyyy, &hc[3*spdimen], spdimen*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2iDerivatives3d*/

void pkn_Comp2iDerivatives4d ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double xuuuu, double yuuuu, double xuuuv,
                               double yuuuv, double xuuvv, double yuuvv,
                               double xuvvv, double yuvvv,
                               double xvvvv, double yvvvv,
                               int spdimen,
                               const double *hu, const double *hv,
                               const double *huu, const double *huv,
                               const double *hvv,
                               const double *huuu, const double *huuv,
                               const double *huvv, const double *hvvv,
                               const double *huuuu, const double *huuuv,
                               const double *huuvv, const double *huvvv,
                               const double *hvvvv,
                               double *gx, double *gy,
                               double *gxx, double *gxy, double *gyy,
                               double *gxxx, double *gxxy,
                               double *gxyy, double *gyyy,
                               double *gxxxx, double *gxxxy, double *gxxyy,
                               double *gxyyy, double *gyyyy )
{
  void   *sp;
  double *A41, *A42, *A43, *A44, *hc;
  int    i;

  sp = pkv_GetScratchMemTop ();
  pkn_Comp2iDerivatives3d ( xu, yu, xv, yv,  xuu, yuu, xuv, yuv, xvv, yvv,
                            xuuu, yuuu, xuuv, yuuv, xuvv, yuvv, xvvv, yvvv,
                            spdimen, hu, hv, huu, huv, hvv,
                            huuu, huuv, huvv, hvvv,
                            gx, gy, gxx, gxy, gyy,
                            gxxx, gxxy, gxyy, gyyy );

  A41 = pkv_GetScratchMemd ( 70+5*spdimen );
  if ( !A41 )
    pkv_SignalError ( LIB_PKNUM, 10, ERRMSG_0 );
  A42 = &A41[10];  A43 = &A42[15];  A44 = &A43[20];  hc = &A44[25];

  pkn_Setup2DerA41Matrixd ( xuuuu, yuuuu, xuuuv, yuuuv, xuuvv, yuuvv,
                            xuvvv, yuvvv, xvvvv, yvvvv, A41);
  pkn_Setup2DerA42Matrixd ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            xuuu, yuuu, xuuv, yuuv, xuvv, yuvv, xvvv, yvvv,
                            A42 );
  pkn_Setup2DerA43Matrixd ( xu, yu, xv, yv,
                            xuu, yuu, xuv, yuv, xvv, yvv, A43 );
  pkn_Setup2DerA44Matrixd ( xu, yu, xv, yv, A44 );

  memcpy ( hc, huuuu, spdimen*sizeof(double) );
  memcpy ( &hc[spdimen], huuuv, spdimen*sizeof(double) );
  memcpy ( &hc[2*spdimen], huuvv, spdimen*sizeof(double) );
  memcpy ( &hc[3*spdimen], huvvv, spdimen*sizeof(double) );
  memcpy ( &hc[4*spdimen], hvvvv, spdimen*sizeof(double) );
  for ( i = 0; i < spdimen; i++ ) {
    hc[i]           -= A41[0]*gx[i] + A41[1]*gy[i] +
               A42[0]*gxx[i] + A42[1]*gxy[i] + A42[2]*gyy[i] +
               A43[0]*gxxx[i] + A43[1]*gxxy[i] + A43[2]*gxyy[i] + A43[3]*gyyy[i];
    hc[spdimen+i]   -= A41[2]*gx[i] + A41[3]*gy[i] +
               A42[3]*gxx[i] + A42[4]*gxy[i] + A42[5]*gyy[i] +
               A43[4]*gxxx[i] + A43[5]*gxxy[i] + A43[6]*gxyy[i] + A43[7]*gyyy[i];
    hc[2*spdimen+i] -= A41[4]*gx[i] + A41[5]*gy[i] +
               A42[6]*gxx[i] + A42[7]*gxy[i] + A42[8]*gyy[i] +
               A43[8]*gxxx[i] + A43[9]*gxxy[i] + A43[10]*gxyy[i] + A43[11]*gyyy[i];
    hc[3*spdimen+i] -= A41[6]*gx[i] + A41[7]*gy[i] +
               A42[9]*gxx[i] + A42[10]*gxy[i] + A42[11]*gyy[i] +
               A43[12]*gxxx[i] + A43[13]*gxxy[i] + A43[14]*gxyy[i] + A43[15]*gyyy[i];
    hc[4*spdimen+i] -= A41[8]*gx[i] + A41[9]*gy[i] +
               A42[12]*gxx[i] + A42[13]*gxy[i] + A42[14]*gyy[i] +
               A43[16]*gxxx[i] + A43[17]*gxxy[i] + A43[18]*gxyy[i] + A43[19]*gyyy[i];
  }

  pkn_multiGaussSolveLinEqd ( 5, A44, spdimen, spdimen, hc );
  memcpy ( gxxxx, hc, spdimen*sizeof(double) );
  memcpy ( gxxxy, &hc[spdimen], spdimen*sizeof(double) );
  memcpy ( gxxyy, &hc[2*spdimen], spdimen*sizeof(double) );
  memcpy ( gxyyy, &hc[3*spdimen], spdimen*sizeof(double) );
  memcpy ( gyyyy, &hc[4*spdimen], spdimen*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
} /*pkn_Comp2iDerivatives4d*/

/* //////////////////////////////////////////////////////////////////// */
/* the procedures below evaluate the partial derivatives up to the      */
/* 4 of a function f: R^2 -> R^2, given by two functions, x(u,v) and    */
/* y(u,v). The partial derivatives of f of the first order must be      */
/* linearly independent.                                                */

void pkn_f2iDerivatives1d ( double xu, double yu, double xv, double yv,
                            double *gx, double *gy )
{
  void   *sp;
  double *e1e2;

  sp = pkv_GetScratchMemTop ();
  e1e2 = pkv_GetScratchMemd ( 4 );
  if ( !e1e2 )
    pkv_SignalError ( LIB_PKNUM, 11, ERRMSG_0 );

  e1e2[0] = e1e2[3] = 1.0;
  e1e2[1] = e1e2[2] = 0.0;

  pkn_Comp2iDerivatives1d ( xu, yu, xv, yv, 2, e1e2, &e1e2[2],
                            gx, gy );

  pkv_SetScratchMemTop ( sp );
} /*pkn_f2iDerivatives1d*/

void pkn_f2iDerivatives2d ( double xu, double yu, double xv, double yv,
                            double xuu, double yuu, double xuv,
                            double yuv, double xvv, double yvv,
                            double *gx, double *gy,
                            double *gxx, double *gxy, double *gyy )
{
  void   *sp;
  double *e1e2, *zero;

  sp = pkv_GetScratchMemTop ();
  e1e2 = pkv_GetScratchMemd ( 7 );
  if ( !e1e2 )
    pkv_SignalError ( LIB_PKNUM, 12, ERRMSG_0 );
  zero = &e1e2[4];

  memset ( e1e2, 0, 7*sizeof(double) );
  e1e2[0] = e1e2[3] = 1.0;

  pkn_Comp2iDerivatives2d ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            2, e1e2, &e1e2[2], zero, zero, zero,
                            gx, gy, gxx, gxy, gyy );

  pkv_SetScratchMemTop ( sp );
} /*pkn_f2iDerivatives2d*/

void pkn_f2iDerivatives3d ( double xu, double yu, double xv, double yv,
                            double xuu, double yuu, double xuv,
                            double yuv, double xvv, double yvv,
                            double xuuu, double yuuu, double xuuv, double yuuv,
                            double xuvv, double yuvv, double xvvv, double yvvv,
                            double *gx, double *gy,
                            double *gxx, double *gxy, double *gyy,
                            double *gxxx, double *gxxy,
                            double *gxyy, double *gyyy )
{
  void   *sp;
  double *e1e2, *zero;

  sp = pkv_GetScratchMemTop ();
  e1e2 = pkv_GetScratchMemd ( 8 );
  if ( !e1e2 )
    pkv_SignalError ( LIB_PKNUM, 13, ERRMSG_0 );
  zero = &e1e2[4];

  memset ( e1e2, 0, 8*sizeof(double) );
  e1e2[0] = e1e2[3] = 1.0;

  pkn_Comp2iDerivatives3d ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            xuuu, yuuu, xuuv, yuuv, xuvv, yuvv, xvvv, yvvv,
                            2, e1e2, &e1e2[2], zero, zero, zero,
                            zero, zero, zero, zero,
                            gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy );

  pkv_SetScratchMemTop ( sp );
} /*pkn_f2iDerivatives3d*/

void pkn_f2iDerivatives4d ( double xu, double yu, double xv, double yv,
                            double xuu, double yuu, double xuv,
                            double yuv, double xvv, double yvv,
                            double xuuu, double yuuu, double xuuv, double yuuv,
                            double xuvv, double yuvv, double xvvv, double yvvv,
                            double xuuuu, double yuuuu, double xuuuv,
                            double yuuuv, double xuuvv, double yuuvv,
                            double xuvvv, double yuvvv,
                            double xvvvv, double yvvvv,
                            double *gx, double *gy,
                            double *gxx, double *gxy, double *gyy,
                            double *gxxx, double *gxxy,
                            double *gxyy, double *gyyy,
                            double *gxxxx, double *gxxxy, double *gxxyy,
                            double *gxyyy, double *gyyyy )
{
  void   *sp;
  double *e1e2, *zero;

  sp = pkv_GetScratchMemTop ();
  e1e2 = pkv_GetScratchMemd ( 9 );
  if ( !e1e2 )
    pkv_SignalError ( LIB_PKNUM, 14, ERRMSG_0 );
  zero = &e1e2[4];

  memset ( e1e2, 0, 9*sizeof(double) );
  e1e2[0] = e1e2[3] = 1.0;

  pkn_Comp2iDerivatives4d ( xu, yu, xv, yv, xuu, yuu, xuv, yuv, xvv, yvv,
                            xuuu, yuuu, xuuv, yuuv, xuvv, yuvv, xvvv, yvvv,
                            xuuuu, yuuuu, xuuuv, yuuuv, xuuvv, yuuvv,
                            xuvvv, yuvvv, xvvvv, yvvvv,
                            2, e1e2, &e1e2[2], zero, zero, zero,
                            zero, zero, zero, zero, zero, zero,
                            zero, zero, zero,
                            gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy,
                            gxxxx, gxxxy, gxxyy, gxyyy, gyyyy );

  pkv_SetScratchMemTop ( sp );
} /*pkn_f2iDerivatives4d*/

