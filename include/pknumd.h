
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpknum library of C procedures - numerical methods */

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef PKNUMD_H
#define PKNUMD_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
typedef struct bandm_profile {
  int firstnz;  /* index of the first row with a nonzero  */
                /* coefficient in the column */
  int ind;      /* index of this coefficient in the array */
} bandm_profile;

typedef struct { double x, y; } complexd;
#endif

/* ///////////////////////////////////////////////////////////////////////// */
double pkn_ScalarProductd ( int spdimen, const double *a, const double *b );
double pkn_SecondNormd ( int spdimen, const double *b );

double pkn_detd ( int n, double *a );

void pkn_AddMatrixd ( int nrows, int rowlen,
                      int inpitch1, const double *indata1,
                      int inpitch2, const double *indata2,
                      int outpitch, double *outdata );
void pkn_AddMatrixMd ( int nrows, int rowlen,
                       int inpitch1, const double *indata1,
                       int inpitch2, const double *indata2,
                       double a,
                       int outpitch, double *outdata );
void pkn_SubtractMatrixd ( int nrows, int rowlen,
                           int inpitch1, const double *indata1,
                           int inpitch2, const double *indata2,
                           int outpitch, double *outdata );
void pkn_MatrixMDifferenced ( int nrows, int rowlen,
                              int inpitch1, const double *indata1,
                              int inpitch2, const double *indata2,
                              double a,
                              int outpitch, double *outdata );
void pkn_MatrixLinCombd ( int nrows, int rowlen,
                          int inpitch1, const double *indata1,
                          double a,
                          int inpitch2, const double *indata2,
                          double b,
                          int outpitch, double *outdata );

void pkn_MultMatrixNumd ( int nrows, int rowlen,
                          int inpitch, const double *indata, 
                          double a,
                          int outpitch, double *outdata );
void pkn_MultArrayd ( int nrows, int rowlen, int pitch_a, CONST_ double *a,
                      int pitch_b, CONST_ double *b,
                      int pitch_c, double *c );

void pkn_MultMatrixd ( int nrows_a, int rowlen_a, int pitch_a, CONST_ double *a,
                       int rowlen_b, int pitch_b, CONST_ double *b,
                       int pitch_c, double *c );
void pkn_MultMatrixAddd ( int nrows_a, int rowlen_a, int pitch_a, CONST_ double *a,
                          int rowlen_b, int pitch_b, CONST_ double *b,
                          int pitch_c, double *c );
void pkn_MultMatrixSubd ( int nrows_a, int rowlen_a, int pitch_a, CONST_ double *a,
                          int rowlen_b, int pitch_b, CONST_ double *b,
                          int pitch_c, double *c );

void pkn_MultTMatrixd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                        int rowlen_b, int pitch_b, const double *b,
                        int pitch_c, double *c );
void pkn_MultTMatrixAddd ( int nrows_a, int rowlen_a, int pitch_a,
                           const double *a,
                           int rowlen_b, int pitch_b, const double *b,
                           int pitch_c, double *c );
void pkn_MultTMatrixSubd ( int nrows_a, int rowlen_a, int pitch_a,
                           const double *a,
                           int rowlen_b, int pitch_b, const double *b,
                           int pitch_c, double *c );
void pkn_MultMatrixTd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                        int nrows_b, int pitch_b, const double *b,
                        int pitch_c, double *c );
void pkn_MultMatrixTAddd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                           int nrows_b, int pitch_b, const double *b,
                           int pitch_c, double *c );
void pkn_MultMatrixTSubd ( int nrows_a, int rowlen_a, int pitch_a, const double *a,
                           int nrows_b, int pitch_b, const double *b,
                           int pitch_c, double *c );

boolean pkn_VarSignd ( int n, const double *a );
boolean pkn_OneSignChanged ( int n, const double *a, boolean *nonzero );

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_FindGivensRotationd ( double a, double b, double *c, double *s );
void pkn_FindGivensRotXid ( double a, double b, double *c, double *s, double *xi );
void pkn_FindXiGivensRotd ( double xi, double *c, double *s );
void pkn_ApplyGivensRotationd ( double c, double s, double *a, double *b );
void pkn_ApplySymGivensRotationd ( double c, double s,
                                   double *d, double *e, double *f );

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_MVectorSumd ( int m, int n, double *sum, ... );
void pkn_MVectorLinCombd ( int m, int n, double *sum, ... );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_GaussDecomposePLUQd ( int n, double *a, int *P, int *Q );
void pkn_multiSolvePLUQd ( int n, const double *lu, const int *P, const int *Q,
                           int spdimen, int pitch, double *b );
boolean pkn_multiGaussSolveLinEqd ( int n, const double *a,
                                    int spdimen, int pitch, double *b );
boolean pkn_GaussInvertMatrixd ( int n, double *a );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QRDecomposeMatrixd ( int nrows, int ncols, double *a, double *aa );
void pkn_multiReflectVectord ( int nrows, int ncols,
                               const double *a, const double *aa,
                               int spdimen, int pitch, double *b );
void pkn_multiInvReflectVectord ( int nrows, int ncols,
                                  const double *a, const double *aa,
                                  int spdimen, int pitch, double *b );
void pkn_multiMultUTVectord ( int nrows, const double *a,
                              int spdimen, int bpitch, double *b,
                              int xpitch, double *x );
void pkn_multiMultInvUTVectord ( int nrows, const double *a,
                                 int spdimen, int bpitch, double *b,
                                 int xpitch, double *x );
void pkn_multiMultTrUTVectord ( int nrows, const double *a,
                                int spdimen, int bpitch, double *b,   
                                int xpitch, double *x );
void pkn_multiMultInvTrUTVectord ( int nrows, const double *a,
                                   int spdimen, int bpitch, double *b,   
                                   int xpitch, double *x );
boolean pkn_multiSolveRLSQd ( int nrows, int ncols, double *a,
                              int spdimen, int bpitch, double *b,
                              int xpitch, double *x );
void pkn_QRGetReflectiond ( int nrows, int ncols,
                            const double *a, const double *aa,
                            int nrefl, double *w, double *gamma );


/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
void pkn_BandmFindQRMSizes ( int ncols, const bandm_profile *aprof,
                             int *qsize, int *rsize );
#endif

void pkn_BandmQRDecomposeMatrixd ( int nrows, int ncols,
                                  const bandm_profile *aprof, const double *a,
                                  bandm_profile *qprof, double *q,
                                  bandm_profile *rprof, double *r );
void pkn_multiBandmReflectVectord ( int ncols,
                                    const bandm_profile *qprof, const double *q,
                                    int spdimen, double *b );
void pkn_multiBandmInvReflectVectord ( int ncols,
                                       const bandm_profile *qprof,
                                       const double *q,
                                       int spdimen, double *b );
void pkn_multiBandmMultVectord ( int nrows, int ncols,
                                 const bandm_profile *aprof, const double *a,
                                 int spdimen, const double *x, double *y );
void pkn_multiBandmMultInvUTMVectord ( int nrows,
                                       const bandm_profile *rprof,
                                       const double *r,
                                       int spdimen, const double *x, double *y );
void pkn_multiBandmMultTrVectord ( int ncols,
                                   const bandm_profile *aprof, const double *a,
                                   int spdimen, const double *x, double *y );
void pkn_multiBandmMultInvTrUTMVectord ( int nrows,
                                         const bandm_profile *rprof,
                                         const double *r,
                                         int spdimen, const double *x, double *y );

boolean pkn_multiBandmSolveRLSQd ( int nrows, int ncols,
                                   const bandm_profile *aprof, const double *a,
                                   int nrsides, int spdimen,
                                   int bpitch, const double *b,
                                   int xpitch, double *x );
boolean pkn_multiBandmSolveDLSQd ( int nrows, int ncols,
                                   const bandm_profile *atprof, const double *at,
                                   int nrsides, int spdimen,
                                   int bpitch, const double *b,
                                   int x0pitch, const double *x0,
                                   int xpitch, double *x );
boolean pkn_multiBandmSolveCRLSQd ( int nrows, int ncols,
                                    const bandm_profile *aprof, const double *a,
                                    int nconstr, int cpitch, const double *c,
                                    int nrsides, int spdimen,
                                    int bpitch, const double *b,
                                    int dpitch, const double *d,
                                    int xpitch, double *x );

void pkn_PrintBandmd ( int ncols, const bandm_profile *aprof, const double *a );
void pkn_PrintBandmRowSumd ( int ncols, const bandm_profile *aprof, const double *a );
void pkn_PrintMatd ( int nrows, int ncols, const double *a );

#ifndef PKNUM_H
void pkn_PrintProfile ( int ncols, const bandm_profile *prof );
#endif

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
#define pkn_LowerTrMatIndex(i,j) \
  ( (i)*((i)+1)/2+(j) )
#define pkn_SymMatIndex(i,j) \
  ( (i) >= (j) ? (i)*((i)+1)/2+(j) : (j)*((j)+1)/2+(i) )
#endif

boolean pkn_CholeskyDecompd ( int n, double *a );
void pkn_SymMatrixMultd ( int n, CONST_ double *a, int spdimen,
                          int bpitch, CONST_ double *b,
                          int xpitch, double *x );
void pkn_LowerTrMatrixMultd ( int n, CONST_ double *l, int spdimen,
                              int bpitch, CONST_ double *b,
                              int xpitch, double *x );
void pkn_UpperTrMatrixMultd ( int n, CONST_ double *l, int spdimen,
                              int bpitch, CONST_ double *b,
                              int xpitch, double *x );
void pkn_LowerTrMatrixSolved ( int n, const double *l, int spdimen,
                               int bpitch, const double *b,
                               int xpitch, double *x );
void pkn_UpperTrMatrixSolved ( int n, const double *l, int spdimen,
                               int bpitch, const double *b,
                               int xpitch, double *x );

void pkn_MatrixLowerTrMultd ( int m, int n, int bpitch, CONST_ double *b,
                              CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixUpperTrMultd ( int m, int n, int bpitch, CONST_ double *b,
                              CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixLowerTrSolved ( int m, int n, int bpitch, CONST_ double *b,
                               CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixUpperTrSolved ( int m, int n, int bpitch, CONST_ double *b,
                               CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixLowerTrMultAddd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x );
void pkn_MatrixUpperTrMultAddd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x );  
boolean pkn_MatrixLowerTrSolveAddd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x );  
boolean pkn_MatrixUpperTrSolveAddd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixLowerTrMultSubd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x );  
void pkn_MatrixUpperTrMultSubd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x );  
boolean pkn_MatrixLowerTrSolveSubd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x );  
boolean pkn_MatrixUpperTrSolveSubd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x );  

void pkn_SymToFullMatrixd ( int n, const double *syma,
                            int pitch, double *fulla );
void pkn_FullToSymMatrixd ( int n, int pitch, const double *fulla,
                            double *syma );
#define pkn_FullToLTrMatrixd(n,pitch,fulla,ltra) \
  pkn_FullToSymMatrixd(n,pitch,fulla,ltra)
void pkn_LTrToFullMatrixd ( int n, const double *ltra,
                            int pitch, double *fulla );
void pkn_UTrToFullMatrixd ( int n, const double *utra,
                            int pitch, double *fulla );
void pkn_FullToUTrMatrixd ( int n, int pitch, const double *fulla,
                            double *utra );

boolean pkn_ComputeQSQTd ( int m, const double *s,   
                           int n, const double *a, const double *aa,
                           double *b );
boolean pkn_ComputeQTSQd ( int m, const double *s,
                           int n, const double *a, const double *aa,
                           double *b );

void pkn_SymMatSubAATd ( int n, double *b, int m, int pitch_a, CONST_ double *a );

void pkn_SymMatFindEigenvalueIntervald ( int n, double *a,
                                         double *lmin, double *lmax );
boolean pkn_SymMatFindEigenvaluesd ( int n, double *a, double *eigenval );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
#define pkn_LHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((i)+1),(j))-1)
#define pkn_UHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((j)+1),(i))-1)
#endif

boolean pkn_QRDecompUHessenbergd ( int n, double *ah );
void pkn_multiSolveQRUHessenbergd ( int n, double *qrh,
                                    int spdimen, int pitch, double *b );
boolean pkn_multiSolveUHessenbergLinEqd ( int n, double *ah,
                                          int spdimen, int pitch, double *b );

/* ///////////////////////////////////////////////////////////////////////// */
double pkn_Illinoisd ( double (*f)(void*,double), void *usrptr,
                       double a, double b, double eps, boolean *error );

boolean pkn_SolveSqEqd ( double p, double q, double *x1, double *x2 );

double pkn_GoldenRatd ( double (*f)(void*,double), void *usrptr,
                        double a, double b, double eps, boolean *error );


/* ///////////////////////////////////////////////////////////////////////// */
void pkn_Setup2DerA11Matrixd ( double xu, double yu, double xv, double yv,
                               double *A11 );

void pkn_Setup2DerA21Matrixd ( double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A21 );
void pkn_Setup2DerA22Matrixd ( double xu, double yu, double xv, double yv,
                               double *A22 );

void pkn_Setup2DerA31Matrixd ( double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double *A31 );
void pkn_Setup2DerA32Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A32 );
void pkn_Setup2DerA33Matrixd ( double xu, double yu, double xv, double yv,
                               double *A33 );

void pkn_Setup2DerA41Matrixd ( double xuuuu, double yuuuu, double xuuuv, double yuuuv,
                               double xuuvv, double yuuvv, double xuvvv, double yuvvv,
                               double xvvvv, double yvvvv,
                               double *A41 );
void pkn_Setup2DerA42Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double *A42 );
void pkn_Setup2DerA43Matrixd ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *A43 );
void pkn_Setup2DerA44Matrixd ( double xu, double yu, double xv, double yv,
                               double *A44 );


boolean pkn_Comp2Derivatives1d ( double xu, double yu, double xv, double yv,
                                 int spdimen, const double *gx, const double *gy,
                                 double *hu, double *hv );
boolean pkn_Comp2Derivatives2d ( double xu, double yu, double xv, double yv,
                                 double xuu, double yuu, double xuv,
                                 double yuv, double xvv, double yvv,
                                 int spdimen,
                                 const double *gx, const double *gy,
                                 const double *gxx, const double *gxy,
                                 const double *gyy,
                                 double *huu, double *huv, double *hvv );
boolean pkn_Comp2Derivatives3d ( double xu, double yu, double xv, double yv,
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
                                 double *huvv, double *hvvv );
boolean pkn_Comp2Derivatives4d ( double xu, double yu, double xv, double yv,
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
                                 double *huvvv, double *hvvvv );

boolean _pkn_Comp2iDerivatives1d ( double xu, double yu, double xv, double yv,
                                   int spdimen, const double *hu, const double *hv,
                                   double *gx, double *gy, double *workspace );
boolean pkn_Comp2iDerivatives1d ( double xu, double yu, double xv, double yv,
                                  int spdimen, const double *hu, const double *hv,
                                  double *gx, double *gy );
boolean _pkn_Comp2iDerivatives2d ( double xu, double yu, double xv, double yv,
                                   double xuu, double yuu, double xuv,
                                   double yuv, double xvv, double yvv,
                                   int spdimen,
                                   const double *hu, const double *hv,
                                   const double *huu, const double *huv,
                                   const double *hvv,
                                   double *gx, double *gy,
                                   double *gxx, double *gxy, double *gyy,
                                   double *workspace );
boolean pkn_Comp2iDerivatives2d ( double xu, double yu, double xv, double yv,
                                  double xuu, double yuu, double xuv,
                                  double yuv, double xvv, double yvv,
                                  int spdimen,
                                  const double *hu, const double *hv,
                                  const double *huu, const double *huv,
                                  const double *hvv,
                                  double *gx, double *gy,
                                  double *gxx, double *gxy, double *gyy );
boolean _pkn_Comp2iDerivatives3d ( double xu, double yu, double xv, double yv,
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
                                   double *gxyy, double *gyyy,
                                   double *workspace );
boolean pkn_Comp2iDerivatives3d ( double xu, double yu, double xv, double yv,
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
                                  double *gxyy, double *gyyy );
boolean _pkn_Comp2iDerivatives4d ( double xu, double yu, double xv, double yv,
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
                                   double *gxyyy, double *gyyyy,
                                   double *workspace );
boolean pkn_Comp2iDerivatives4d ( double xu, double yu, double xv, double yv,
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
                                  double *gxyyy, double *gyyyy );


boolean pkn_f2iDerivatives1d ( double xu, double yu, double xv, double yv,
                               double *gx, double *gy );
boolean pkn_f2iDerivatives2d ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double *gx, double *gy,
                               double *gxx, double *gxy, double *gyy );
boolean pkn_f2iDerivatives3d ( double xu, double yu, double xv, double yv,
                               double xuu, double yuu, double xuv,
                               double yuv, double xvv, double yvv,
                               double xuuu, double yuuu, double xuuv, double yuuv,
                               double xuvv, double yuvv, double xvvv, double yvvv,
                               double *gx, double *gy,
                               double *gxx, double *gxy, double *gyy,
                               double *gxxx, double *gxxy,
                               double *gxyy, double *gyyy );
boolean pkn_f2iDerivatives4d ( double xu, double yu, double xv, double yv,
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
                               double *gxyyy, double *gyyyy );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block1ArraySize ( int k, int r, int s );
int pkn_Block1FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block1FindElemPos ( int k, int r, int s, int i, int j );
#endif

boolean pkn_Block1CholeskyDecompMd ( int k, int r, int s, double *A );
void pkn_Block1LowerTrMSolved ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, double *x );
void pkn_Block1UpperTrMSolved ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, double *x );
boolean pkn_Block1SymMatrixMultd ( int k, int r, int s, CONST_ double *A,
                                   int spdimen, int xpitch, double *x,
                                   int ypitch, double *y );


/* ///////////////////////////////////////////////////////////////////////// */
int pkn_FFTPermuted ( int n, int rowlen, int pitch, complexd *a );
boolean pkn_FFTd ( int n, int rowlen, int pitch, complexd *a );
boolean pkn_InvFFTd ( int n, int rowlen, int pitch, complexd *a );


/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block2ArraySize ( int k, int r, int s, int t );
int pkn_Block2FindBlockPos ( int k, int r, int s, int t, int i, int j );
int pkn_Block2FindElemPos ( int k, int r, int s, int t, int i, int j );
#endif

boolean pkn_Block2CholeskyDecompMd ( int k, int r, int s, int t, double *A );
void pkn_Block2LowerTrMSolved ( int k, int r, int s, int t, CONST_ double *L,
                                int spdimen, int xpitch, double *x );
void pkn_Block2UpperTrMSolved ( int k, int r, int s, int t, CONST_ double *L,
                                int spdimen, int xpitch, double *x );
void pkn_Block2SymMatrixMultd ( int k, int r, int s, int t, CONST_ double *A,
                                int spdimen, int xpitch, CONST_ double *x,
                                int ypitch, double *y );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block3ArraySize ( int k, int r, int s );
int pkn_Block3FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block3FindElemPos ( int k, int r, int s, int i, int j );
#endif

boolean pkn_Block3CholeskyDecompMd ( int k, int r, int s, double *A );
void pkn_Block3LowerTrMSolved ( int k, int r, int s, CONST_ double *L,
                                int spdimen, int xpitch, double *x );
void pkn_Block3UpperTrMSolved ( int k, int r, int s, CONST_ double *L,
                                int spdimen, int xpitch, double *x );

void pkn_Block3SymMatrixMultd ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, CONST_ double *x,
                                int ypitch, double *y );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_NRBArraySize ( int n, const int *prof );
#endif

boolean pkn_NRBFindRowsd ( int n, const int *prof, CONST_ double *a,
                           double **row );
boolean pkn_NRBSymCholeskyDecompd ( int n, const int *prof, double *a,
                                    double **row, boolean *abort );
boolean pkn_NRBSymMultd ( int n, const int *prof, CONST_ double *a,
                          double **row,
                          int spdimen, int xpitch, CONST_ double *x,
                          int ypitch, double *y );
boolean pkn_NRBLowerTrMultd ( int n, const int *prof, CONST_ double *a,
                              double **row,
                              int spdimen, int xpitch, CONST_ double *x,
                              int ypitch, double *y );
boolean pkn_NRBUpperTrMultd ( int n, const int *prof, CONST_ double *a,
                              double **row,
                              int spdimen, int xpitch, CONST_ double *x,
                              int ypitch, double *y );
boolean pkn_NRBLowerTrSolved ( int n, const int *prof, CONST_ double *l,
                               double **row,
                               int spdimen, int bpitch, CONST_ double *b,
                               int xpitch, double *x );
boolean pkn_NRBUpperTrSolved ( int n, const int *prof, CONST_ double *l,
                               double **row,
                               int spdimen, int bpitch, CONST_ double *b,
                               int xpitch, double *x );

boolean pkn_NRBSymFindEigenvalueIntervald ( int n, const int *prof,
                                            double *a, double **row,
                                            double *amin, double *amax );

boolean pkn_NRBComputeQTSQd ( int n, int *prof, double *Amat, double **Arows,
                              int w, double *Bmat, double *bb,
                              int *qaprof, double **QArows );
boolean pkn_NRBComputeQSQTd ( int n, int *prof, double *Amat, double **Arows,
                              int w, double *Bmat, double *bb,
                              int *qaprof, double **QArows );

boolean pkn_NRBComputeQTSQbld ( int n, int *prof, double *Amat, double **Arows,
                                int w, double *Bmat, double *bb,
                                int *qa11prof, double **QA11rows,
                                int *qa22prof, double **QA22rows,
                                double **QA21 );
boolean pkn_NRBComputeQSQTbld ( int n, int *prof, double *Amat, double **Arows,
                                int w, double *Bmat, double *bb,
                                int *qa11prof, double **QA11rows,
                                int *qa22prof, double **QA22rows,
                                double **QA21 );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_PCGd ( int n, void *usrdata, double *b, double *x,
            boolean (*multAx)( int n, void *usrdata, const double *x, double *Ax ),
            boolean (*multQIx)( int n, void *usrdata, const double *x, double *Qix ),
            int maxit, double eps, double delta, int *itm );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_MultSPMVectord ( int nrows, int ncols, int nnz,
                             const index2 *ai, const double *ac,
                             int spdimen, const double *x,
                             double *y );
boolean pkn_MultSPMTVectord ( int nrows, int ncols, int nnz,
                              const index2 *ai, const double *ac,
                              int spdimen, const double *x,
                              double *y );

boolean pkn_MultSPsubMVectord ( int nrows, int ncols, int nnz,
                                const index3 *ai, const double *ac,
                                int spdimen, const double *x,
                                double *y );
boolean pkn_MultSPsubMTVectord ( int nrows, int ncols, int nnz,
                                 const index3 *ai, const double *ac,
                                 int spdimen, const double *x,
                                 double *y );

/* fast sparse matrix multiplication - uses data produced by */
/* one of the pkn_SPMFind*nnz* or pkn_SPsubMFind*nnz* procedures */
void pkn_SPMFastMultMMd ( double *ac, double *bc,
                          int nnzab, int *abpos, index2 *aikbkj,
                          double *abc );

/* slower (but still decent) matrix multiplication procedures */
boolean pkn_SPMmultMMCd ( int nra, int nca, int ncb,
                          unsigned int nnza, index2 *ai, double *ac,
                          unsigned int *apermut, int *acols, boolean ca,
                          unsigned int nnzb, index2 *bi, double *bc,
                          unsigned int *bpermut, int *bcols, boolean cb,
                          index2 *abi, double *abc );
boolean pkn_SPMmultMMTCd ( int nra, int nca, int nrb,
                           unsigned int nnza, index2 *ai, double *ac,
                           unsigned int *apermut, int *acols, boolean ca,
                           unsigned int nnzb, index2 *bi, double *bc,
                           unsigned int *bpermut, int *brows, boolean rb,
                           index2 *abi, double *abc );
boolean pkn_SPMmultMTMCd ( int nra, int nca, int ncb,
                           unsigned int nnza, index2 *ai, double *ac,
                           unsigned int *apermut, int *arows, boolean ra,
                           unsigned int nnzb, index2 *bi, double *bc,
                           unsigned int *bpermut, int *bcols, boolean ba,
                           index2 *abi, double *abc );

boolean pkn_SPsubMmultMMCd ( int nra, int nca, int ncb,
                             unsigned int nnza, index3 *ai, double *ac,
                             unsigned int *apermut, int *acols, boolean ca,
                             unsigned int nnzb, index3 *bi, double *bc,
                             unsigned int *bpermut, int *bcols, boolean cb,
                             index2 *abi, double *abc );
boolean pkn_SPsubMmultMMTCd ( int nra, int nca, int nrb,
                              unsigned int nnza, index3 *ai, double *ac,
                              unsigned int *apermut, int *acols, boolean ca,
                              unsigned int nnzb, index3 *bi, double *bc,
                              unsigned int *bpermut, int *brows, boolean rb,
                              index2 *abi, double *abc );
boolean pkn_SPsubMmultMTMCd ( int nra, int nca, int ncb,
                              unsigned int nnza, index3 *ai, double *ac,
                              unsigned int *apermut, int *arows, boolean ra,
                              unsigned int nnzb, index3 *bi, double *bc,
                              unsigned int *bpermut, int *bcols, boolean cb,
                              index2 *abi, double *abc );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QuadRectanglesd ( double a, double b, int n,
                              double *qknots, double *qcoeff );
boolean pkn_QuadSimpsond ( double a, double b, int n,
                           double *qknots, double *qcoeff );

boolean pkn_QuadGaussLegendre4d ( double a, double b, int n,
                                  double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre6d ( double a, double b, int n,
                                  double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre8d ( double a, double b, int n,
                                  double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre10d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre12d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre14d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre16d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre18d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );
boolean pkn_QuadGaussLegendre20d ( double a, double b, int n,
                                   double *qknots, double *qcoeff );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKN_LMT_CONTINUE
/* return values of pkn_NLMIterd */
#define PKN_LMT_ERROR           -1
#define PKN_LMT_CONTINUE_LM      0
#define PKN_LMT_CONTINUE_LM_P    1
#define PKN_LMT_CONTINUE_N       2
#define PKN_LMT_FOUND_MINIMUM    3
#define PKN_LMT_FOUND_ZEROGRAD   4
#define PKN_LMT_FOUND_ZEROGRAD_P 5
#define PKN_LMT_FOUND_BARRIER    6
#define PKN_LMT_CROSSED_LIMIT    7
#define PKN_LMT_NO_PROGRESS      8
#endif

#ifndef PKN_SD_CONTINUE
/* return values of pkn_SDIterd */
#define PKN_SD_ERROR            -1
#define PKN_SD_CONTINUE          0
#define PKN_SD_FOUND_ZEROGRAD    1
#define PKN_SD_FOUND_BARRIER     2
#define PKN_SD_CROSSED_LIMIT     3
#define PKN_SD_NO_PROGRESS       4
#endif

typedef boolean (*pkn_NLMTevalfuncd)( int n, void *usrdata,
                                      double *x, double *f );
typedef boolean (*pkn_NLMTevalfuncgd)( int n, void *usrdata,
                                       double *x, double *f, double *g );
typedef boolean (*pkn_NLMTevalfuncghd)( int n, void *usrdata,
                                        double *x, double *f, double *g, double *h );
typedef boolean (*pkn_NLMTtransfuncd)( int n, void *usrdata, double *x );
typedef boolean (*pkn_NLMTtunnelfuncd)( int n, void *usrdata,
                                        double *x0, double *x1, boolean *went_out );

int pkn_NLMIterd ( int n, void *usrdata, double *x,
                   pkn_NLMTevalfuncd funcf, pkn_NLMTevalfuncgd funcfg,
                   pkn_NLMTevalfuncghd funcfgh,
                   pkn_NLMTtransfuncd trans, pkn_NLMTtunnelfuncd tunnel,
                   double lowerbound, double eps, double delta,
                   double *nu );

int pkn_SDIterd ( int n, void *usrdata, double *x,
                  pkn_NLMTevalfuncd funcf, pkn_NLMTevalfuncgd funcfg,
                  pkn_NLMTtransfuncd trans, pkn_NLMTtunnelfuncd tunnel,
                  double lowerbound, double eps, double delta,
                  double *nu );

boolean _pkn_DivideIntervald ( double *ga, double *gc, double *gd, double *gb,
                               double *fa, double *fc, double *fd, double *fb );

#ifdef __cplusplus
}
#endif

#endif /* PKNUMD_H*/

