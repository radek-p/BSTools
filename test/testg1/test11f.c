
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "psout.h"
#include "multibs.h"  
#include "eg1holef.h"

#include "datagenf.h"
#include "drawitf.h"
#include "writdatf.h"

#define _TESTMATRIX

#define HOLE_K 5

#define NK 2
#define M1 2
#define M2 2

#define QUAD_FACTOR 10

char fn1[] = "g1q2_splmatrixf.ps";
char fn2[] = "g1q2splsurff.ps";

GHoleDomainf *domain;

int nfinalp;
point3f FinalCPoints[HOLE_K][121];
int nfunc_a, nfunc_b, nfunc_c, nfunc_d;

int     nfpatches = 0;
int     fpdeg, fplkn;
float   *fpknots;
point3f *fpcp[HOLE_K];

float iparams[3] = {0.5,0.5,0.5};

float q2const = 1.0/*20.0*/;

/* ///////////////////////////////////////////////////////////////////////// */
void PaintMatrix ( int k, int r, int s, float *A, float *B )
{
  void  *sp;
  int   n, m, i, j, p;
  byte  *buf;
  char  ss[40];
  float scale;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  m = 6*HOLE_K+1;
  buf = pkv_GetScratchMem ( max(n,m) );
  if ( buf ) {
    scale = 2000.0/(float)n;
    sprintf ( ss, "100 100 translate %6.2f dup scale", scale );
    ps_GSave ();
    ps_Write_Command ( ss );
    ps_Init_BitmapP ( n, n, 0, 0 );
    for ( i = 0; i < n; i++ ) {
      memset ( buf, 0xDF, n );
      for ( j = 0; j < n; j++ ) {
        p = pkn_Block3FindElemPos ( k, r, s, i, j );
        if ( p >= 0 ) {
          if ( A[p] )
            buf[j] = 0;
        }
        else buf[j] = 0xAF;
      }
      ps_Out_LineP ( buf );
    }
    ps_Init_BitmapP ( m, n, n+10, 0 );
    for ( i = p = 0; i < n; i++ ) {
      memset ( buf, 0xDF, m );
      for ( j = 0; j < m; j++, p++ )
        if ( B[p] ) buf[j] = 0;
      ps_Out_LineP ( buf );
    }
    ps_GRestore ();
  }
  pkv_SetScratchMemTop ( sp );
} /*PaintMatrix*/

void PaintMatrixError ( int k, int r, int s, float *A, float *a,
                        float *B, float *b, float eps )
{
  void  *sp;
  int   n, m, i, j, p;
  byte  *buf;
  char  ss[40];
  float scale;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  m = 6*HOLE_K+1;
  buf = pkv_GetScratchMem ( 3*max(n,m) );
  if ( buf ) {
    scale = 2000.0/(float)n;
    sprintf ( ss, "100 100 translate %6.2f dup scale", scale );
    ps_GSave ();
    ps_Write_Command ( ss );
    ps_Init_BitmapRGBP ( n, n, 0, 0 );
    for ( i = 0; i < n; i++ ) {
      memset ( buf, 0xDF, 3*n );
      for ( j = 0; j < 3*n; j += 3 ) {
        p = pkn_Block3FindElemPos ( k, r, s, i, j/3 );
        if ( p >= 0 ) {
          if ( A[p] || a[p] ) {
            if ( fabs(A[p]-a[p]) < eps )
              buf[j] = buf[j+1] = buf[j+2] = 0;
            else if ( !A[p] )
              buf[j] = 255, buf[j+1] = buf[j+2] = 0;
            else
              buf[j] = buf[j+1] = 0, buf[j+2] = 255;
          }
        }
        else buf[j] = buf[j+1] = buf[j+2] = 0xAF;
      }
      ps_Out_LineRGBP ( buf );
    }
    ps_Init_BitmapRGBP ( m, n, n+10, 0 );
    for ( i = p = 0; i < n; i++ ) {
      memset ( buf, 0xDF, 3*m );
      for ( j = 0; j < 3*m; j += 3, p++ ) {
        if ( B[p] || b[p] ) {
          if ( fabs(B[p]-b[p]) < eps )
            buf[j] = buf[j+1] = buf[j+2] = 0;
          else if ( !B[p] )
            buf[j] = 255, buf[j+1] = buf[j+2] = 0;
          else
            buf[j] = buf[j+1] = 0, buf[j+2] = 255;
        }
      }
      ps_Out_LineRGBP ( buf );
    }
    ps_GRestore ();
  }
  pkv_SetScratchMemTop ( sp );
} /*PaintMatrixError*/

void FindCondNumber ( int k, int r, int s, const float *a )
{
  void  *sp;
  float *b, *c, *x, *y;
  float lambda_1, lambda_n, pr, yl, ca;
  int   i, j, n, size;

  sp = pkv_GetScratchMemTop ();
  size = pkn_Block3ArraySize ( k, r, s );
  n = k*r+s;

  /* power method of finding the greatest eigenvalue */
  j = -1;
  yl = 0.0;
  for ( i = 0; i < size; i++ )
    if ( fabs(a[i]) > yl ) {
      yl = max ( yl, fabs(a[i]) );
      j = i;
    }
  printf ( "a_max = %g\n", yl );

  b = pkv_GetScratchMemf ( size );
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  if ( !b || !x || !y )
    exit ( 1 );
  c = &b[n*n];
  memcpy ( b, a, size*sizeof(float) );

  for ( x[0] = 1.0, i = 1;  i < n;  i++ )
    x[i] = 1.01*x[i-1];
  yl = pkn_SecondNormf ( n, x );
  pkn_MultMatrixNumf ( 1, n, 0, x, 1.0/yl, 0, x );

  for ( i = 0; i < 100; i++ ) {
    pkn_Block3SymMatrixMultf ( k, r, s, b, 1, 1, x, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
/*
  printf ( "ca = %f, yl = %g\n", ca, yl );
*/
  lambda_1 = yl;

  pkn_Block3CholeskyDecompMf ( k, r, s, b );
  for ( x[0] = 1.0, i = 1;  i < n;  i++ )
    x[i] = 1.01*x[i-1];
  yl = pkn_SecondNormf ( n, x );
  pkn_MultMatrixNumf ( 1, n, 0, x, 1.0/yl, 0, x );
  for ( i = 0; i < 100; i++ ) {
    memcpy ( y, x, n*sizeof(float) );
    pkn_Block3LowerTrMSolvef ( k, r, s, b, 1, 1, y );
    pkn_Block3UpperTrMSolvef ( k, r, s, b, 1, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
/*
  printf ( "ca = %f, yl = %g\n", ca, yl );
*/
  lambda_n = 1.0/yl;
  printf ( "l_1 = %g, l_n = %g, cond A = %g\n",
           lambda_1, lambda_n, lambda_1/lambda_n );

  pkv_SetScratchMemTop ( sp );
} /*FindCondNumber*/

void DrawMatrices ( int k, int r, int s, float *amat, float *bmat )
{
  PaintMatrix ( k, r, s, amat, bmat );
  FindCondNumber ( k, r, s, amat );
} /*DrawMatrices*/

void MakePicture1 ( void )
{
  ps_OpenFile ( fn1, 600 );
  g1h_Q2DrawSplMatricesf ( domain, DrawMatrices );
  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( float w, float h, float x, float y )
{
/*
  vector3f v;
*/
  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
/*
  SetVector3f ( &v, 0.0, 0.0, -25.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, 0.25*PI );
  CameraZoomf ( &CPos, 5.0 );
*/
  SetPoint3f ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 4.0;
  CameraSetMappingf ( &CPos );
} /*SetupCamera*/

void GetHoleSurrndPatch ( int i, int j, point3f *dombcp )
{
  int     k;   
  int     *ind;
  point3f *q; 
  void    *sp;
  float   *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q )
    exit ( 1 );

  gh_GetBspInd ( HOLE_K, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];

  ukn = &knots[11*((i+HOLE_K-1) % HOLE_K)+3];
  vkn = &knots[11*i+j];
  mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)dombcp );

  pkv_SetScratchMemTop ( sp );
} /*GetHoleSurrndPatch*/

void DrawBezPatches ( void )
{
  int     i, j;  
  point3f cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3f ( 3, 3, cp, 8, 8 );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ( void )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfpatches; i++ )
    DrawBSPatch3f ( fpdeg, fplkn, fpknots, fpdeg, fplkn, fpknots,
                    fpcp[i], 4, 4 );
} /*DrawFinalPatches*/

void MakePicture2 ( void )
{
  ps_OpenFile ( fn2, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn2 );
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ( void )
{
  MakePicture1 ();
  MakePicture2 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
int m_k, m_r, m_s, m_size, nf_a, nf_b, nf_c;
float *_amat = NULL, *_bmat = NULL;
float *Amat = NULL, *Bmat = NULL;
int ndip, ndsp, nf1, nf2;
point2f dipatch[HOLE_K][36];
point2f dspatch[3*HOLE_K][16];
float func1p[HOLE_K][400];
float func2p[4*HOLE_K][400];
int   bflkn1, bflkn2;
float bfkn1[28], bfkn2[28];

void SaveMatrices ( int k, int r, int s, float *amat, float *bmat )
{
  m_size = pkn_Block3ArraySize ( k, r, s );
  m_k = k;  m_r = r;  m_s = s;
  nf_a = k*r+s;
  nf_b = 6*HOLE_K+1;
  nf_c = (k+1)*r;
  _amat = malloc ( m_size*sizeof(float) );
  _bmat = malloc ( nf_a*nf_b*sizeof(float) );
  if ( !_amat || !_bmat )
    exit ( 1 );
  memcpy ( _amat, amat, m_size*sizeof(float) );
  memcpy ( _bmat, bmat, nf_a*nf_b*sizeof(float) );
} /*SaveMatrices*/

void OutDiPatches ( int n, int m, const point2f *cp )
{
  memcpy ( dipatch[ndip], cp, 36*sizeof(point2f) );
  ndip ++;
} /*OutDiPatches*/

void OutDomSPatches ( int n, int m, const point2f *cp )
{
  memcpy ( dspatch[ndsp], cp, 16*sizeof(point2f) );
  ndsp ++;
} /*OutDomSPatches*/

void OutBasf1 ( int n, int m, const point3f *f )
{
  float zeroone[2] = {0.0,1.0};

  pkv_Selectf ( 36, 1, 3, 1, &f[0].z, func1p[nf1] );
  mbs_SetKnotPatternf ( 1, zeroone, n+1, &bflkn1, bfkn1 );
  nf1 ++;
} /*OutBasf1*/

void OutBasf2 ( int n, int m, const point3f *f )
{
  float zeroone[2] = {0.0,1.0};

  pkv_Selectf ( 36, 1, 3, 1, &f[0].z, func2p[nf2] );
  mbs_SetKnotPatternf ( 1, zeroone, n+1, &bflkn2, bfkn2 );
  nf2 ++;
} /*OutBasf2*/

void OutBasf1s ( int n, int lknu, const float *knu,
                 int m, int lknv, const float *knv, const point3f *f )
{
  int ncp;

  ncp = lknu-n;  ncp *= ncp;
  pkv_Selectf ( ncp, 1, 3, 1, &f[0].z, func1p[nf1] );
  bflkn1 = lknu;
  memcpy ( bfkn1, knu, (lknu+1)*sizeof(float) );
  nf1 ++;
} /*OutBasf1s*/

void OutBasf2s ( int n, int lknu, const float *knu,
                 int m, int lknv, const float *knv, const point3f *f )
{
  int ncp;

  ncp = lknu-n;  ncp *= ncp;
  pkv_Selectf ( ncp, 1, 3, 1, &f[0].z, func2p[nf2] );
  bflkn2 = lknu;
  memcpy ( bfkn2, knu, (lknu+1)*sizeof(float) );
  nf2 ++;
} /*OutBasf2s*/

boolean IsNonZero ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const float *f )
{
  int ncp, i;

  ncp = (lknu-n)*(lknv-m);
  for ( i = 0; i < ncp; i++ )
    if ( f[i] )
      return true;
  return false;
} /*IsNonZero*/

float CompIntegral1 ( void )
{
  void     *sp;
  double   s, ss, sss;
  int      i, j, k, d1, d2, nquad;
  vector2f di, diu, div, diuu, diuv, divv, diuuu, diuuv, diuvv, divvv;
  float    f1[4], f1v[4], f1vv[4], f1vvv[4];
  float    f2[4], f2v[4], f2vv[4], f2vvv[4];
  float    f1x, f1y, f1xx, f1xy, f1yy, f1xxx, f1xxy, f1xyy, f1yyy;
  float    f2x, f2y, f2xx, f2xy, f2yy, f2xxx, f2xxy, f2xyy, f2yyy;
  float    *ap1, *ap2;
  float    u, v, jac;
  vector2f lgr1, lgr2;

  sp = pkv_GetScratchMemTop ();
  ap1 = pkv_GetScratchMemf ( 4*(d1 = bflkn1-G1H_FINALDEG) );
  ap2 = pkv_GetScratchMemf ( 4*(d2 = bflkn2-G1H_FINALDEG) );
  if ( !ap1 || !ap2 )
    exit ( 1 );
  nquad = QUAD_FACTOR*(NK+1);
  sss = 0.0;
  for ( i = 0; i < HOLE_K; i++ ) {
    if ( IsNonZero ( G1H_FINALDEG, bflkn1, bfkn1, G1H_FINALDEG, bflkn1, bfkn1,
                     func1p[i] ) &&
         IsNonZero ( G1H_FINALDEG, bflkn2, bfkn2, G1H_FINALDEG, bflkn2, bfkn2,
                     func2p[3*HOLE_K+i] ) ) {
      ss = 0.0;
      for ( j = 0; j < nquad; j++ ) {
        u = (float)(j+j+1)/(float)(2*nquad);
        mbs_multideBoorDer3f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
               u, ap1, &ap1[d1], &ap1[2*d1], &ap1[3*d1] );
        mbs_multideBoorDer3f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
               u, ap2, &ap2[d2], &ap2[2*d2], &ap2[3*d2] );
        s = 0.0;
        for ( k = 0; k < nquad; k++ ) {
          v = (float)(k+k+1)/(float)(2*nquad);
          mbs_BCHornerDer3Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                 u, v, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x,
                 &diuuu.x, &diuuv.x, &diuvv.x, &divvv.x );
          jac = fabs ( diu.x*div.y-div.x*diu.y );
          mbs_multideBoorDer3f ( G1H_FINALDEG, bflkn1, bfkn1, 4, 1, d1, ap1,
                 v, f1, f1v, f1vv, f1vvv );
          mbs_multideBoorDer3f ( G1H_FINALDEG, bflkn2, bfkn2, 4, 1, d2, ap2,
                 v, f2, f2v, f2vv, f2vvv );
          pkn_Comp2iDerivatives3f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y, diuuu.x, diuuu.y, diuuv.x, diuuv.y,
                 diuvv.x, diuvv.y, divvv.x, divvv.y,
                 1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                 &f1[3], &f1v[2], &f1vv[1], &f1vvv[0],
                 &f1x, &f1y, &f1xx, &f1xy, &f1yy, &f1xxx, &f1xxy, &f1xyy, &f1yyy );
          lgr1.x = f1xxx + f1xyy;
          lgr1.y = f1xxy + f1yyy;
          pkn_Comp2iDerivatives3f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y, diuuu.x, diuuu.y, diuuv.x, diuuv.y,
                 diuvv.x, diuvv.y, divvv.x, divvv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2[3], &f2v[2], &f2vv[1], &f2vvv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy, &f2xxx, &f2xxy, &f2xyy, &f2yyy );
          lgr2.x = f2xxx + f2xyy;
          lgr2.y = f2xxy + f2yyy;
          s += (lgr1.x*lgr2.x+lgr1.y*lgr2.y)*jac;
        }
        ss += s;
      }
      sss += ss;
    }
  }
  pkv_SetScratchMemTop ( sp );
  return sss/(float)(nquad*nquad);
} /*CompIntegral1*/

float CompIntegral2 ( void )
{
  void     *sp;
  int      nquad, i, j, k, l, ii, jj, d1, d2;
  double   s, ss;
  float    t, u, jac, lj1, lj2;
  float    *ap1, *ap2;
  float    f1[3], f1v[3], f1vv[3], f2[3], f2v[3], f2vv[3];
  float    f1x, f1y, f1xx, f1xy, f1yy, f2x, f2y, f2xx, f2xy, f2yy;
  vector2f di, diu, div, diuu, diuv, divv;

  sp = pkv_GetScratchMemTop ();
  ap1 = pkv_GetScratchMemf ( 3*(d1 = bflkn1-G1H_FINALDEG) );
  ap2 = pkv_GetScratchMemf ( 3*(d2 = bflkn2-G1H_FINALDEG) );
  if ( !ap1 || !ap2 )
    exit ( 1 );
  nquad = QUAD_FACTOR*(NK+1);
  ss = 0.0;
  for ( i = 0; i < HOLE_K; i++ ) {
    ii = (i+HOLE_K-1) % HOLE_K;
    jj = (i+1) % HOLE_K;
    s = 0.0;
          /* integrate along v == 0 */
    for ( j = 0; j < nquad; j++ ) {
      t = (float)(j+j+1)/(float)(2*nquad);
      mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                  t, 0.0, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      jac = sqrt ( diu.x*diu.x+diu.y*diu.y );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
                             t, ap1, &ap1[d1], &ap1[2*d1] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                             0.0, f1, f1v, f1vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                 &f1x, &f1y, &f1xx, &f1xy, &f1yy );
      lj1 = f1xx+f1yy;
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                             t, ap2, &ap2[d2], &ap2[2*d2] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 3, 1, d2, ap2,
                             0.0, f2, f2v, f2vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 = f2xx+f2yy;
      mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[ii],
                  0.0, t, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[ii],
                             0.0, ap1, &ap1[d1], &ap1[2*d1] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                             t, f1, f1v, f1vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                 &f1x, &f1y, &f1xx, &f1xy, &f1yy );
      lj1 -= f1xx+f1yy;
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+ii],
                             0.0, ap2, &ap2[d2], &ap2[2*d2] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 3, 1, d2, ap2,
                             t, f2, f2v, f2vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 -= f2xx+f2yy;
      s += lj1*lj2*jac;
    }
          /* integrate along v == 1 */
    for ( j = 0; j < nquad; j++ ) {
      t = (float)(j+j+1)/(float)(2*nquad);
      mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                  t, 1.0, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      jac = sqrt ( diu.x*diu.x+diu.y*diu.y );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
                             t, ap1, &ap1[d1], &ap1[2*d1] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                             1.0, f1, f1v, f1vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                 &f1x, &f1y, &f1xx, &f1xy, &f1yy );
      lj1 = -(f1xx+f1yy);
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                             t, ap2, &ap2[d2], &ap2[2*d2] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 3, 1, d2, ap2,
                             1.0, f2, f2v, f2vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 = f2xx+f2yy;
      mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)dspatch[3*jj+1],
                     0.0, t, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      mbs_BCHornerDer2Pf ( 3, 3, 1, func2p[3*jj+1], 0.0, t,
                   &f2[0], &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0] );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 = f2xx+f2yy-lj2;
      s += lj1*lj2*jac;
    }
          /* integrate along u == 1 */
    for ( j = 0; j < nquad; j++ ) {
      t = (float)(j+j+1)/(float)(2*nquad);
      mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                  1.0, t, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      jac = sqrt ( div.x*div.x+div.y*div.y );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
                             1.0, ap1, &ap1[d1], &ap1[2*d1] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                             t, f1, f1v, f1vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                 &f1x, &f1y, &f1xx, &f1xy, &f1yy );
      lj1 = -(f1xx+f1yy);
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                             1.0, ap2, &ap2[d2], &ap2[2*d2] );
      mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 3, 1, d2, ap2,
                             t, f2, f2v, f2vv );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 = f2xx+f2yy;
      mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)dspatch[3*i+2],
                     0.0, t, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
      mbs_BCHornerDer2Pf ( 3, 3, 1, func2p[3*i+2], 0.0, t,
                   &f2[0], &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0] );
      pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                 diuv.x, diuv.y, divv.x, divv.y,
                 1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                 &f2x, &f2y, &f2xx, &f2xy, &f2yy );
      lj2 = f2xx+f2yy-lj2;
      s += lj1*lj2*jac;
    }
    if ( (M1 > 1 || M2 > 3) &&
         bflkn1 > 2*G1H_FINALDEG+1 && bflkn2 > 2*G1H_FINALDEG+1 ) {
          /* integrate along "inner" curves v == const */
      for ( k = 0; k < NK; k++ ) {
        u = (float)(k+1)/(float)(NK+1);
        for ( j = 0; j < nquad; j++ ) {
          t = (float)(j+j+1)/(float)(2*nquad);
          mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                      t, u, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
          jac = sqrt ( diu.x*diu.x+diu.y*diu.y );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
                                 t, ap1, &ap1[d1], &ap1[2*d1] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                                 u, f1, f1v, f1vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                     &f1x, &f1y, &f1xx, &f1xy, &f1yy );
          lj1 = f1xx+f1yy;
          l = mbs_FindKnotIntervalf ( G1H_FINALDEG, bflkn1, bfkn1, u-0.01, &l );
          mbs_multideBoorDer2f ( G1H_FINALDEG, l+G1H_FINALDEG+1, bfkn1,
                                 3, 1, d1, ap1,
                                 u, f1, f1v, f1vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                     &f1x, &f1y, &f1xx, &f1xy, &f1yy );
          lj1 -= f1xx+f1yy;
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                                 t, ap2, &ap2[d2], &ap2[2*d2] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 3, 1, d2, ap2,
                                 u, f2, f2v, f2vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                     &f2x, &f2y, &f2xx, &f2xy, &f2yy );
          lj2 = f2xx+f2yy;
          l = mbs_FindKnotIntervalf ( G1H_FINALDEG, bflkn2, bfkn2, u-0.01, &l );
          mbs_multideBoorDer2f ( G1H_FINALDEG, l+G1H_FINALDEG+1, bfkn2,
                                 3, 1, d2, ap2,
                                 u, f2, f2v, f2vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                     &f2x, &f2y, &f2xx, &f2xy, &f2yy );
          lj2 -= f2xx+f2yy;
          s += lj1*lj2*jac;
        }
      }
          /* integrate along "inner" curves u == const */
      for ( k = 0; k < NK; k++ ) {
        u = (float)(k+1)/(float)(NK+1);
        for ( j = 0; j < nquad; j++ ) {
          t = (float)(j+j+1)/(float)(2*nquad);
          mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)dipatch[i],
                      u, t, &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x );
          jac = sqrt ( div.x*div.x+div.y*div.y );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 1, d1, 0, func1p[i],
                                 u, ap1, &ap1[d1], &ap1[2*d1] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1, 3, 1, d1, ap1,
                                 t, f1, f1v, f1vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                     &f1x, &f1y, &f1xx, &f1xy, &f1yy );
          lj1 = f1xx+f1yy;
          l = mbs_FindKnotIntervalf ( G1H_FINALDEG, bflkn1, bfkn1, u-0.01, &l );
          mbs_multideBoorDer2f ( G1H_FINALDEG, l+G1H_FINALDEG+1, bfkn1, 1, d1, 0, func1p[i],
                                 u, ap1, &ap1[d1], &ap1[2*d1] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn1,
                                 3, 1, d1, ap1,
                                 t, f1, f1v, f1vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f1[1], &f1v[0], &f1[2], &f1v[1], &f1vv[0],
                     &f1x, &f1y, &f1xx, &f1xy, &f1yy );
          lj1 -= f1xx+f1yy;
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                                 u, ap2, &ap2[d2], &ap2[2*d2] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn1, bfkn2, 3, 1, d2, ap2,
                                 t, f2, f2v, f2vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                     &f2x, &f2y, &f2xx, &f2xy, &f2yy );
          lj2 = f2xx+f2yy;
          l = mbs_FindKnotIntervalf ( G1H_FINALDEG, bflkn2, bfkn2, u-0.01, &l );
          mbs_multideBoorDer2f ( G1H_FINALDEG, l+G1H_FINALDEG+1, bfkn2, 1, d2, 0, func2p[3*HOLE_K+i],
                                 u, ap2, &ap2[d2], &ap2[2*d2] );
          mbs_multideBoorDer2f ( G1H_FINALDEG, bflkn2, bfkn2,
                                 3, 1, d2, ap2,
                                 t, f2, f2v, f2vv );
          pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                     diuv.x, diuv.y, divv.x, divv.y,
                     1, &f2[1], &f2v[0], &f2[2], &f2v[1], &f2vv[0],
                     &f2x, &f2y, &f2xx, &f2xy, &f2yy );
          lj2 -= f2xx+f2yy;
          s += lj1*lj2*jac;
        }
      }
    }
    ss += s;
  }
  pkv_SetScratchMemTop ( sp );
  return q2const*ss/(float)nquad;
} /*CompIntegral2*/

float CompCoeff ( int fni, int fnj )
{
  float int1, int2;

  nf1 = 0;
  g1h_DrawSplBasFunctionf ( domain, fni, OutBasf1s );
  if ( fnj >= nfunc_a && fnj < nfunc_a+nfunc_b ) {
    nf2 = 0;
    g1h_DrawBasBFunctionf ( domain, fnj-nfunc_a, OutBasf2 );
  }
  else
    memset ( func2p, 0, sizeof(func2p) );
  nf2 = 3*HOLE_K;
  g1h_DrawSplBasFunctionf ( domain, fnj, OutBasf2s );
  int1 = CompIntegral1 ();
  int2 = CompIntegral2 ();
  return int1 + int2;
} /*CompCoeff*/

void ComputeMatrices ( void )
{
  int   i, j, fni, fnj, k;

  ndip = ndsp = 0;
  g1h_DrawDiPatchesf ( domain, OutDiPatches );
  gh_DrawDomSurrndPatchesf ( domain, OutDomSPatches );

  memset ( Amat, 0, m_size*sizeof(float) );
  memset ( Bmat, 0, nf_a*nf_b*sizeof(float) );

CompCoeff ( 0, 0 );
/*exit ( 0 ); */

  for ( i = 0; i < nfunc_a+nfunc_c+nfunc_d; i++ ) {
    if ( i < nfunc_c+nfunc_d ) fni = nfunc_a+nfunc_b+i;
    else fni = i-nfunc_c-nfunc_d;
    printf ( "%4d %4d\b\b\b\b\b\b\b\b\b", i, fni );
    nf1 = 0;
    g1h_DrawSplBasFunctionf ( domain, fni, OutBasf1s );
    for ( j = 0; j <= i; j++ ) {
      if ( j < nfunc_c+nfunc_d ) fnj = nfunc_a+nfunc_b+j;
      else fnj = j-nfunc_c-nfunc_d;
      k = pkn_Block3FindElemPos ( HOLE_K-1, nfunc_c/HOLE_K,
                                  nfunc_c/HOLE_K+nfunc_d+nfunc_a, i, j );
      if ( k >= 0 ) {
        memset ( func2p, 0, sizeof(func2p) );
        nf2 = 3*HOLE_K;
        g1h_DrawSplBasFunctionf ( domain, fnj, OutBasf2s );
        Amat[k] = CompIntegral1 () + CompIntegral2 ();
      }
    }
    for ( j = 0; j < nfunc_b; j++ ) {
      nf2 = 0;
      g1h_DrawBasBFunctionf ( domain, j, OutBasf2 );
      nf2 = 3*HOLE_K;
      g1h_DrawSplBasFunctionf ( domain, nfunc_a+j, OutBasf2s );
      Bmat[nfunc_b*i+j] = CompIntegral1 () + CompIntegral2 ();
    }
  }
  printf ( "\n" );
} /*ComputeMatrices*/

void CompareMatrices ( void )
{
  FILE  *f;
  int   i, j, k, i0, j0;
  float e, maxe;

  f = fopen ( "splmatrf.txt", "w+" );
  maxe = 0.0;  i0 = j0 = 0;
  for ( i = 0; i < nfunc_a+nfunc_c+nfunc_d; i++ )
    for ( j = 0; j <= i; j++ ) {
      k = pkn_Block3FindElemPos ( HOLE_K-1, nfunc_c/HOLE_K,
                                  nfunc_c/HOLE_K+nfunc_d+nfunc_a, i, j );
      if ( k >= 0 )
        fprintf ( f, "%3d %3d %4d: %10.4f %10.4f %10.4f\n",
                  i, j, k, _amat[k], Amat[k], e = _amat[k]-Amat[k] );
      if ( (e = fabs(e)) > maxe ) {
        maxe = e;
        i0 = i;
        j0 = j;
      }
    }
  k = pkn_Block3FindElemPos ( HOLE_K-1, nfunc_c/HOLE_K,
                              nfunc_c/HOLE_K+nfunc_d+nfunc_a, i0, j0 );
  printf ( "%3d %3d %4d: %10.4f %10.4f %10.4f\n",
           i0, j0, k, _amat[k], Amat[k], _amat[k]-Amat[k] );
  
  fprintf ( f, "\n" );
  maxe = 0.0;  i0 = 0;
  for ( i = 0; i < nf_a*nf_b; i++ ) {
    fprintf ( f, "%4d: %10.4f %10.4f %10.4f\n",
              i, _bmat[i], Bmat[i], e = _bmat[i]-Bmat[i] );
    if ( (e = fabs(e)) > maxe ) {
      maxe = e;
      i0 = i;
    }
  }
  printf ( "%4d: %10.4f %10.4f %10.4f\n",
           i0, _bmat[i0], Bmat[i0], _bmat[i0]-Bmat[i0] );
  
  fclose ( f );
  printf ( "splmatrf.txt\n" );
} /*CompareMatrices*/

void TestMatrix ( void )
{
  g1h_Q2DrawSplMatricesf ( domain, SaveMatrices );
  Amat = malloc ( m_size*sizeof(float) );
  Bmat = malloc ( nf_a*nf_b*sizeof(float) );
  if ( !Amat || !Bmat )
    exit ( 1 );
  ComputeMatrices ();
  CompareMatrices ();
  ps_OpenFile ( fn1, 600 );
  PaintMatrixError ( m_k, m_r, m_s, _amat, Amat, _bmat, Bmat, 0.05 );
  ps_CloseFile ();
  free ( Amat );
  free ( Bmat );
  free ( _amat );
  free ( _bmat );
} /*TestMatrix*/

/* ///////////////////////////////////////////////////////////////////////// */
int GetOption ( GHoleDomainf *domain, int query, int qn,
                int *ndata, int **idata, float **fdata )
{
  switch ( query ) {
/*
case G1HQUERY_BASIS:
    return G1H_USE_RESTRICTED_BASIS;
*/
case G1HQUERY_Q2_FORM_CONSTANT:
    *ndata = 1;
    *fdata = &q2const;
    return G1H_Q2_USE_SUPPLIED_CONSTANT;
default:
    return G1H_DEFAULT;
  }
} /*GetOption*/

void OutputSplPatch ( int degu, int lknu, const float *knotsu,
                      int degv, int lknv, const float *knotsv,
                      const float *cp, void *usrptr )
{
  int s;

  if ( nfpatches == 0 ) {
    fpknots = malloc ( (lknu+1)*sizeof(float) );
    if ( !fpknots )
      exit ( 1 );
    memcpy ( fpknots, knotsu, (lknu+1)*sizeof(float) );
    fpdeg = degu;
    fplkn = lknu;
    printf ( "%d %d %d %d\n", degu, lknu, degv, lknv );
  }
  s = (lknu-degu)*(lknu-degu);
  fpcp[nfpatches] = malloc ( s*sizeof(point3f) );
  if ( !fpcp[nfpatches] )
    exit ( 1 );
  memcpy ( fpcp[nfpatches], cp, s*sizeof(point3f) );
  nfpatches ++;
} /*OutputSplPatch*/

int main ( void )
{
  struct tms start, stop;
  float time;
  float param[2] = {0.0,0.0};
  float *acoeff;
  
  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 67108864 ); /* 64 MB */
  acoeff = pkv_GetScratchMemf ( 90 );
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    g1h_SetOptionProcf ( domain, GetOption );
    if ( g1h_ComputeSplBasisf ( domain, NK, M1, M2 ) ) {
      g1h_DrawSplBasFuncNumf ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      printf ( "a=%d, b=%d, c=%d, d=%d\n", nfunc_a, nfunc_b, nfunc_c, nfunc_d );
      printf ( "    A    B    C    D\n" );
      printf ( "%5d%5d%5d%5d%5d\n", 0, nfunc_a, nfunc_a+nfunc_b,
               nfunc_a+nfunc_b+nfunc_c, nfunc_a+nfunc_b+nfunc_c+nfunc_d );
      printf ( "%5d     %5d%5d\n", nfunc_c+nfunc_d, 0, nfunc_c );
      if ( g1h_Q2SplComputeFormMatrixf ( domain ) ) {
        nfinalp = 0;
        InitHole ( HOLE_K, iparams );
        if ( !g1h_Q2SplFillHolef ( domain, 3, (float*)Bpt, NULL, NULL, OutputSplPatch ) )
          exit ( 1 );
        times ( &stop );
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
        MakePictures ();
#ifdef TESTMATRIX
        TestMatrix ();
#endif
      }
    }
    gh_DestroyDomainf ( domain );
  }
  pkv_DestroyScratchMem ();
  printf ( "\a\n" );
  exit ( 0 );
} /*main*/

