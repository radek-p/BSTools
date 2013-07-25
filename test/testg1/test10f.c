
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

char fn1[] = "g1q2_extmatrixf.ps";
char fn2[] = "g1q2extsurff.ps";

GHoleDomainf *domain;

int nfinalp;
point3f FinalCPoints[HOLE_K][121];

float iparams[3] = {0.5,0.5,0.5};

float q2const = 1.0/*20.0*/;

/* ///////////////////////////////////////////////////////////////////////// */
void PaintMatrix ( int k, int r, int s, float *A, float *B )
{
  void *sp;
  int  n, m, i, j, p;
  byte *buf;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  m = 6*HOLE_K+1;
  buf = pkv_GetScratchMem ( max(n,m) );
  if ( buf ) {
    ps_GSave ();
    ps_Write_Command ( "100 100 translate 40 dup scale" );
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
  g1h_Q2DrawExtMatricesf ( domain, DrawMatrices );
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
  for ( i = 0; i < nfinalp; i++ )
    DrawBezPatch3f ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i], 8, 8 );
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
float func1p[HOLE_K][36];
float func2p[4*HOLE_K][36];

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
  pkv_Selectf ( 36, 1, 3, 1, &f[0].z, func1p[nf1] );
  nf1 ++;
} /*OutBasf1*/

void OutBasf2 ( int n, int m, const point3f *f )
{
  pkv_Selectf ( 36, 1, 3, 1, &f[0].z, func2p[nf2] );
  nf2 ++;
} /*OutBasf2*/

void ComputeMatrices ( void )
{
#define NQUAD 16
  void   *sp;
  int    f1, f2, i, j, k, l, m;
  float  u, v;
  vector2f dp, dpu, dpv, dpuu, dpuv, dpvv, dpuuu, dpuuv, dpuvv, dpvvv;
  vector2f lgr1, lgr2;
  float   f1f, f1u, f1v, f1uu, f1uv, f1vv, f1uuu, f1uuv, f1uvv, f1vvv;
  float   g1u, g1v, g1uu, g1uv, g1vv, g1uuu, g1uuv, g1uvv, g1vvv;
  float   f2f, f2u, f2v, f2uu, f2uv, f2vv, f2uuu, f2uuv, f2uvv, f2vvv;
  float   g2u, g2v, g2uu, g2uv, g2vv, g2uuu, g2uuv, g2uvv, g2vvv;
  float   lap1, lap2, jac;
  double s, t;

  sp = pkv_GetScratchMemTop ();
  ndip = ndsp = 0;
  g1h_DrawDiPatchesf ( domain, OutDiPatches );
  gh_DrawDomSurrndPatchesf ( domain, OutDomSPatches );

  for ( f1 = 0; f1 < nf_a; f1++ ) {
    printf ( "%4d\b\b\b\b", f1 );
    if ( f1 >= nf_c ) {
      nf1 = 0;
      g1h_DrawBasAFunctionf ( domain, f1-nf_c, OutBasf1 );
    }
    else {
      memset ( func1p, 0, 36*HOLE_K*sizeof(float) );
      i = f1 / ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
      j = f1 % ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
      k = j / (G1H_FINALDEG-3) + 2;
      j = j % (G1H_FINALDEG-3) + 2;
      func1p[i][k*(G1H_FINALDEG+1)+j] = 1.0;
      nf1 = HOLE_K;
    }
    for ( f2 = 0; f2 <= f1; f2++ ) {
      if ( f2 >= nf_c ) {
        nf2 = 0;
        g1h_DrawBasAFunctionf ( domain, f2-nf_c, OutBasf2 );
      }
      else {
        memset ( func2p, 0, 36*HOLE_K*sizeof(float) );
        i = f2 / ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
        j = f2 % ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
        k = j / (G1H_FINALDEG-3) + 2;
        j = j % (G1H_FINALDEG-3) + 2;
        func2p[i][k*(G1H_FINALDEG+1)+j] = 1.0;
        nf2 = HOLE_K;
      }
      l = pkn_Block3FindElemPos ( m_k, m_r, m_s, f1, f2 );
      if ( l < 0 )
        continue;
      s = t = 0.0;
        /* integrate the Laplacian gradient */
      for ( k = 0; k < HOLE_K; k++ ) {
        for ( i = 0; i < NQUAD; i++ ) {
          u = (float)(2*i+1)/(float)(2*NQUAD);
          for ( j = 0; j < NQUAD; j++ ) {
            v = (float)(2*j+1)/(float)(2*NQUAD);
            mbs_BCHornerDer3Pf ( 5, 5, 2, (float*)dipatch[k], u, v,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv,
                (float*)&dpuuu, (float*)&dpuuv, (float*)&dpuvv, (float*)&dpvvv );
            jac = fabs(dpu.x*dpv.y-dpv.x*dpu.y);
            mbs_BCHornerDer3Pf ( 5, 5, 1, func1p[k], u, v,
                &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv, &f1uuu, &f1uuv, &f1uvv, &f1vvv );
            mbs_BCHornerDer3Pf ( 5, 5, 1, func2p[k], u, v,
                &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv, &f2uuu, &f2uuv, &f2uvv, &f2vvv );
            pkn_Comp2iDerivatives3f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, dpuuu.x, dpuuu.y, dpuuv.x, dpuuv.y,
                dpuvv.x, dpuvv.y, dpvvv.x, dpvvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &f1uuu, &f1uuv, &f1uvv, &f1vvv,
                &g1u, &g1v, &g1uu, &g1uv, &g1vv, &g1uuu, &g1uuv, &g1uvv, &g1vvv );
            pkn_Comp2iDerivatives3f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, dpuuu.x, dpuuu.y, dpuuv.x, dpuuv.y,
                dpuvv.x, dpuvv.y, dpvvv.x, dpvvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &f2uuu, &f2uuv, &f2uvv, &f2vvv,
                &g2u, &g2v, &g2uu, &g2uv, &g2vv, &g2uuu, &g2uuv, &g2uvv, &g2vvv );
            SetVector2f ( &lgr1, g1uuu+g1uvv, g1uuv+g1vvv );
            SetVector2f ( &lgr2, g2uuu+g2uvv, g2uuv+g2vvv );
            s += (lgr1.x*lgr2.x+lgr1.y*lgr2.y)*jac;
          }
        }
      }
        /* integrate the Laplacian jump */
      for ( k = 0; k < HOLE_K; k++ ) {
        j = (k+HOLE_K-1) % HOLE_K;
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], u, 0.0,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpu.x*dpu.x + dpu.y*dpu.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], u, 0.0,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[k], u, 0.0,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = g1uu + g1vv;
          lap2 = g2uu + g2vv;
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[j], 0.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[j], 0.0, u,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[j], 0.0, u,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 -= g1uu + g1vv;
          lap2 -= g2uu + g2vv;
          t += lap1*lap2*jac;
        }
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], u, 1.0,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpu.x*dpu.x + dpu.y*dpu.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], u, 1.0,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[k], u, 1.0,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = -(g1uu + g1vv);
          lap2 = -(g2uu + g2vv);
          t += lap1*lap2*jac;
        }
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], 1.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpv.x*dpv.x + dpv.y*dpv.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], 1.0, u,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[k], 1.0, u,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = -(g1uu + g1vv);
          lap2 = -(g2uu + g2vv);
          t += lap1*lap2*jac;
        }
      }
      Amat[l] = s/(float)(NQUAD*NQUAD) + q2const*t/(float)NQUAD;
    }
  }
  printf ( "\n" );

  for ( f1 = l = 0;  f1 < nf_a;  f1++ ) {
    printf ( "%4d\b\b\b\b", f1 );
    if ( f1 >= nf_c ) {
      nf1 = 0;
      g1h_DrawBasAFunctionf ( domain, f1-nf_c, OutBasf1 );
    }
    else {
      memset ( func1p, 0, 36*HOLE_K*sizeof(float) );
      i = f1 / ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
      j = f1 % ((G1H_FINALDEG-3)*(G1H_FINALDEG-3));
      k = j / (G1H_FINALDEG-3) + 2;
      j = j % (G1H_FINALDEG-3) + 2;
      func1p[i][k*(G1H_FINALDEG+1)+j] = 1.0;
      nf1 = HOLE_K;
    }
    for ( f2 = 0;  f2 < nf_b;  f2++, l++ ) {
      nf2 = 0;
      g1h_DrawBasBFunctionf ( domain, f2, OutBasf2 );
      s = t = 0.0;
        /* Integrate the Laplacian gradient */
      for ( k = 0; k < HOLE_K; k++ ) {
        for ( i = 0; i < NQUAD; i++ ) {
          u = (float)(2*i+1)/(float)(2*NQUAD);
          for ( j = 0; j < NQUAD; j++ ) {
            v = (float)(2*j+1)/(float)(2*NQUAD);
            mbs_BCHornerDer3Pf ( 5, 5, 2, (float*)dipatch[k], u, v,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv,
                (float*)&dpuuu, (float*)&dpuuv, (float*)&dpuvv, (float*)&dpvvv );
            jac = fabs(dpu.x*dpv.y-dpv.x*dpu.y);
            mbs_BCHornerDer3Pf ( 5, 5, 1, func1p[k], u, v,
                &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv, &f1uuu, &f1uuv, &f1uvv, &f1vvv );
            mbs_BCHornerDer3Pf ( 5, 5, 1, func2p[3*HOLE_K+k], u, v,
                &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv, &f2uuu, &f2uuv, &f2uvv, &f2vvv );
            pkn_Comp2iDerivatives3f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, dpuuu.x, dpuuu.y, dpuuv.x, dpuuv.y,
                dpuvv.x, dpuvv.y, dpvvv.x, dpvvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &f1uuu, &f1uuv, &f1uvv, &f1vvv,
                &g1u, &g1v, &g1uu, &g1uv, &g1vv, &g1uuu, &g1uuv, &g1uvv, &g1vvv );
            pkn_Comp2iDerivatives3f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, dpuuu.x, dpuuu.y, dpuuv.x, dpuuv.y,
                dpuvv.x, dpuvv.y, dpvvv.x, dpvvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &f2uuu, &f2uuv, &f2uvv, &f2vvv,
                &g2u, &g2v, &g2uu, &g2uv, &g2vv, &g2uuu, &g2uuv, &g2uvv, &g2vvv );
            SetVector2f ( &lgr1, g1uuu+g1uvv, g1uuv+g1vvv );
            SetVector2f ( &lgr2, g2uuu+g2uvv, g2uuv+g2vvv );
            s += (lgr1.x*lgr2.x+lgr1.y*lgr2.y)*jac;
          }
        }
      }
        /* Integrate the Laplacian jump */
      for ( k = 0; k < HOLE_K; k++ ) {
        j = (k+HOLE_K-1) % HOLE_K;
        m = (k+1) % HOLE_K;
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], u, 0.0,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpu.x*dpu.x + dpu.y*dpu.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], u, 0.0,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[3*HOLE_K+k], u, 0.0,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = g1uu + g1vv;
          lap2 = g2uu + g2vv;
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[j], 0.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[j], 0.0, u,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[3*HOLE_K+j], 0.0, u,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 -= g1uu + g1vv;
          lap2 -= g2uu + g2vv;
          t += lap1*lap2*jac;
        }
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], u, 1.0,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpu.x*dpu.x + dpu.y*dpu.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], u, 1.0,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[3*HOLE_K+k], u, 1.0,
                &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = -(g1uu + g1vv);
          lap2 = -(g2uu + g2vv);
          mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)dspatch[3*m+1], 0.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          mbs_BCHornerDer2Pf ( 3, 3, 1, func2p[3*m+1], 0.0, u,
                &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap2 += g2uu + g2vv;
          t += lap1*lap2*jac;
        }
        for ( i = 0; i < NQUAD; i++ ) {
          u = (2*i+1)/(float)(2*NQUAD);
          mbs_BCHornerDer2Pf ( 5, 5, 2, (float*)dipatch[k], 1.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          jac = sqrt ( dpv.x*dpv.x + dpv.y*dpv.y );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func1p[k], 1.0, u,
                 &f1f, &f1u, &f1v, &f1uu, &f1uv, &f1vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f1u, &f1v, &f1uu, &f1uv, &f1vv, &g1u, &g1v, &g1uu, &g1uv, &g1vv );
          mbs_BCHornerDer2Pf ( 5, 5, 1, func2p[3*HOLE_K+k], 1.0, u,
                 &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap1 = -(g1uu + g1vv);
          lap2 = -(g2uu + g2vv);
          mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)dspatch[3*k+2], 0.0, u,
                (float*)&dp, (float*)&dpu, (float*)&dpv,
                (float*)&dpuu, (float*)&dpuv, (float*)&dpvv );
          mbs_BCHornerDer2Pf ( 3, 3, 1, func2p[3*k+2], 0.0, u,
                &f2f, &f2u, &f2v, &f2uu, &f2uv, &f2vv );
          pkn_Comp2iDerivatives2f ( dpu.x, dpu.y, dpv.x, dpv.y, dpuu.x, dpuu.y,
                dpuv.x, dpuv.y, dpvv.x, dpvv.y, 1,
                &f2u, &f2v, &f2uu, &f2uv, &f2vv, &g2u, &g2v, &g2uu, &g2uv, &g2vv );
          lap2 += g2uu + g2vv;
          t += lap1*lap2*jac;
        }
      }
      Bmat[l] = s/(float)(NQUAD*NQUAD) + q2const*t/(float)NQUAD;
    }
  }
  pkv_SetScratchMemTop ( sp );
#undef NQUAD
} /*ComputeMatrices*/

void CompareMatrices ( void )
{
  FILE *f;
  int i;

  f = fopen ( "extmatr.txt", "w+" );
  for ( i = 0; i < m_size; i++ )
    fprintf ( f, "%3d: %10.4f %10.4f %10.4f\n", i, _amat[i], Amat[i], _amat[i]-Amat[i] );
  fprintf ( f, "\n" );
  for ( i = 0; i < nf_a*nf_b; i++ )
    fprintf ( f, "%3d: %10.4f %10.4f %10.4f\n", i, _bmat[i], Bmat[i], _bmat[i]-Bmat[i] );
  fclose ( f );
} /*CompareMatrices*/

void TestMatrix ( void )
{
  g1h_Q2DrawExtMatricesf ( domain, SaveMatrices );
  Amat = malloc ( m_size*sizeof(float) );
  Bmat = malloc ( nf_a*nf_b*sizeof(float) );
  if ( !Amat || !Bmat )
    exit ( 1 );
  ComputeMatrices ();
  CompareMatrices ();
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

void OutPatch ( int n, int m, const float *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, (n+1)*(m+1)*sizeof(point3f) );
  nfinalp ++;
} /*OutPatch*/

int main ( void )
{
  float param[2] = {0.0,0.0};
  float *acoeff;
  
  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 67108864 ); /* 64 MB */
  acoeff = pkv_GetScratchMemf ( 90 );
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    g1h_SetOptionProcf ( domain, GetOption );
    if ( g1h_ComputeBasisf ( domain ) ) {
      if ( g1h_Q2ExtComputeFormMatrixf ( domain ) ) {
        nfinalp = 0;
        InitHole ( HOLE_K, iparams );
        if ( !g1h_Q2ExtFillHolef ( domain, 3, (float*)Bpt, NULL, NULL, OutPatch ) )
          exit ( 1 );
        MakePictures ();
#ifdef TESTMATRIX
        TestMatrix ();
#endif
      }
    }
    gh_DestroyDomainf ( domain );
  }
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

