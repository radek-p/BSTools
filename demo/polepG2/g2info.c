
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "raybez.h"
#include "eg2holef.h"

#include "oldxgedit.h"
#include "datagenf.h"
#include "edg2hole.h"
#include "g2ekernel.h"

/* ///////////////////////////////////////////////////////////////////////// */
static void WritePartitionInfo ()
{
  void  *sp;
  int   hole_k, hole_m;
  float *partition;
  int   i;

  sp = pkv_GetScratchMemTop ();
  partition = pkv_GetScratchMemf ( domain->hole_k*sizeof(float) );
  if ( partition ) {
    g2h_ExtractPartitionf ( domain, &hole_k, &hole_m, partition,
                            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
    printf ( "k = %d, m = %d\n", hole_k, hole_m );
    printf ( "partition:\n" );
    for ( i = 0; i < hole_k; i++ )
      printf ( "%7.4f (%7.2f)\n", partition[i], 180.0/PI*partition[i] );
  }
  pkv_SetScratchMemTop ( sp );
} /*WritePartitionInfo*/

/* ///////////////////////////////////////////////////////////////////////// */
void InitAVector ( int n, float *x )
{
  int   i;
  float xl;

  x[0] = 1.0;
  for ( i = 1; i < n; i++ )
    x[i] = -0.97*x[i-1];
  xl = 1.0/pkn_SecondNormf ( n, x );
  for ( i = 0; i < n; i++ )
    x[i] *= xl;
} /*InitVector*/

void FindCondNumber ( int n, const float *a )
{
  void  *sp;
  float *b, *x, *y, *z;
  float lambda_1, lambda_n, pr, yl, ca;
  int   i;

  sp = pkv_GetScratchMemTop ();
  /* power method of finding the greatest eigenvalue */

  b = pkv_GetScratchMemf ( n*(n+1)/2 );
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  z = pkv_GetScratchMemf ( n );
  if ( !b || !x || !y || !z )
    exit ( 1 );
  memcpy ( b, a, (n*(n+1)/2)*sizeof(float) );

  InitAVector ( n, x );
  for ( i = 0; i < 100; i++ ) {
    pkn_SymMatrixMultf ( n, b, 1, 1, x, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
  lambda_1 = yl;

  pkn_CholeskyDecompf ( n, b );
  InitAVector ( n, x );
  for ( i = 0; i < 100; i++ ) {
    pkn_LowerTrMatrixSolvef ( n, b, 1, 1, x, 1, z );
    pkn_UpperTrMatrixSolvef ( n, b, 1, 1, z, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
  lambda_n = 1.0/yl;
  printf ( "l_1 = %g,  l_n = %g,  cond_2(A) = %g\n",
           lambda_1, lambda_n, lambda_1*yl );

  pkv_SetScratchMemTop ( sp );
} /*FindCondNumber*/

static void DrawMatrix ( int nfa, int nfb, float *amat, float *bmat )
{
  printf ( "Coons space: dim V_0 = %d\n", nfa );
  FindCondNumber ( nfa, amat );
} /*DrawMatrix*/

static void WriteMatrixInfo ()
{
  g2h_DrawMatricesf ( domain, DrawMatrix );
} /*WriteMatrixInfo*/

static void EmptyOutput ( int n, int m, const float *cp, void *usrptr )
{
} /*EmptyOutput*/

static void WriteFuncValue ()
{
  void  *sp;
  float fval, *acoeff;

  sp = pkv_GetScratchMemTop ();
  acoeff = pkv_GetScratchMemf ( 90 );
  if ( acoeff ) {
    g2h_FillHolef ( domain, 3, (float*)surfcp, acoeff, NULL, EmptyOutput );
    fval = g2h_FunctionalValuef ( domain, 3, (float*)surfcp, acoeff );
    printf ( "F(u) = %f\n", fval );
  }
  pkv_SetScratchMemTop ( sp );
} /*WriteFuncValue*/

static void WriteConstrInfo ()
{
  void     *sp;
  float    fval, *acoeff, *rhsf;
  vector3f *rhs;

  sp = pkv_GetScratchMemTop ();
  acoeff = pkv_GetScratchMemf ( 90 );
  if ( acoeff ) {
    if ( !swNormalConstr ) {
      if ( !(rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) )) )
        goto wayout;
      SetupConstraintsRHS ( rhs );
      g2h_FillHoleConstrf ( domain, 3, (float*)surfcp,
                            nconstr, (float*)rhs, acoeff, NULL, EmptyOutput );
    }
    else {
      if ( !(rhsf = pkv_GetScratchMemf ( nconstr )) )
        goto wayout;
      SetupNConstraintsRHS ( rhsf );
      g2h_FillHoleAltConstrf ( domain, 3, (float*)surfcp,
                               nconstr, rhsf, acoeff, NULL, EmptyOutput );
    }
    fval = g2h_FunctionalValuef ( domain, 3, (float*)surfcp, acoeff );
    printf ( "%d constraints, F(u) = %f\n", nconstr, fval );
  }
wayout:
  pkv_SetScratchMemTop ( sp );
} /*WriteConstrInfo*/

void FindExtCondNumber ( int k, int r, int s, float *Aii )
{
#define ITER 50
  void  *sp;
  int   n, i;
  float *Lii, *x, *y;
  float ca, pr, yl, lmin, lmax, size;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  Lii = pkv_GetScratchMemf ( (size = pkn_Block1ArraySize ( k, r, s )) );
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  if ( Lii && x && y ) {
    memcpy ( Lii, Aii, size*sizeof(float) );
    if ( !pkn_Block1CholeskyDecompMf ( k, r, s, Lii ) )
      exit ( 1 );

    InitAVector ( n, x );
    for ( i = 0; i < ITER; i++ ) {
      pkn_Block1SymMatrixMultf ( k, r, s, Aii, 1, 1, x, 1, y );
      pr = pkn_ScalarProductf ( n, x, y );
      yl = pkn_SecondNormf ( n, y );
      ca = pr/yl;
      pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmax = yl;

    InitAVector ( n, x );
    for ( i = 0; i < ITER; i++ ) {
      memcpy ( y, x, n*sizeof(float) );
      pkn_Block1LowerTrMSolvef ( k, r, s, Lii, 1, 1, y );
      pkn_Block1UpperTrMSolvef ( k, r, s, Lii, 1, 1, y );
      pr = pkn_ScalarProductf ( n, x, y );
      yl = pkn_SecondNormf ( n, y );
      ca = pr/yl;pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmin = 1.0/yl;
    printf ( "l_1 = %g,  l_n = %g,  cond_2(A) = %g\n", lmax, lmin, lmax*yl );
  }

  pkv_SetScratchMemTop ( sp );
#undef ITER
} /*FindExtCondNumber*/

static void DrawExtMatrix ( int k, int r, int s, float *Aii, float *Bi )
{
  printf ( "Bezier space: dim V_0 = %d\n", k*r+s );
  FindExtCondNumber ( k, r, s, Aii );
} /*DrawExtMatrix*/

static void WriteExtMatrixInfo ()
{
  g2h_DrawExtMatricesf ( domain, DrawExtMatrix );
} /*WriteExtMatrixInfo*/

static void WriteExtFuncValue ()
{
  void  *sp;
  float fval, *acoeff;

  sp = pkv_GetScratchMemTop ();
  acoeff = pkv_GetScratchMemf ( 3*158 );
  if ( acoeff ) {
    g2h_ExtFillHolef ( domain, 3, (float*)surfcp, acoeff, NULL, EmptyOutput );
    fval = g2h_ExtFunctionalValuef ( domain, 3, (float*)surfcp, acoeff );
    printf ( "F(u) = %f\n", fval );
  }
  pkv_SetScratchMemTop ( sp );
} /*WriteExtFuncValue*/

static void WriteExtConstrInfo ()
{
  void     *sp;
  float    fval, *acoeff, *rhsf;
  vector3f *rhs;

  sp = pkv_GetScratchMemTop ();
  acoeff = pkv_GetScratchMemf ( 3*158 );
  if ( acoeff ) {
    if ( !swNormalConstr ) {
      if ( !(rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) )) )
        goto wayout;
      SetupConstraintsRHS ( rhs );
      g2h_ExtFillHoleConstrf ( domain, 3, (float*)surfcp,
                               nconstr, (float*)rhs, acoeff,
                               NULL, EmptyOutput );
    }
    else {
      if ( !(rhsf = pkv_GetScratchMemf ( nconstr )) )
        goto wayout;
      SetupNConstraintsRHS ( rhsf );
      g2h_ExtFillHoleAltConstrf ( domain, 3, (float*)surfcp,
                                  nconstr, rhsf, acoeff, NULL,
                                  EmptyOutput );
    }
    fval = g2h_ExtFunctionalValuef ( domain, 3, (float*)surfcp, acoeff );
    printf ( "%d constraints, F(u) = %f\n", nconstr, fval );
  }
wayout:
  pkv_SetScratchMemTop ( sp );
} /*WriteExtConstrInfo*/

/* ///////////////////////////////////////////////////////////////////////// */
void WriteInfo ()
{
  WritePartitionInfo ();
  if ( swCoonsPatches ) {
    WriteMatrixInfo ();
    WriteFuncValue ();
    if ( swConstraintsOn && nconstr )
      WriteConstrInfo ();
  }
  else {
    WriteExtMatrixInfo ();
    WriteExtFuncValue ();
    if ( swConstraintsOn && nconstr )
      WriteExtConstrInfo ();
  }
} /*WriteInfo*/

