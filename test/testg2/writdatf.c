
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>

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

#include "writdatf.h"

/* ///////////////////////////////////////////////////////////////////////// */
static FILE *f;

void WriteKnots ( int hole_k, const float knots[][11] )
{
  int i, j;

  fprintf ( f, "# knots\n" );
  fprintf ( f, "k %d\n", hole_k );
  for ( i = 0;  i < hole_k;  i++ ) {
    for ( j = 0; j < 11; j++ )
      fprintf ( f, "%f ", knots[i][j] );
    fprintf ( f, "\n" );
  }
} /*WriteKnots*/

void WriteDomCP ( int hole_np, const point2f *domcp )
{
  int i;

  fprintf ( f, "# domain control points\n" );
  fprintf ( f, "d %d\n", hole_np );
  for ( i = 0; i < hole_np; i++ )
    fprintf ( f, "%f %f\n", domcp[i].x, domcp[i].y );
} /*WriteDomCP*/

void WriteSurfCP ( int hole_np, const point3f *surfcp )
{
  int i;

  fprintf ( f, "# surface control points\n" );
  fprintf ( f, "s %d\n", hole_np );
  for ( i = 0; i < hole_np; i++ )
    fprintf ( f, "%f %f %f\n", surfcp[i].x, surfcp[i].y, surfcp[i].z );
} /*WriteSurfCP*/

void WriteCamera ( const CameraRecf *CPos )
{
  fprintf ( f, "# camera parameters\n" );
  fprintf ( f, "c %f %f %f %f %f %f %f %d %d\n",
    CPos->position.x, CPos->position.y, CPos->position.z,
    CPos->psi, CPos->theta, CPos->phi,
    CPos->vd.persp.f, CPos->width, CPos->height );
} /*WriteCamera*/

void WriteLightDir ( const vector3f *ld )
{
  fprintf ( f, "# Light direction\n" );
  fprintf ( f, "l %f %f %f\n", ld->x, ld->y, ld->z );
} /*WriteLightDir*/

void WriteBezPatch ( int n, int m, const float *p )
{
  int i;
  point3f *cp;

  cp = (point3f*)p;
  fprintf ( f, "# Bezier patch\n" );
  fprintf ( f, "b %d %d\n", n, m );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    fprintf ( f, "%f %f %f\n", cp[i].x, cp[i].y, cp[i].z );
} /*WriteBezPatch*/

void WriteBSPatch ( int n, int lknu, const float *knotsu,
                    int m, int lknv, const float *knotsv, const point3f *cp )
{
  void  *sp;
  int   ku, kv, pitch1, pitch2, pitch3;
  float *b, *c;
  int   i, j, start;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knotsu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knotsv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  if ( b && c ) {
    fprintf ( f, "# B-spline patch\n" );
    mbs_BSPatchToBezf ( 3, n, lknu, knotsu, m, lknv, knotsv, pitch1, (float*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        WriteBezPatch ( n, m, c );
      }
  }

  pkv_SetScratchMemTop ( sp );
} /*WriteBSPatch*/

static float *GetKnotSequencef ( int hole_k, float knots[][11], int i )
{
  if ( i < 0 ) i += hole_k;
    else if ( i >= hole_k ) i -= hole_k;

  return knots[i];
} /*GetKnotSequencef*/

void WriteSurfPatches ( int hole_k, float knots[][11], const point3f *surfcp )
{
  void    *sp;
  int     *ind, i, j, k;
  point3f *cp, *cq;
  float   *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp = pkv_GetScratchMem ( 32*sizeof(point3f) );
  if ( ind && cp ) {
    cq = &cp[16];
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ ) {
        ukn = GetKnotSequencef ( hole_k, knots, i-1 );  ukn += 3;
        vkn = GetKnotSequencef ( hole_k, knots, i );    vkn += j;
        gh_GetBspInd ( hole_k, i, j, ind );
        for ( k = 0; k < 16; k++ )
          cp[k] = surfcp[ind[k]];
        mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)cp,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)cq );
        WriteBezPatch ( 3, 3, (float*)cq );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*WriteSurfPatches*/

void OpenDataFile ( const char *fn )
{
  f = fopen ( fn, "w+" );
} /*OpenDataFile*/

void CloseDataFile ( void )
{
  fprintf ( f, "e\n" );
  fclose ( f );
} /*CloseDataFile*/

