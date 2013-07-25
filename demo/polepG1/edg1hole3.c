
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>
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
#include "eg1holef.h"

#include "oldxgedit.h"
#include "g1ekernel.h"
#include "edg1hole.h"

/* ///////////////////////////////////////////////////////////////////////// */
static FILE *f;

void WriteKnots ()
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

void WriteDomCP ()
{
  int i;

  fprintf ( f, "# domain control points\n" );
  fprintf ( f, "d %d\n", hole_np );
  for ( i = 0; i < hole_np; i++ )
    fprintf ( f, "%f %f\n", domcp[i].x, domcp[i].y );
} /*WriteDomCP*/

void WriteSurfCP ()
{
  int i;

  fprintf ( f, "# surface control points\n" );
  fprintf ( f, "s %d\n", hole_np );
  for ( i = 0; i < hole_np; i++ )
    fprintf ( f, "%f %f %f\n", surfcp[i].x, surfcp[i].y, surfcp[i].z );
} /*WriteSurfCP*/

void WriteCamera ()
{
  fprintf ( f, "# camera parameters\n" );
  fprintf ( f, "c %f %f %f %f %f %f %f %d %d\n",
    CPos.position.x, CPos.position.y, CPos.position.z,
    CPos.psi, CPos.theta, CPos.phi,
    CPos.vd.persp.f, 2*CPos.width, 2*CPos.height );
} /*WriteCamera*/

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

static float *GetKnotSequencef ( int i )
{
  if ( i < 0 ) i += hole_k;
    else if ( i >= hole_k ) i -= hole_k;

  return knots[i];
} /*GetKnotSequencef*/

void WriteSurfPatches ()
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
        ukn = GetKnotSequencef ( i-1 );  ukn += 3;
        vkn = GetKnotSequencef ( i );    vkn += j;
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

void WriteFinalPatches ()
{
  int i;
  fprintf ( f, "# hole filling patches\n" );
  for ( i = 0; i < hole_k; i++ )
    WriteBezPatch ( finaldegu, finaldegv, (float*)FinalCP[i] );
} /*WriteFinalPatches*/

void WriteNLFinalPatches ()
{
  int i;
  fprintf ( f, "# hole filling patches\n" );
  for ( i = hole_k; i < 2*hole_k; i++ )
    WriteBezPatch ( finaldegu, finaldegv, (float*)FinalCP[i] );
} /*WriteNLFinalPatches*/

void WriteFile ()
{
  char *fn;

  f = fopen ( fn = "dane.dat", "w+" );
  WriteCamera ();
  WriteKnots ();
  WriteDomCP ();
  WriteSurfCP ();
  WriteSurfPatches ();
  if ( FinalSurfValid )
    WriteFinalPatches ();
  if ( NLFinalSurfValid )
    WriteNLFinalPatches ();

  fprintf ( f, "e\n" );
  fclose ( f );
  printf ( "%s\n", fn );
} /*WriteFile*/

