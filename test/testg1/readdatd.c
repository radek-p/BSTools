
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"

void SkipComment ( FILE *f )
{
  char c;

  do {
    fscanf ( f, "%c", &c );
  } while ( c != '\n' );
} /*SkipComment*/

void SkipLines ( FILE *f )
{
  int n, i;

  fscanf ( f, "%d", &n );
  for ( i = 0; i <= n; i++ )
    SkipComment ( f );
} /*SkipLines*/

void SkipNFloats ( FILE *f, int n )
{
  int   i;
  double x;

  for ( i = 0; i < n; i++ )
    fscanf ( f, "%lf", &x );
} /*SkipNFloats*/

void ReadKnots ( FILE *f, int *k, double *knots )
{
  int i;

  fscanf ( f, "%d", k );
  for ( i = 0; i < 11*(*k); i++ )
    fscanf ( f, "%lf", &knots[i] );
} /*ReadKnots*/

void ReadCamera ( FILE *f, CameraRecd *CPos )
{
  if ( CPos ) {
    CameraInitFramed ( CPos, false, true, 100, 100, 0, 0, 1.0, 4 );
    CameraInitPosd ( CPos );
    fscanf ( f, " %lf %lf %lf %lf %lf %lf %lf %hd %hd",
      &CPos->position.x, &CPos->position.y, &CPos->position.z,
      &CPos->psi, &CPos->theta, &CPos->phi, &CPos->vd.persp.f,
      &CPos->width, &CPos->height );
    CameraSetMappingd ( CPos );
  }
  else
    SkipNFloats ( f, 9 );
} /*ReadCamera*/

void ReadBezPatch ( FILE *f, int *n, int *m, point3d *cp )
{
  int   i;

  fscanf ( f, "%d %d", n, m );
  if ( cp ) {
    for ( i = 0; i < ((*n)+1)*((*m)+1); i++ )
      fscanf ( f, "%lf %lf %lf", &cp[i].x, &cp[i].y, &cp[i].z );
  }
  else
    SkipNFloats ( f, 3*((*n)+1)*((*m)+1) );
} /*ReadBezPatch*/

void ReadLightDir ( FILE *f, vector3d *lightdir )
{
  if ( lightdir ) {
    fscanf ( f, "%lf %lf %lf", &lightdir->x, &lightdir->y, &lightdir->z );
    NormalizeVector3d ( lightdir );
  }
  else
    SkipNFloats ( f, 3 );
} /*ReadLightDir*/

void ReadHighlight ( FILE *f, point3d *rp0, vector3d *rv1, vector3d *rv2 )
{
  if ( rp0 )
    fscanf ( f, "%lf %lf %lf", &rp0->x, &rp0->y, &rp0->z );
  else SkipNFloats ( f, 3 );
  if ( rv1 )
    fscanf ( f, "%lf %lf %lf", &rv1->x, &rv1->y, &rv1->z );
  else SkipNFloats ( f, 3 );
  if ( rv2 )
    fscanf ( f, "%lf %lf %lf", &rv2->x, &rv2->y, &rv2->z );
  else SkipNFloats ( f, 3 );
} /*ReadHighlight*/

void ReadDomCP ( FILE *f, point2d *dcp )
{
  int i, n;

  fscanf ( f, "%d", &n );
  if ( dcp ) {
    for ( i = 0; i < n; i++ )
      fscanf ( f, "%lf %lf", &dcp[i].x, &dcp[i].y );
  }
  else
    SkipNFloats ( f, 2*n );
} /*ReadDomCP*/

void ReadCP ( FILE *f, point3d *cp )
{
  int i, n;

  fscanf ( f, "%d", &n );
  if ( cp ) {
    for ( i = 0; i < n; i++ )
      fscanf ( f, "%lf %lf %lf", &cp[i].x, &cp[i].y, &cp[i].z );
  }
  else
    SkipNFloats ( f, 3*n );
} /*ReadCP*/

void ReadDatFile ( const char *fn, int *hole_k, double *knots,
                   point2d *dcp, point3d *cp, CameraRecd *CPos )
{
  FILE *f;
  char c;
  int  n, m;

  printf ( "Reading file %s\n", fn );
  f = fopen ( fn, "r+" );
  do {
    fscanf ( f, "%c", &c );
    switch ( c ) {
  case '#': SkipComment ( f );                      break;
  case 'd': ReadDomCP ( f, dcp );                   break;
  case 's': ReadCP ( f, cp );                       break;
  case 'k': ReadKnots ( f, hole_k, knots );         break;
  case 'c': ReadCamera ( f, CPos );                 break;
  case 'b': ReadBezPatch ( f, &n, &m, NULL );       break;
  case 'l': ReadLightDir ( f, NULL );               break;
  case 'h': ReadHighlight ( f, NULL, NULL, NULL );  break;
   default: break;
    }
  } while ( c != 'e' );
  fclose ( f );
} /*ReadDatFile*/

