
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkgeom.h"
#include "multibs.h"
#include "convh.h"
#include "camerad.h"
#include "xgedit.h"

#include "spl3d.h"


static void PovWritePoint3d ( FILE *f, point3d *p )
{
  fprintf ( f, "<%f,%f,%f>", p->x, p->y, p->z );
} /*PovWritePoint3d*/

static void PovWriteCamera ( FILE *f )
{
  point3d p, q, r;
  double   a, w;

  fprintf ( f, "camera {\n" );
  a = (double)swind.CPos[3].height/((double)swind.CPos[3].width*swind.CPos[3].aspect);
  w = (double)sqrt(1.0/(1.0+a*a));
  a = (double)(2.0*atan2 ( w, 2.0*swind.CPos[3].vd.persp.f ));  /* radians */
  a *= (double)(180.0/PI); /* to degrees */

  fprintf ( f, "  angle %f\n", a );
  fprintf ( f, "  location " );
  PovWritePoint3d ( f, &swind.CPos[3].position );
  fprintf ( f, "\n" );
  SetPoint3d ( &p, (double)(swind.CPos[3].width/2),
               (double)(swind.CPos[3].height/2), 1.0 );
  CameraUnProjectPoint3d ( &swind.CPos[3], &p, &q );
  fprintf ( f, "  look_at " );
  PovWritePoint3d ( f, &swind.CPos[3].g_centre );
  fprintf ( f, "\n" );
  SetPoint3d ( &p, (double)(swind.CPos[3].width/2), 0.0, 1.0 );
  CameraUnProjectPoint3d ( &swind.CPos[3], &p, &r );
  SubtractPoints3d ( &r, &q, &r );
  NormalizeVector3d ( &r );
  fprintf ( f, "  sky " );
  PovWritePoint3d ( f, &r );
  fprintf ( f, "\n" );
  fprintf ( f, "  right <-1.3333,0,0>\n" );
  fprintf ( f, "}\n\n" );
} /*PovWriteCamera*/

static void PovWriteBicubicPatch ( FILE *f, point3d *p )
{
  int i, j;

  fprintf ( f, "  bicubic_patch {\n" );
  fprintf ( f, "    type 1 u_steps 4 v_steps 4\n" );
  for ( i = 0; i < 4; i++ ) {
    fprintf ( f, "    " );
    for ( j = 0; j < 4; j++ )
      PovWritePoint3d ( f, &p[4*i+j] );
    fprintf ( f, "\n" );
  }
  fprintf ( f, "  }\n" );
} /*PovWriteBicubicPatch*/

static boolean PovWriteBezierPatch ( FILE *f, int n, int m, point4d *cp )
{
  /* the trouble is that PovRay does not accept NURBS patches, */
  /* nor it accepts Bezier patches other than bicubic. */
  /* A conversion is theferore necessary, with a possible loss */
  /* of the original shape */

  void    *sp;
  point3d *p, duu, dvv;

  sp = pkv_GetScratchMemTop ();
  p = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !p )
    goto failure;

    /* compute points and derivatives at corners */
  mbs_BCHornerDer2P3Rd ( n, m, cp, 0.0, 0.0,
                         &p[0], &p[4], &p[1], &duu, &p[5], &dvv );
  mbs_BCHornerDer2P3Rd ( n, m, cp, 1.0, 0.0,
                         &p[12], &p[8], &p[13], &duu, &p[9], &dvv );
  mbs_BCHornerDer2P3Rd ( n, m, cp, 0.0, 1.0,
                         &p[3], &p[7], &p[2], &duu, &p[6], &dvv );
  mbs_BCHornerDer2P3Rd ( n, m, cp, 1.0, 1.0,
                         &p[15], &p[11], &p[14], &duu, &p[10], &dvv );
    /* compute the bicubic Bezier representation */
  MultVector3d ( 1.0/3.0, &p[1], &p[1] );
  MultVector3d ( -1.0/3.0, &p[2], &p[2] );
  MultVector3d ( 1.0/3.0, &p[4], &p[4] );
  MultVector3d ( 1.0/9.0, &p[5], &p[5] );
  MultVector3d ( -1.0/9.0, &p[6], &p[6] );
  MultVector3d ( 1.0/3.0, &p[7], &p[7] );
  MultVector3d ( -1.0/3.0, &p[8], &p[8] );
  MultVector3d ( -1.0/9.0, &p[9], &p[9] );
  MultVector3d ( 1.0/9.0, &p[10], &p[10] );
  MultVector3d ( -1.0/3.0, &p[11], &p[11] );
  MultVector3d ( 1.0/3.0, &p[13], &p[13] );
  MultVector3d ( -1.0/3.0, &p[14], &p[14] );

  AddVector3d ( &p[0], &p[1], &p[1] );
  AddVector3d ( &p[4], &p[5], &p[5] );
  AddVector3d ( &p[1], &p[5], &p[5] );
  AddVector3d ( &p[0], &p[4], &p[4] );

  AddVector3d ( &p[3], &p[2], &p[2] );
  AddVector3d ( &p[7], &p[6], &p[6] );
  AddVector3d ( &p[2], &p[6], &p[6] );
  AddVector3d ( &p[3], &p[7], &p[7] );

  AddVector3d ( &p[12], &p[13], &p[13] );
  AddVector3d ( &p[8], &p[9], &p[9] );
  AddVector3d ( &p[13], &p[9], &p[9] );
  AddVector3d ( &p[12], &p[8], &p[8] );

  AddVector3d ( &p[15], &p[14], &p[14] );
  AddVector3d ( &p[11], &p[10], &p[10] );
  AddVector3d ( &p[14], &p[10], &p[10] );
  AddVector3d ( &p[15], &p[11], &p[11] );

  PovWriteBicubicPatch ( f, p );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*PovWriteBezierPatch*/

boolean PovRayExport ( void )
{
  void  *sp;
  FILE  *f;
  int   ku, kv, pitch1, pitch2, pitch3;
  double *b, *c;
  int   i, j, start;

  sp = pkv_GetScratchMemTop ();
  f = fopen ( "pomnij.pov", "w+" );
  if ( !f )
    return false;

  ku = mbs_NumKnotIntervalsd ( degree_u, lastknot_u, knots_u );
  kv = mbs_NumKnotIntervalsd ( degree_v, lastknot_v, knots_v );
  pitch1 = (lastknot_v-degree_v)*4;
  pitch2 = (degree_v+1)*4*kv;
  pitch3 = (degree_v+1)*4;
  b = pkv_GetScratchMemd ( pitch2*ku*(degree_u+1)*4 );
  c = pkv_GetScratchMemd ( (degree_u+1)*pitch3 );
  if ( b && c ) {
    fprintf ( f, "#include \"colors.inc\"\n\n" );
    fprintf ( f, "background { White }\n" );
    fprintf ( f, "light_source { <20,20,20> colour White }\n\n" );
    PovWriteCamera ( f );

    mbs_BSPatchToBezd ( 4, degree_u, lastknot_u, knots_u,
                        degree_v, lastknot_v, knots_v,
                        pitch1, (double*)cpoints,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    fprintf ( f, "union {\n" );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(degree_u+1)*pitch2+pitch3*j;
        pkv_Selectd ( degree_u+1, pitch3, pitch2, pitch3, &b[start], c );
        if ( !PovWriteBezierPatch ( f, degree_u, degree_v,
                                       (point4d*)c ) )
          goto failure;
      }
    fprintf ( f, "  pigment { Yellow }\n" );
    fprintf ( f, "  finish { ambient 0.3 phong 0.75 }\n" );
    fprintf ( f, "}\n\n" );
  }
  else
    goto failure;

  fclose ( f );
  printf ( "%s\n", "pomnij.pov" );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  fclose ( f );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*PovRayExport*/

