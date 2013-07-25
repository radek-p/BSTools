
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


void EqMerSetQuarterCircleArc ( void )
{
  xge_KnotWind *ckwind;
  double        kn[2] = {0.0, 1.0};

  ckwind = eq_ckwind.er->data0;
  ckwind->degree = 2;
  eqmer_nurbs = true;
  if ( equator )
    equator_nurbs = true;
  else
    meridian_nurbs = true;
  eqmer_closed = false;
  arc_angle = (double)(0.5*PI);
  mbs_SetKnotPatternd ( 1, kn, 3, &ckwind->lastknot, ckwind->knots );
  SetPoint3d ( &eqmer_cpoints[0], 1.0, 0.0, 1.0 );
  SetPoint3d ( &eqmer_cpoints[1], (double)(0.5*SQRT2), (double)(0.5*SQRT2),
               (double)(0.5*SQRT2) );
  SetPoint3d ( &eqmer_cpoints[2], 0.0, 1.0, 1.0 );
  neqmerpoints = 3;
  ProjectEqMerCurve ();
} /*EqMerSetQuarterCircleArc*/

void EqMerSetHalfCircleArc ( void )
{
  xge_KnotWind *ckwind;
  double        kn[3] = {0.0, 0.5, 1.0};

  ckwind = eq_ckwind.er->data0;
  ckwind->degree = 2;
  eqmer_nurbs = true;
  eqmer_closed = false;
  arc_angle = PI;
  ckwind->knots[0] = 0.0;
  mbs_SetKnotPatternd ( 2, kn, 2, &ckwind->lastknot, &ckwind->knots[1] );
  ckwind->knots[7] = 1.0;
  ckwind->lastknot = 7;
  if ( equator ) {
    equator_nurbs = true;
    SetPoint3d ( &eqmer_cpoints[0], 1.0, 0.0, 1.0 );
    SetPoint3d ( &eqmer_cpoints[1], (double)(0.5*SQRT2), (double)(0.5*SQRT2),
                 (double)(0.5*SQRT2) );
    SetPoint3d ( &eqmer_cpoints[2], 0.0, 1.0, 1.0 );
    SetPoint3d ( &eqmer_cpoints[3], (double)(-0.5*SQRT2), (double)(0.5*SQRT2),
                 (double)(0.5*SQRT2) );
    SetPoint3d ( &eqmer_cpoints[4], -1.0, 0.0, 1.0 );
  }
  else {
    meridian_nurbs = true;
    SetPoint3d ( &eqmer_cpoints[0], 0.0, -1.0, 1.0 );
    SetPoint3d ( &eqmer_cpoints[1], (double)(0.5*SQRT2), (double)(-0.5*SQRT2),
                 (double)(0.5*SQRT2) );
    SetPoint3d ( &eqmer_cpoints[2], 1.0, 0.0, 1.0 );
    SetPoint3d ( &eqmer_cpoints[3], (double)(0.5*SQRT2), (double)(0.5*SQRT2),
                 (double)(0.5*SQRT2) );
    SetPoint3d ( &eqmer_cpoints[4], 0.0, 1.0, 1.0 );
  }
  neqmerpoints = 5;
  ProjectEqMerCurve ();
} /*EqMerSetHalfCircleArc*/

void EqMerSetFullCircleArc ( void )
{
  xge_KnotWind *ckwind;
  double        kn[5] = {0.0, 0.25, 0.5, 0.75, 1.0};

  ckwind = eq_ckwind.er->data0;
  ckwind->degree = 2;
  eqmer_nurbs = true;
  if ( equator )
    equator_nurbs = true;
  else
    meridian_nurbs = true;
  eqmer_closed = false;
  arc_angle = (double)(2.0*PI);
  ckwind->knots[0] = 0.0;
  mbs_SetKnotPatternd ( 4, kn, 2, &ckwind->lastknot, &ckwind->knots[1] );
  ckwind->knots[11] = 1.0;
  ckwind->lastknot = 11;
  SetPoint3d ( &eqmer_cpoints[0], 1.0, 0.0, 1.0 );
  SetPoint3d ( &eqmer_cpoints[1], (double)(0.5*SQRT2), (double)(0.5*SQRT2),
               (double)(0.5*SQRT2) );
  SetPoint3d ( &eqmer_cpoints[2], 0.0, 1.0, 1.0 );
  SetPoint3d ( &eqmer_cpoints[3], (double)(-0.5*SQRT2), (double)(0.5*SQRT2),
               (double)(0.5*SQRT2) );
  SetPoint3d ( &eqmer_cpoints[4], -1.0, 0.0, 1.0 );
  SetPoint3d ( &eqmer_cpoints[5], (double)(-0.5*SQRT2), (double)(-0.5*SQRT2),
               (double)(0.5*SQRT2) );
  SetPoint3d ( &eqmer_cpoints[6], 0.0, -1.0, 1.0 );
  SetPoint3d ( &eqmer_cpoints[7], (double)(0.5*SQRT2), (double)(-0.5*SQRT2),
               (double)(0.5*SQRT2) );
  SetPoint3d ( &eqmer_cpoints[8], 1.0, 0.0, 1.0 );
  neqmerpoints = 9;
  ProjectEqMerCurve ();
} /*EqMerSetFullCircleArc*/

void EqMerSetCircleArc ( double angle )
{
  xge_KnotWind *ckwind;
  int          narcs, i;
  double        *kn, a, w;
  trans2d      tr;
  double   kn1[2] = {0.0, 1.0};
  double   kn2[3] = {0.0, 0.5, 1.0};
  double   kn3[4] = {0.0, 0.333333333333333, 0.6666666666666666, 1.0};
  double   kn4[5] = {0.0, 0.25, 0.5, 0.75, 1.0};

  ckwind = eq_ckwind.er->data0;

  if ( angle < 0.0 ) angle += (double)(2.0*PI);
  angle = (double)(max ( 0.0, angle ));
  angle = (double)(min ( 2.0*PI, angle ));
  if ( angle <= 0.5*PI ) { narcs = 1;  kn = kn1; }
  else if ( angle <= PI ) { narcs = 2;  kn = kn2; }
  else if ( angle <= 1.5*PI ) { narcs = 3;  kn = kn3; }
  else { narcs = 4;  kn = kn4; }
  ckwind->degree = 2;
  eqmer_nurbs = true;
  eqmer_closed = false;
  arc_angle = angle;
  ckwind->knots[0] = 0.0;
  mbs_SetKnotPatternd ( narcs, kn, 2, &ckwind->lastknot, &ckwind->knots[1] );
  ckwind->lastknot += 2;
  ckwind->knots[ckwind->lastknot] = 1.0;
  a = (double)(angle/(2.0*narcs));
  w = (double)(cos ( a ));
  SetPoint3d ( &eqmer_cpoints[0], 1.0, 0.0, 1.0 );
  IdentTrans2d ( &tr );
  RotTrans2d ( &tr, a );
  for ( i = 0; i < 2*narcs; i++ )
    Trans2Point3d ( &tr, &eqmer_cpoints[i], &eqmer_cpoints[i+1] );
  for ( i = 1; i < 2*narcs; i += 2 )
    eqmer_cpoints[i].z = w;
  neqmerpoints = 2*narcs+1;
  ProjectEqMerCurve ();
} /*EqMerSetCircleArc*/

