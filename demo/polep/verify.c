
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

static void Verify ( const float *b0, const float *c0,
                     const float *f0, const float *g0,
                     const point3f *q0, const vector3f *q0v,
                     const point3f *r0, const vector3f *r0v )
{
  point3f  p;
  vector3f qu, qv, quu, quv, ru, rv, ruu, ruv, nv, x, y, z;
  float    b, c, f, g, db, dc, df, dg;

  printf ( "\n" );
  p = *q0;
  printf ( "Punkt srodkowy: %f %f %f\n", p.x, p.y, p.z );
  printf ( "Blad pozycyjny: %g %g %g\n", p.x-r0->x, p.y-r0->y, p.z-r0->z );

  SubtractPoints3f ( &q0[1], &q0[0], &qu );
  MultVector3f ( 3.0, &qu, &qu );
  qv = q0v[0];
  SubtractPoints3f ( &q0v[1], &q0v[0], &quv );
  MultVector3f ( 3.0, &quv, &quv );
  AddVector3f ( &q0[0], &q0[2], &quu );
  AddVector3Mf ( &quu, &q0[1], -2.0, &quu );
  MultVector3f ( 6.0, &quu, &quu );
  printf ( "Qu : %f, %f, %f\n", qu.x, qu.y, qu.z );
  printf ( "Qv : %f, %f, %f\n", qv.x, qv.y, qv.z );
  printf ( "Quu: %f, %f, %f\n", quu.x, quu.y, quu.z );
  printf ( "Quv: %f, %f, %f\n", quv.x, quv.y, quv.z );

  SubtractPoints3f ( &r0[1], &r0[0], &ru );
  MultVector3f ( 3.0, &ru, &ru );
  rv = r0v[0];
  SubtractPoints3f ( &r0v[1], &r0v[0], &ruv );
  MultVector3f ( 3.0, &ruv, &ruv );
  AddVector3f ( &r0[0], &r0[2], &ruu );
  AddVector3Mf ( &ruu, &r0[1], -2.0, &ruu );
  MultVector3f ( 6.0, &ruu, &ruu );
  printf ( "Ru : %f, %f, %f\n", ru.x, ru.y, ru.z );
  printf ( "Rv : %f, %f, %f\n", rv.x, rv.y, rv.z );
  printf ( "Ruu: %f, %f, %f\n", ruu.x, ruu.y, ruu.z );
  printf ( "Ruv: %f, %f, %f\n", ruv.x, ruv.y, ruv.z );
  CrossProduct3f ( &qu, &ru, &nv );

  b = b0[0];  c = c0[0];  f = f0[0];  g = g0[0];
  db = 3.0*(b0[1]-b0[0]);  dc = c0[2]-c0[0];
  df = 3.0*(f0[1]-f0[0]);  dg = g0[2]-g0[0];
  printf ( " b %f, db %f, c %f, dc %f\n f %f, df %f, g %f, dg %f\n",
           b, db, c, dc, f, df, g, dg );

  MultVector3f ( b, &qu, &x );
  AddVector3Mf ( &x, &qv, c, &x );
  printf ( "Blad rv: %g, %g, %g\n", x.x-ru.x, x.y-ru.y, x.z-rv.z );
  MultVector3f ( f, &ru, &x );
  AddVector3Mf ( &x, &rv, g, &x );
  printf ( "Blad qv: %g, %g, %g\n", x.x-qu.x, x.y-qu.y, x.z-qv.z );

  MultVector3f ( -b, &quu, &x );
  AddVector3Mf ( &x, &quv, -c, &x );
  AddVector3Mf ( &x, &ruu, f, &x );
  AddVector3Mf ( &x, &ruv, g, &x );
  printf ( "Blad v00: %g\n", DotProduct3f(&x, &nv) );


  MultVector3f ( -b, &quu, &z );
  AddVector3Mf ( &z, &quv, -c, &z );
  AddVector3Mf ( &z, &ruu, f, &z );
  AddVector3Mf ( &z, &ruv, g, &z );
  AddVector3Mf ( &z, &qv, -dc, &z );
  AddVector3Mf ( &z, &rv, dg, &z );

printf ( "V00: %f, %f, %f\n", z.x, z.y, z.z );

/*
  x = qu;
  y = ru;
  Solve3x2 ( (float*)&x, (float*)&y, (float*)&z );
  db = z.x;
  df = -z.y;
  printf ( "db %f, df %f\n", db, df );
*/

  MultVector3f ( db, &qu, &x );
  AddVector3Mf ( &x, &qv, dc, &x );
  AddVector3Mf ( &x, &quu, b, &x );
  AddVector3Mf ( &x, &quv, c, &x );
  MultVector3f ( df, &ru, &y );
  AddVector3Mf ( &y, &rv, dg, &y );
  AddVector3Mf ( &y, &ruu, f, &y );
  AddVector3Mf ( &y, &ruv, g, &y );
  printf ( "x: %f, %f, %f\n", x.x, x.y, x.z );
  printf ( "y: %f, %f, %f\n", y.x, y.y, y.z );
  printf ( "Blad x-y: %g, %g, %g\n", x.x-y.x, x.y-y.y, x.z-y.z );

  printf ( "\n" );

exit ( 0 );

} /*Verify*/

