
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindJFunctionsDif (
     const vector2f *omcbci, const vector2f *omcbcdi, const vector2f *omcbcddi,
     const vector2f *omcbcj, const vector2f *omcbcdj, const vector2f *omcbcddj,
     const vector2f *surrpcbci,
     const vector2f *surrpcbcj,
     float *b01, float *c01, float *f01, float *g01,
     float *b11, float *c11, float *f11, float *g11,
     float *b02, float *c02, float *f02, float *g02,
     float *b12, float *c12, float *f12, float *g12 )
{
/* the following #definitions are introduced for convenient accessing */
/* the partial derivatives of auxiliary and surrounding patches */
#define r0u0    omcbci[1]
#define r0uu0   omcbci[2]
#define r0uuu0  omcbci[3]
#define r0uuuu0 omcbci[4]
#define r0s0    omcbcdi[0]
#define r0us0   omcbcdi[1]
#define r0uus0  omcbcdi[2]
#define r0uuus0 omcbcdi[3]
#define r0ss0   omcbcddi[0]
#define r0uss0  omcbcddi[1]
#define r0uuss0 omcbcddi[2]
#define r0u1    omcbci[6]
#define r0uu1   omcbci[7]
#define r0s1    omcbcdi[4]
#define r0us1   omcbcdi[5]
#define r0uus1  omcbcdi[6]
#define r0ss1   omcbcddi[3]
#define r0uss1  omcbcddi[4]
#define r0uuss1 omcbcddi[5]
#define q0v0    omcbcj[1]
#define q0vv0   omcbcj[2]
#define q0vvv0  omcbcj[3]
#define q0vvvv0 omcbcj[4]
#define q0t0    omcbcdj[0]
#define q0vt0   omcbcdj[1]
#define q0vvt0  omcbcdj[2]
#define q0vvvt0 omcbcdj[3]
#define q0tt0   omcbcddj[0]
#define q0vtt0  omcbcddj[1]
#define q0vvtt0 omcbcddj[2]
#define q0v1    omcbcj[6]
#define q0t1    omcbcdj[4]
#define q0vv1   omcbcj[7]
#define q0vt1   omcbcdj[5]
#define q0tt1   omcbcddj[3]
#define q0vvt1  omcbcdj[6]
#define q0vtt1  omcbcddj[4]
#define q0vvtt1 omcbcddj[5]
#define q1v0    surrpcbci[3]
#define q1t0    surrpcbci[1]
#define q1vv0   surrpcbci[6]
#define q1vt0   surrpcbci[4]
#define q1tt0   surrpcbci[2]
#define q1vvt0  surrpcbci[7]
#define q1vtt0  surrpcbci[5]
#define q1vvtt0 surrpcbci[8]
#define q1v1    surrpcbci[9+3]
#define q1t1    surrpcbci[9+1]
#define q1vv1   surrpcbci[9+6]
#define q1vt1   surrpcbci[9+4]
#define q1tt1   surrpcbci[9+2]
#define q1vvt1  surrpcbci[9+7]
#define q1vtt1  surrpcbci[9+5]
#define q1vvtt1 surrpcbci[9+8]
#define r1u0    surrpcbcj[3]
#define r1s0    surrpcbcj[1]
#define r1uu0   surrpcbcj[6]
#define r1us0   surrpcbcj[4]
#define r1ss0   surrpcbcj[2]
#define r1uus0  surrpcbcj[7]
#define r1uss0  surrpcbcj[5]
#define r1uuss0 surrpcbcj[8]
#define r1u1    surrpcbcj[9+3]
#define r1s1    surrpcbcj[9+1]
#define r1uu1   surrpcbcj[9+6]
#define r1us1   surrpcbcj[9+4]
#define r1ss1   surrpcbcj[9+2]
#define r1uus1  surrpcbcj[9+7]
#define r1uss1  surrpcbcj[9+5]
#define r1uuss1 surrpcbcj[9+8]

  void     *sp;
  float    b010, b01d0, b01dd0, c010, c01d0, c01dd0,
           b011, b01d1, b01dd1, c011, c01d1, c01dd1,
           f010, f01d0, f01dd0, g010, g01d0, g01dd0,
           f011, f01d1, f01dd1, g011, g01d1, g01dd1,
           b110, b11d0, b11dd0, c110, c11d0, c11dd0,
           b111, b11d1, b11dd1, c111, c11d1, c11dd1,
           f110, f11d0, f11dd0, g110, g11d0, g11dd0,
           f111, f11d1, f11dd1, g111, g11d1, g11dd1,
           b020, b02d0, b02dd0, c020, c02d0, c02dd0,
           b021, b02d1, b02dd1, c021, c02d1, c02dd1,
           f020, f02d0, f02dd0, g020, g02d0, g02dd0,
           f021, f02d1, f02dd1, g021, g02d1, g02dd1,
           b120, b12d0, b12dd0, c120, c12d0, c12dd0,
           b121, b12d1, b12dd1, c121, c12d1, c12dd1,
           f120, f12d0, f12dd0, g120, g12d0, g12dd0,
           f121, f12d1, f12dd1, g121, g12d1, g12dd1;

  sp = pkv_GetScratchMemTop ();

        /* solve the compatibility equations */
            /* the construction of curves implies zero values of some */
            /* functions, so they are assigned to cancel rounding errors; */
            /* perhaps this might be optimized */
          /* function values at (0,0) */
  SolveCompatibilityEq1f ( &r0u0, &r0s0, &q0v0, &b010, &c010 );
  SolveCompatibilityEq1f ( &q0v0, &q0t0, &r0u0, &f010, &g010 );
          /* function values at (0,1) */
  SolveCompatibilityEq1f ( &r1u0, &r1s0, &q0v1, &b110, &c110 );
  SolveCompatibilityEq1f ( &q0v1, &q0t1, &r1u0, &f011, &g011 );
  b110 = f011 = 0.0;
          /* function values at (1,0) */
  SolveCompatibilityEq1f ( &r0u1, &r0s1, &q1v0, &b011, &c011 );
  SolveCompatibilityEq1f ( &q1v0, &q1t0, &r0u1, &f110, &g110 );
  b011 = f110 = 0.0;
          /* function values at (1,1) */
  SolveCompatibilityEq1f ( &r1u1, &r1s1, &q1v1, &b111, &c111 );
  SolveCompatibilityEq1f ( &q1v1, &q1t1, &r1u1, &f111, &g111 );
  b111 = f111 = 0.0;

            /* second order function values at the corners */
            /* are the unique solutions of the second order derivatives */
            /* compatibility conditions */
          /* function values at (0,0) */
  SolveCompatibilityEq3f ( &r0u0, &r0s0, &r0uu0, &r0us0, &r0ss0, &q0vv0,
                           b010, c010, &b020, &c020 );
  SolveCompatibilityEq3f ( &q0v0, &q0t0, &q0vv0, &q0vt0, &q0tt0, &r0uu0,
                           f010, g010, &f020, &g020 );
          /* function values at (0,1) */
  SolveCompatibilityEq3f ( &r1u0, &r1s0, &r1uu0, &r1us0, &r1ss0, &q0vv1,
                           b110, c110, &b120, &c120 );
  SolveCompatibilityEq3f ( &q0v1, &q0t1, &q0vv1, &q0vt1, &q0tt1, &r1uu0,
                           f011, g011, &f021, &g021 );
  b120 = f021 = 0.0;
          /* function values at (1,0) */
  SolveCompatibilityEq3f ( &r0u1, &r0s1, &r0uu1, &r0us1, &r0ss1, &q1vv0,
                           b011, c011, &b021, &c021 );
  SolveCompatibilityEq3f ( &q1v0, &q1t0, &q1vv0, &q1vt0, &q1tt0, &r0uu1,
                           f110, g110, &f120, &g120 );
  b021 = f120 = 0.0;
          /* function values at (1,1) */
  SolveCompatibilityEq3f ( &r1u1, &r1s1, &r1uu1, &r1us1, &r1ss1, &q1vv1,
                           b111, c111, &b121, &c121 );
  SolveCompatibilityEq3f ( &q1v1, &q1t1, &q1vv1, &q1vt1, &q1tt1, &r1uu1,
                           f111, g111, &f121, &g121 );
  b121 = f121 = 0.0;

          /* by assumption, c01 and g01 are polynomials of degree 1, */
          /* hence their derivatives are the following: */
  c01d0 = c01d1 = c011-c010;  c01dd0 = c01dd1 = 0.0;
  g01d0 = g01d1 = g011-g010;  g01dd0 = g01dd1 = 0.0;

          /* the values of derivatives of b01 and f01 at 0 are now */
          /* determined by the mixed partial derivatives compatibility */
          /* condition */
  SolveCompatibilityEq2af (
        &r0u0, &r0s0, &r0uu0, &r0us0, &q0v0, &q0t0, &q0vv0, &q0vt0,
        b010, c010, f010, g010, &b01d0, c01d0, &f01d0, g01d0 );

          /* by assumption, b01 and f01 are polynomials of degree 2, */
          /* hence their derivatives are the following: */
  b01d1 = (float)(2.0*(b011-b010)-b01d0);  b01dd0 = b01dd1 = b01d1-b01d0;
  f01d1 = (float)(2.0*(f011-f010)-f01d0);  f01dd0 = f01dd1 = f01d1-f01d0;

          /* now derivatives of b11, c11 and f11, g11 at 0 are determined */
          /* by mixed partial derivatives compatibility conditions */
          /* at (0,1) and (1,0) respectively */
  SolveCompatibilityEq2bf (
        &q0v1, &q0t1, &q0vv1, &q0vt1, &r1u0, &r1s0, &r1uu0, &r1us0,
        f011, g011, b110, c110, f01d1, g01d1, &b11d0, &c11d0 );
  SolveCompatibilityEq2bf (
        &r0u1, &r0s1, &r0uu1, &r0us1, &q1v0, &q1t0, &q1vv0, &q1vt0,
        b011, c011, f110, g110, b01d1, c01d1, &f11d0, &g11d0 );

          /* derivatives at the corners are now unique solutions of */
          /* the third order mixed partial derivatives compatibility */
          /* equations */
          /* at (0,0) */
  SolveCompatibilityEq4af (
            &q0v0, &q0t0, &q0vv0, &q0vt0, &q0vvv0, &q0vvt0,
            &r0u0, &r0s0, &r0uu0, &r0us0, &r0ss0, &r0uuu0, &r0uus0, &r0uss0,
            b010, b01d0, c010, c01d0, f010, f01d0, f01dd0, g010, g01d0, g01dd0,
            b020, &b02d0, c020, &c02d0 );
  SolveCompatibilityEq4af (
            &r0u0, &r0s0, &r0uu0, &r0us0, &r0uuu0, &r0uus0,
            &q0v0, &q0t0, &q0vv0, &q0vt0, &q0tt0, &q0vvv0, &q0vvt0, &q0vtt0,
            f010, f01d0, g010, g01d0, b010, b01d0, b01dd0, c010, c01d0, c01dd0,
            f020, &f02d0, g020, &g02d0 );
          /* at (0,1) */
  SolveCompatibilityEq4bf (
            &q0v1, &q0t1, &q0vv1, &q0vt1, &q0vvt1,
            &r1u0, &r1s0, &r1us0, &r1ss0, &r1uss0,
            b11d0, c110, c11d0, f01d1, f01dd1, g011, g01d1, g01dd1,
            &b12d0, c120, &c12d0 );
          /* at (1,0) */
  SolveCompatibilityEq4bf (
            &r0u1, &r0s1, &r0uu1, &r0us1, &r0uus1,
            &q1v0, &q1t0, &q1vt0, &q1tt0, &q1vtt0,
            f11d0, g110, g11d0, b01d1, b01dd1, c011, c01d1, c01dd1,
            &f12d0, g120, &g12d0 );
  
          /* by assumption, c02 and g02 are polynomials of degree 2, */
          /* hence their derivatives are the following: */
  c02d1 = (float)(2.0*(c021-c020)-c02d0);  c02dd0 = c02dd1 = c02d1-c02d0;
  g02d1 = (float)(2.0*(g021-g020)-g02d0);  g02dd0 = g02dd1 = g02d1-g02d0;

          /* the second order derivatives of b02 and f02 at 0 are */
          /* determined by the fourth order derivatives compatibility */
          /* condition */
  SolveCompatibilityEq5af (
            &q0v0, &q0t0, &q0vv0, &q0vt0, &q0tt0,
            &q0vvv0, &q0vvt0, &q0vtt0, &q0vvvv0, &q0vvvt0, &q0vvtt0,
            &r0u0, &r0s0, &r0uu0, &r0us0, &r0ss0,
            &r0uuu0, &r0uus0, &r0uss0, &r0uuuu0, &r0uuus0, &r0uuss0,
            b010, b01d0, b01dd0, c010, c01d0, c01dd0,
            f010, f01d0, f01dd0, g010, g01d0, g01dd0,
            b020, b02d0, &b02dd0, c020, c02d0, c02dd0,
            f020, f02d0, &f02dd0, g020, g02d0, g02dd0 );

          /* by assumption, b02 and f02 are polynomials of degree 3, */
          /* hence their derivatives are the following: */
  b02d1  = (float)(3.0*(b021-b020)-2.0*b02d0-0.5*b02dd0);
  b02dd1 = (float)(6.0*(b021-b020-b02d0)-2.0*b02dd0);
  f02d1  = (float)(3.0*(f021-f020)-2.0*f02d0-0.5*f02dd0);
  f02dd1 = (float)(6.0*(f021-f020-f02d0)-2.0*f02dd0);

          /* now the second order derivatives of b11, c11, b12, c12, */
          /* f11, g11, f12, g12 at 0 are determined by the third and */
          /* fourth order derivatives compatibility conditions */
          /* at (0,1) */
  SolveCompatibilityEq4cf (
            &r1u0, &r1s0, &r1uu0, &r1us0, &r1uus0,
            &q0v1, &q0t1, &q0vt1, &q0tt1, &q0vtt1,
            f01d1, g011, g01d1, b11d0, &b11dd0, c110, c11d0, &c11dd0,
            f02d1, g021, g02d1 );
  SolveCompatibilityEq5bf (
    &q0v1, &q0t1, &q0vv1, &q0vt1, &q0tt1, &q0vvt1, &q0vtt1, &q0vvtt1,
    &r1u0, &r1s0, &r1uu0, &r1us0, &r1ss0, &r1uus0, &r1uss0, &r1uuss0,
    b11d0, b11dd0, c110, c11d0, c11dd0, f01d1, f01dd1, g011, g01d1, g01dd1,
    b12d0, &b12dd0, c120, c12d0, &c12dd0, f02d1, f02dd1, g021, g02d1, g02dd1 );
          /* at (1,0) */
  SolveCompatibilityEq4cf (
            &q1v0, &q1t0, &q1vv0, &q1vt0, &q1vvt0,
            &r0u1, &r0s1, &r0us1, &r0ss1, &r0uss1,
            b01d1, c011, c01d1, f11d0, &f11dd0, g110, g11d0, &g11dd0,
            b02d1, c021, c02d1 );
  SolveCompatibilityEq5bf (
    &r0u1, &r0s1, &r0uu1, &r0us1, &r0ss1, &r0uus1, &r0uss1, &r0uuss1,
    &q1v0, &q1t0, &q1vv0, &q1vt0, &q1tt0, &q1vvt0, &q1vtt0, &q1vvtt0,
    f11d0, f11dd0, g110, g11d0, g11dd0, b01d1, b01dd1, c011, c01d1, c01dd1,
    f12d0, &f12dd0, g120, g12d0, &g12dd0, b02d1, b02dd1, c021, c02d1, c02dd1 );

          /* by assumption, c11 and g11 are polynomials of degree 3, */
          /* hence their derivatives are the following: */
  c11d1  = (float)(3.0*(c111-c110)-2.0*c11d0-0.5*c11dd0);
  c11dd1 = (float)(6.0*(c111-c110-c11d0)-2.0*c11dd0);
  g11d1  = (float)(3.0*(g111-g110)-2.0*g11d0-0.5*g11dd0);
  g11dd1 = (float)(6.0*(g111-g110-g11d0)-2.0*g11dd0);

          /* the derivatives of b11 and f11 at 1 must satisfy the */
          /* mixed partial derivatives compatibility condition at (1,1) */
  SolveCompatibilityEq2af (
        &r1u1, &r1s1, &r1uu1, &r1us1, &q1v1, &q1t1, &q1vv1, &q1vt1,
        b111, c111, f111, g111, &b11d1, c11d1, &f11d1, g11d1 );

          /* by assumption, b11 and f11 are polynomials of degree 4, */
          /* hence their second order derivatives at 1 are the following: */
  b11dd1 = (float)(12.0*(b110-b111)+6.0*(b11d0+b11d1)+b11dd0);
  f11dd1 = (float)(12.0*(f110-f111)+6.0*(f11d0+f11d1)+f11dd0);

          /* the derivatives of b12, c12, f12, g12 at 1 must satisfy the */
          /* third order derivatives compatibility condition at (1,1) */
  SolveCompatibilityEq4bf (
            &r1u1, &r1s1, &r1uu1, &r1us1, &r1uus1,
            &q1v1, &q1t1, &q1vt1, &q1tt1, &q1vtt1,
            f11d1, g111, g11d1, b11d1, b11dd1, c111, c11d1, c11dd1,
            &f12d1, g121, &g12d1 );
  SolveCompatibilityEq4bf (
            &q1v1, &q1t1, &q1vv1, &q1vt1, &q1vvt1,
            &r1u1, &r1s1, &r1us1, &r1ss1, &r1uss1,
            b11d1, c111, c11d1, f11d1, f11dd1, g111, g11d1, g11dd1,
            &b12d1, c121, &c12d1 );

          /* by assumption, c12 and g12 are polynomials of degree 4, */
          /* hence their second order derivatives at 1 are the following: */
  c12dd1 = (float)(12.0*(c120-c121)+6.0*(c12d0+c12d1)+c12dd0);
  g12dd1 = (float)(12.0*(g120-g121)+6.0*(g12d0+g12d1)+g12dd0);

          /* the second order derivatives of b12 and f12 now are determined */
          /* by the compatibility equations of mixed partial derivatives */
          /* of the fourth order */
  SolveCompatibilityEq5cf (
    &q1v1, &q1t1, &q1vv1, &q1vt1, &q1tt1, &q1vvt1, &q1vtt1, &q1vvtt1,
    &r1u1, &r1s1, &r1uu1, &r1us1, &r1ss1, &r1uus1,  &r1uss1, &r1uuss1,
    b11d1, b11dd1, c111, c11d1, c11dd1, f11d1, f11dd1, g111, g11d1, g11dd1,
    b12d1, &b12dd1, c121, c12d1, c12dd1, f12d1, &f12dd1, g121, g12d1, g12dd1 );


        /* now compute the coefficients of the polynomials in Bernstein bases */
  b01[0] = b010;  b01[1] = (float)(b010+0.5*b01d0);  b01[2] = b011;
  c01[0] = c010;  c01[1] = c011;
  f01[0] = f010;  f01[1] = (float)(f010+0.5*f01d0);  f01[2] = f011;
  g01[0] = g010;  g01[1] = g011;
  b11[0] = b110;  b11[1] = (float)(b110+0.25*b11d0);
  b11[2] = (float)(b110+0.5*b11d0+b11dd0/12.0);
  b11[3] = (float)(b111-0.25*b11d1);  b11[4] = b111;
  c11[0] = c110;  c11[1] = (float)(c110+c11d0/3.0);
  c11[2] = (float)(c111-c11d1/3.0);   c11[3] = c111;
  f11[0] = f110;  f11[1] = (float)(f110+0.25*f11d0);
  f11[2] = (float)(f110+0.5*f11d0+f11dd0/12.0);
  f11[3] = (float)(f111-0.25*f11d1);  f11[4] = f111;
  g11[0] = g110;  g11[1] = (float)(g110+g11d0/3.0);
  g11[2] = (float)(g111-g11d1/3.0);   g11[3] = g111;

  b02[0] = b020;  b02[1] = (float)(b020+b02d0/3.0);
  b02[2] = (float)(b021-b02d1/3.0);  b02[3] = b021;
  c02[0] = c020;  c02[1] = (float)(c020+0.5*c02d0);  c02[2] = c021;
  f02[0] = f020;  f02[1] = (float)(f020+f02d0/3.0);
  f02[2] = (float)(f021-f02d1/3.0);  f02[3] = f021;
  g02[0] = g020;  g02[1] = (float)(g020+0.5*g02d0);  g02[2] = g021;
  b12[0] = b120;  b12[1] = (float)(b120+0.2*b12d0);
  b12[2] = (float)(b120+0.4*b12d0+0.05*b12dd0);
  b12[3] = (float)(b121-0.4*b12d1+0.05*b12dd1);
  b12[4] = (float)(b121-0.2*b12d1);  b12[5] = b121;
  c12[0] = c120;  c12[1] = (float)(c120+0.25*c12d0);
  c12[2] = (float)(c120+0.5*c12d0+c12dd0/12.0);
  c12[3] = (float)(c121-0.25*c12d1);  c12[4] = c121;
  f12[0] = f120;  f12[1] = (float)(f120+0.2*f12d0);
  f12[2] = (float)(f120+0.4*f12d0+0.05*f12dd0);
  f12[3] = (float)(f121-0.4*f12d1+0.05*f12dd1);
  f12[4] = (float)(f121-0.2*f12d1);  f12[5] = f121;
  g12[0] = g120;  g12[1] = (float)(g120+0.25*g12d0);
  g12[2] = (float)(g120+0.5*g12d0+g12dd0/12.0);
  g12[3] = (float)(g121-0.25*g12d1);  g12[4] = g121;
/*
printf ( "b12: %f %f %f %f %f %f\n", b120, b12d0, b12dd0, b121, b12d1, b12dd1 );
printf ( "f12: %f %f %f %f %f %f\n", f120, f12d0, f12dd0, f121, f12d1, f12dd1 );
printf ( "c12: %f %f %f %f %f %f\n", c120, c12d0, c12dd0, c121, c12d1, c12dd1 );
printf ( "g12: %f %f %f %f %f %f\n\n", g120, g12d0, g12dd0, g121, g12d1, g12dd1 );

printf ( "b12: %f %f %f %f %f %f\n", b12[0], b12[1], b12[2], b12[3], b12[4], b12[5] );
printf ( "f12: %f %f %f %f %f %f\n", f12[0], f12[1], f12[2], f12[3], f12[4], f12[5] );
printf ( "c12: %f %f %f %f %f\n", c12[0], c12[1], c12[2], c12[3], c12[4] );
printf ( "g12: %f %f %f %f %f\n", g12[0], g12[1], g12[2], g12[3], g12[4] );
*/
  pkv_SetScratchMemTop ( sp );
  return true;

/*
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
*/
#undef r0u0
#undef r0uu0
#undef r0uuu0
#undef r0uuuu0
#undef r0s0
#undef r0us0
#undef r0uus0
#undef r0uuus0
#undef r0ss0
#undef r0uss0
#undef r0uuss0
#undef r0u1
#undef r0uu1
#undef r0s1
#undef r0us1
#undef r0uus1
#undef r0ss1
#undef r0uss1
#undef r0uuss1
#undef q0v0
#undef q0vv0
#undef q0vvv0
#undef q0vvvv0
#undef q0t0
#undef q0vt0
#undef q0vvt0
#undef q0vvvt0
#undef q0tt0
#undef q0vtt0
#undef q0vvtt0
#undef q0v1
#undef q0vv1
#undef q0t1
#undef q0vt1
#undef q0vvt1
#undef q0tt1
#undef q0vtt1
#undef q0vvtt1
#undef q1v0
#undef q1vv0
#undef q1t0
#undef q1vt0
#undef q1vvt0
#undef q1tt0
#undef q1vtt0
#undef q1vvtt0
#undef r1uu0
#undef r1s0
#undef r1uus0
#undef r1ss0
#undef r1uuss0
#undef r1uu1
#undef r1s1
#undef r1uus1
#undef r1ss1
#undef r1uuss1
#undef q1v1
#undef q1vv1
#undef q1t1
#undef q1vt1
#undef q1vvt1
#undef q1tt1
#undef q1vtt1
#undef q1vvtt1
#undef r1u0
#undef r1us0
#undef r1uss0
#undef r1u1
#undef r1us1
#undef r1uss1
} /*FindJFunctionsDif*/

static void MultJFunctionsf ( int hole_k,
               const float *b01, const float *c01,
               float *b01b01, float *twob01c01, float *c01c01,
               const float *b11, const float *c11,
               float *b11b11, float *twob11c11, float *c11c11 )
{
/* compute the products of polynomials, which will be used many times later */
  int degprod;

  mbs_multiMultBezCf ( hole_k, G2_BF01DEG, G2_BF01DEG+1, b01,
                    1, hole_k, G2_BF01DEG, G2_BF01DEG+1, b01,
                    &degprod, 2*G2_BF01DEG+1, b01b01 );
  mbs_multiMultBezCf ( hole_k, G2_BF01DEG, G2_BF01DEG+1, b01,
                    1, hole_k, G2_CG01DEG, G2_CG01DEG+1, c01,
                    &degprod, G2_BF01DEG+G2_CG01DEG+1, twob01c01 );
  pkn_MultMatrixNumf ( 1, hole_k*(G2_BF01DEG+G2_CG01DEG+1), 0, twob01c01, 2.0,
                       0, twob01c01 );
  mbs_multiMultBezCf ( hole_k, G2_CG01DEG, G2_CG01DEG+1, c01,
                    1, hole_k, G2_CG01DEG, G2_CG01DEG+1, c01,
                    &degprod, 2*G2_CG01DEG+1, c01c01 );

  mbs_multiMultBezCf ( hole_k, G2_BF11DEG, G2_BF11DEG+1, b11,
                    1, hole_k, G2_BF11DEG, G2_BF11DEG+1, b11,
                    &degprod, 2*G2_BF11DEG+1, b11b11 );
  mbs_multiMultBezCf ( hole_k, G2_BF11DEG, G2_BF11DEG+1, b11,
                    1, hole_k, G2_CG11DEG, G2_CG11DEG+1, c11,
                    &degprod, G2_BF11DEG+G2_CG11DEG+1, twob11c11 );
  pkn_MultMatrixNumf ( 1, hole_k*(G2_BF11DEG+G2_CG11DEG+1), 0, twob11c11, 2.0,
                       0, twob11c11 );
  mbs_multiMultBezCf ( hole_k, G2_CG11DEG, G2_CG11DEG+1, c11,
                    1, hole_k, G2_CG11DEG, G2_CG11DEG+1, c11,
                    &degprod, 2*G2_CG11DEG+1, c11c11 );
} /*MultJFunctionsf*/

static boolean FindJFunctionsf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf *privateG;
  G2HolePrivateRecf *privateG2;
  int      hole_k;
  int      i, j;
  vector2f *omcbc, *omcbcd, *omcbcdd, *surrpcbc;
  float    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
           *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
           *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01,
           *b11b11, *twob11c11, *c11c11, *f11f11, *twof11g11, *g11g11;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  i = 2*(4*(G2_BF01DEG+G2_CG01DEG+G2_BF11DEG+G2_CG11DEG)+
         G2_BF02DEG+G2_CG02DEG+G2_BF12DEG+G2_CG12DEG);
  privateG2->jfunc = malloc ( hole_k*(28+i)*sizeof(float) );
  if ( !privateG2->jfunc ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    goto failure;
  }
  G2GetPolynomialAddresses ( privateG2->jfunc,
      b01, c01, f01, g01, b11, c11, f11, g11,
      b02, c02, f02, g02, b12, c12, f12, g12, b01b01, twob01c01, c01c01,
      f01f01, twof01g01, g01g01, b11b11, twob11c11, c11c11,
      f11f11, twof11g11, g11g11 );

  omcbc = privateG2->omcbc;
  omcbcd = &omcbc[(G2H_OMCDEG+1)*hole_k];
  omcbcdd = &omcbcd[G2H_OMCDEG*hole_k];
  surrpcbc = privateG->surrpcbc;

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    if ( !FindJFunctionsDif (
            &omcbc[(G2H_OMCDEG+1)*i], &omcbcd[G2H_OMCDEG*i], &omcbcdd[(G2H_OMCDEG-1)*i],
            &omcbc[(G2H_OMCDEG+1)*j], &omcbcd[G2H_OMCDEG*j], &omcbcdd[(G2H_OMCDEG-1)*j],
            &surrpcbc[18*(2*i+1)], &surrpcbc[36*j],
            &b01[i*(G2_BF01DEG+1)], &c01[i*(G2_CG01DEG+1)],
            &f01[i*(G2_BF01DEG+1)], &g01[i*(G2_CG01DEG+1)],
            &b11[i*(G2_BF11DEG+1)], &c11[i*(G2_CG11DEG+1)],
            &f11[i*(G2_BF11DEG+1)], &g11[i*(G2_CG11DEG+1)],
            &b02[i*(G2_BF02DEG+1)], &c02[i*(G2_CG02DEG+1)],
            &f02[i*(G2_BF02DEG+1)], &g02[i*(G2_CG02DEG+1)],
            &b12[i*(G2_BF12DEG+1)], &c12[i*(G2_CG12DEG+1)],
            &f12[i*(G2_BF12DEG+1)], &g12[i*(G2_CG12DEG+1)] ) )
      goto failure;
  }
  if ( !_g2h_VerifyJunctionFunctionsf ( domain ) ) {
    domain->error_code = G2H_ERROR_INVALID_JUNC_FUNC;
    goto failure;
  }
  MultJFunctionsf ( hole_k, b01, c01, b01b01, twob01c01, c01c01,
                            b11, c11, b11b11, twob11c11, c11c11 );
  MultJFunctionsf ( hole_k, f01, g01, f01f01, twof01g01, g01g01,
                            f11, g11, f11f11, twof11g11, g11g11 );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindJFunctionsf*/

