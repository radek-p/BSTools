
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindJFunctionsDif (
     const vector2f *omcbci, const vector2f *omcbcdi,
     const vector2f *omcbcj, const vector2f *omcbcdj,
     const vector2f *surrpcbci, const vector2f *surrpcbcj,
     float *b01, float *c01, float *f01, float *g01,
     float *b11, float *c11, float *f11, float *g11 )
{
/* the following #definitions are introduced for convenient accessing */
/* the partial derivatives of auxiliary and surrounding patches */
#define r0u0    omcbci[1]
#define r0uu0   omcbci[2]
#define r0s0    omcbcdi[0]
#define r0us0   omcbcdi[1]
#define r0u1    omcbci[4]
#define r0s1    omcbcdi[2]
#define r0us1   omcbcdi[3]
#define q0v0    omcbcj[1]
#define q0vv0   omcbcj[2]
#define q0t0    omcbcdj[0]
#define q0vt0   omcbcdj[1]
#define q0v1    omcbcj[4]
#define q0t1    omcbcdj[2]
#define q0vt1   omcbcdj[3]
#define q1v0    surrpcbci[3]
#define q1t0    surrpcbci[1]
#define q1vv0   surrpcbci[6]
#define q1vt0   surrpcbci[4]
#define q1v1    surrpcbci[9+3]
#define q1t1    surrpcbci[9+1]
#define q1vt1   surrpcbci[9+4]
#define r1u0    surrpcbcj[3]
#define r1s0    surrpcbcj[1]
#define r1uu0   surrpcbcj[6]
#define r1us0   surrpcbcj[4]
#define r1u1    surrpcbcj[9+3]
#define r1s1    surrpcbcj[9+1]
#define r1us1   surrpcbcj[9+4]

  void     *sp;
  float    b010, b01d0, c010, c01d0, b011, b01d1, c011, c01d1,
           f010, f01d0, g010, g01d0, f011, f01d1, g011, g01d1,
           b110, b11d0, c110, c11d0, b111, b11d1, c111, c11d1,
           f110, f11d0, g110, g11d0, f111, f11d1, g111, g11d1;

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

          /* by assumption, c01 and g01 are polynomials of degree 1, */
          /* hence their derivatives are the following: */
  c01d0 = c01d1 = c011-c010;
  g01d0 = g01d1 = g011-g010;

          /* the values of derivatives of b01 and f01 at 0 are now */
          /* determined by the mixed partial derivatives compatibility */
          /* condition */
  SolveCompatibilityEq2af (
        &r0u0, &r0s0, &r0uu0, &r0us0, &q0v0, &q0t0, &q0vv0, &q0vt0,
        b010, c010, f010, g010, &b01d0, c01d0, &f01d0, g01d0 );

          /* by assumption, b01 and f01 are polynomials of degree 2, */
          /* hence their derivatives are the following: */
  b01d1 = (float)(2.0*(b011-b010)-b01d0);
  f01d1 = (float)(2.0*(f011-f010)-f01d0);

          /* now derivatives of b11, c11 and f11, g11 at 0 are determined */
          /* by mixed partial derivatives compatibility conditions */
          /* at (0,1) and (1,0) respectively */
  SolveCompatibilityEq2df ( &q0v1, &q0t1, &q0vt1, &r1u0, &r1s0, &r1us0,
                            g011, c110, f01d1, g01d1, &b11d0, &c11d0 );
  SolveCompatibilityEq2df ( &r0u1, &r0s1, &r0us1, &q1v0, &q1t0, &q1vt0,
                            c011, g110, b01d1, c01d1, &f11d0, &g11d0 );

          /* by assumption, c11 and g11 are polynomials of degree 2, */
          /* hence their derivatives are the following: */
  c11d1  = (float)(2.0*(c111-c110)-c11d0);
  g11d1  = (float)(2.0*(g111-g110)-g11d0);

          /* the derivatives of b11 and f11 at 1 must satisfy the */
          /* mixed partial derivatives compatibility condition at (1,1) */
  SolveCompatibilityEq2cf ( &r1u1, &r1s1, &r1us1, &q1v1, &q1t1, &q1vt1,
                            c111, g111, &b11d1, c11d1, &f11d1, g11d1 );

        /* now compute the coefficients of the polynomials in Bernstein bases */
  b01[0] = b010;  b01[1] = (float)(b010+0.5*b01d0);  b01[2] = b011;
  c01[0] = c010;  c01[1] = c011;
  f01[0] = f010;  f01[1] = (float)(f010+0.5*f01d0);  f01[2] = f011;
  g01[0] = g010;  g01[1] = g011;
  b11[0] = b110;  b11[1] = (float)(b110+b11d0/3.0);
  b11[2] = (float)(b111-b11d1/3.0);  b11[3] = b111;
  c11[0] = c110;  c11[1] = (float)(c110+0.5*c11d0);  c11[2] = c111;
  f11[0] = f110;  f11[1] = (float)(f110+f11d0/3.0);
  f11[2] = (float)(f111-f11d1/3.0);  f11[3] = f111;
  g11[0] = g110;  g11[1] = (float)(g110+0.5*g11d0);  g11[2] = g111;

  pkv_SetScratchMemTop ( sp );
  return true;

/*
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
*/
#undef r0u0
#undef r0uu0
#undef r0s0
#undef r0us0
#undef r0u1
#undef r0s1
#undef r0us1
#undef q0v0
#undef q0vv0
#undef q0t0
#undef q0vt0
#undef q0v1
#undef q0t1
#undef q0vt1
#undef q1v0
#undef q1t0
#undef q1vt0
#undef r1s0
#undef r1s1
#undef q1v1
#undef q1t1
#undef q1vt1
#undef r1u0
#undef r1us0
#undef r1u1
#undef r1us1
} /*FindJFunctionsDif*/

static boolean FindJFunctionsf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf *privateG;
  G1HolePrivateRecf *privateG1;
  int      hole_k;
  int      i, j;
  vector2f *omcbc, *omcbcd, *surrpcbc;
  float    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  i = 2*(G1_BF01DEG+G1_CG01DEG+G1_BF11DEG+G1_CG11DEG+4);
  privateG1->jfunc = malloc ( hole_k*i*sizeof(float) );
  if ( !privateG1->jfunc ) {
    domain->error_code = G1H_ERROR_CANNOT_MALLOC;
    goto failure;
  }
  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );

  omcbc = privateG1->omcbc;
  omcbcd = &omcbc[(G1H_OMCDEG+1)*hole_k];
  surrpcbc = privateG->surrpcbc;

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    if ( !FindJFunctionsDif ( &omcbc[(G1H_OMCDEG+1)*i], &omcbcd[G1H_OMCDEG*i],
                              &omcbc[(G1H_OMCDEG+1)*j], &omcbcd[G1H_OMCDEG*j],
                              &surrpcbc[18*(2*i+1)], &surrpcbc[36*j],
                              &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
                              &f01[i*(G1_BF01DEG+1)], &g01[i*(G1_CG01DEG+1)],
                              &b11[i*(G1_BF11DEG+1)], &c11[i*(G1_CG11DEG+1)],
                              &f11[i*(G1_BF11DEG+1)], &g11[i*(G1_CG11DEG+1)] ) )
      goto failure;
  }
  if ( !_g1h_VerifyJunctionFunctionsf ( domain ) ) {
    domain->error_code = G1H_ERROR_INVALID_JUNC_FUNC;
    goto failure;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindJFunctionsf*/

