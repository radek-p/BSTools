
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean AnalyzePartitionf ( GHoleDomainf *domain )
{
  G2HolePrivateRecf *privateG2;
  int            auxp_k;
  unsigned char  *spdimen;

  privateG2 = domain->privateG2;
  spdimen = &privateG2->spdimen[0];

  if ( _gh_AnalyzePartitionf ( domain, G2H_OMCDEG, privateG2->omcbc,
                      &privateG2->hole_m, &auxp_k, &privateG2->partition,
                      &privateG2->spartition, &privateG2->spart_alpha0 ) ) {
    spdimen[0] = (unsigned char)privateG2->hole_m;     /* cubic half-polynomials */
    spdimen[1] = (unsigned char)max ( 0, auxp_k-4 );   /* cubic B-splines */
    spdimen[2] = (unsigned char)(2*privateG2->hole_m); /* quartic half-polynomials */
    spdimen[3] = (unsigned char)max ( 0, 2*auxp_k-5 ); /* quartic B-splines */

    privateG2->nfunc_a = 15 + 3*privateG2->hole_m + spdimen[1] + spdimen[3];
/*
printf ( "spdimen: %d %d %d %d %d\n",
         spdimen[0], spdimen[1], spdimen[2], spdimen[3], privateG2->nfunc_a );
*/
    return true;
  }
  else
    return false;
} /*AnalyzePartitionf*/

static void ComputeAuxiMat ( GHoleDomainf *domain )
{
  G2HolePrivateRecf *privateG2;
  int      i, hole_k;
  vector2f diu, div, diuu, diuv, divv, diuuu, diuuv, diuvv,
           diuuuu, diuuuv, diuuvv;
  vector2f *omcbc, *omcbcd, *omcbcdd;
  float    *A11, *A21, *A22, *A31, *A32, *A33, *A41, *A42, *A43, *A44;

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  omcbc = privateG2->omcbc;
  omcbcd = &omcbc[(G2H_OMCDEG+1)*hole_k];
  omcbcdd = &omcbcd[G2H_OMCDEG*hole_k];

  for ( i = 0; i < hole_k; i++ ) {

        /* setup the addresses of the transformation matrices */
    A11 = &privateG2->AuxiMat[i*125];
    A21 = &A11[4];   A22 = &A21[6];
    A31 = &A22[9];   A32 = &A31[8];   A33 = &A32[12];
    A41 = &A33[16];  A42 = &A41[10];  A43 = &A42[15];  A44 = &A43[20];

        /* extract the derivatives of the domain patches */
    diu = omcbc[(G2H_OMCDEG+1)*i+1];     div = omcbcd[G2H_OMCDEG*i];
    diuu = omcbc[(G2H_OMCDEG+1)*i+2];    diuv = omcbcd[G2H_OMCDEG*i+1];    divv = omcbcdd[(G2H_OMCDEG-1)*i];
    diuuu = omcbc[(G2H_OMCDEG+1)*i+3];   diuuv = omcbcd[G2H_OMCDEG*i+2];   diuvv = omcbcdd[(G2H_OMCDEG-1)*i+1];
    diuuuu = omcbc[(G2H_OMCDEG+1)*i+4];  diuuuv = omcbcd[G2H_OMCDEG*i+3];  diuuvv = omcbcdd[(G2H_OMCDEG-1)*i+2];

        /* compute the derivative transformation matrices */
    pkn_Setup2DerA11Matrixf ( diu.x, diu.y, div.x, div.y, A11 );
    pkn_Setup2DerA21Matrixf ( diuu.x, diuu.y, diuv.x,
                              diuv.y, divv.x, divv.y, A21 );
    pkn_Setup2DerA22Matrixf ( diu.x, diu.y, div.x, div.y, A22 );
    pkn_Setup2DerA31Matrixf ( diuuu.x, diuuu.y, diuuv.y, diuuv.y,
                              diuvv.x, diuvv.y, 0.0, 0.0, A31 );
    pkn_Setup2DerA32Matrixf ( diu.x, diu.y, div.x, div.y,
                              diuu.x, diuu.y, diuv.x,
                              diuv.y, divv.x, divv.y, A32 );
    pkn_Setup2DerA33Matrixf ( diu.x, diu.y, div.x, div.y, A33 );
    pkn_Setup2DerA41Matrixf ( diuuuu.x, diuuuu.y, diuuuv.x, diuuuv.y,
                              diuuvv.x, diuuvv.y, 0.0, 0.0, 0.0, 0.0, A41 );
    pkn_Setup2DerA42Matrixf ( diu.x, diu.y, div.x, div.y,
                              diuu.x, diuu.y, diuv.x,
                              diuv.y, divv.x, divv.y,
                              diuuu.x, diuuu.y, diuuv.y, diuuv.y,
                              diuvv.x, diuvv.y, 0.0, 0.0, A42 );
    pkn_Setup2DerA43Matrixf ( diu.x, diu.y, div.x, div.y,
                              diuu.x, diuu.y, diuv.x,
                              diuv.y, divv.x, divv.y, A43 );
    pkn_Setup2DerA44Matrixf ( diu.x, diu.y, div.x, div.y, A44 );
  }
} /*ComputeAuxiMat*/

/* ///////////////////////////////////////////////////////////////////////// */
static void FindDivDiff4Coeff ( float a, float b, float c, float d, float e,
                                float *dd )
{
    /* compute the coefficients of the divided difference of 4th order  */
    /* for single knots a,b,c,d,e, i.e. the numbers such that           */
    /* f[a,b,c,d,e] = dd[0]f(a)+dd[1]f(b)+dd[2]f(c)+dd[3]f(d)+dd[4]f(e) */
  dd[0] = (float)(1.0/((a-b)*(a-c)*(a-d)*(a-e)));
  dd[1] = (float)(1.0/((b-a)*(b-c)*(b-d)*(b-e)));
  dd[2] = (float)(1.0/((c-a)*(c-b)*(c-d)*(c-e)));
  dd[3] = (float)(1.0/((d-a)*(d-b)*(d-c)*(d-e)));
  dd[4] = (float)(1.0/((e-a)*(e-b)*(e-c)*(e-d)));
} /*FindDivDiff4Coeff*/

static void FindDivDiff5Coeff1 ( float a, float b, float c,
                                 float *dd )
{
    /* compute the coefficients of the divided difference of 5th order  */
    /* for double knots a,b,c, i.e. the numbers such that               */
    /* f[a,a,b,b,c,c] = dd[0]f(a)+dd[1]f'(a)+dd[2]f(b)+                 */
    /*                  dd[3]f'(b)+dd[4]f(c)+dd[5]f'(c).                */
  float ba, ca, cb;

  ba = (float)(1.0/(b-a));  ca = (float)(1.0/(c-a));  cb = (float)(1.0/(c-b));

  dd[1] = ba*ba*ca*ca;
  dd[0] = (float)(2.0*(ca+ba)*dd[1]);
  dd[3] = cb*cb*ba*ba;
  dd[2] = (float)(2.0*(cb-ba)*dd[3]);
  dd[5] = cb*cb*ca*ca;
  dd[4] = (float)(-2.0*(cb+ca)*dd[5]);
} /*FindDivDiff5Coeff1*/

static void FindDivDiff5Coeff2 ( float a, float b, float c, float d,
                                 float *dd )
{
    /* compute the coefficients of the divided difference of 5th order  */
    /* for knots a,b,b,c,c,d i.e. the numbers such that                 */
    /* f[a,b,b,c,c,d] = dd[0]f(a)+dd[1]f(b)+dd[2]f'(b)+                 */
    /*                  dd[3]f(c)+dd[4]f'(c)+dd[5]f(d).                 */
  float ba, ca, cb, da, db, dc;

  ba = (float)(1.0/(b-a));  ca = (float)(1.0/(c-a));  cb = (float)(1.0/(c-b));
  da = (float)(1.0/(d-a));  db = (float)(1.0/(d-b));  dc = (float)(1.0/(d-c));

  dd[0] = -ba*ba*ca*ca*da;
  dd[1] = (float)(-(db+2.0*cb-ba)*ba*db*cb*cb);
  dd[2] = -ba*db*cb*cb;
  dd[3] = (float)(-(dc-2.0*cb-ca)*dc*ca*cb*cb);
  dd[4] = -dc*ca*cb*cb;
  dd[5] = da*db*db*dc*dc;
} /*FindDivDiff5Coeff2*/

static void ReparHomoPoly1Der ( float *p1, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc, float *br0cr2bc )
{
  float q1[2], q2[3], q3[4], q4[5];
  float *A11, *A21, *A31, *A41;

  A11 = &AuxiMat[i*125];  A21 = &A11[4];  A31 = &A21[6+9];  A41 = &A31[8+12+16];
  pkn_MultMatrixf ( 2, 2, 2, A11, 1, 1, p1, 1, q1 );
  pkn_MultMatrixf ( 3, 2, 2, A21, 1, 1, p1, 1, q2 );
  pkn_MultMatrixf ( 4, 2, 2, A31, 1, 1, p1, 1, q3 );
  pkn_MultMatrixf ( 5, 2, 2, A41, 1, 1, p1, 1, q4 );
  br0bc[(G2H_OMCDEG+1)*i+1]    = q1[0];
  br0bc[(G2H_OMCDEG+1)*i+2]    = q2[0];
  br0bc[(G2H_OMCDEG+1)*i+3]    = q3[0];
  br0bc[(G2H_OMCDEG+1)*i+4]    = q4[0];
  br0cr1bc[G2H_OMCDEG*i]   = q1[1];
  br0cr1bc[G2H_OMCDEG*i+1] = q2[1];
  br0cr1bc[G2H_OMCDEG*i+2] = q3[1];
  br0cr1bc[G2H_OMCDEG*i+3] = q4[1];
  br0cr2bc[(G2H_OMCDEG-1)*i]   = q2[2];
  br0cr2bc[(G2H_OMCDEG-1)*i+1] = q3[2];
  br0cr2bc[(G2H_OMCDEG-1)*i+2] = q4[2];
} /*ReparHomoPoly1Der*/

static void ReparHomoPoly2Der ( float *p2, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc, float *br0cr2bc )
{
  float q2[3], q3[4], q4[5];
  float *A22, *A32, *A42;

  A22 = &AuxiMat[i*125+4+6];  A32 = &A22[9+8];  A42 = &A32[12+16+10];
  pkn_MultMatrixf ( 3, 3, 3, A22, 1, 1, p2, 1, q2 );
  pkn_MultMatrixf ( 4, 3, 3, A32, 1, 1, p2, 1, q3 );
  pkn_MultMatrixf ( 5, 3, 3, A42, 1, 1, p2, 1, q4 );
  br0bc[(G2H_OMCDEG+1)*i+2]    = q2[0];
  br0bc[(G2H_OMCDEG+1)*i+3]    = q3[0];
  br0bc[(G2H_OMCDEG+1)*i+4]    = q4[0];
  br0cr1bc[G2H_OMCDEG*i+1] = q2[1];
  br0cr1bc[G2H_OMCDEG*i+2] = q3[1];
  br0cr1bc[G2H_OMCDEG*i+3] = q4[1];
  br0cr2bc[(G2H_OMCDEG-1)*i]   = q2[2];
  br0cr2bc[(G2H_OMCDEG-1)*i+1] = q3[2];
  br0cr2bc[(G2H_OMCDEG-1)*i+2] = q4[2];
} /*ReparHomoPoly2Der*/

static void ReparHomoPoly3Der ( float *p3, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc, float *br0cr2bc,
                                char sw )
{
  float q3[4], q4[5];
  float *A33, *A43;

  A33 = &AuxiMat[i*125+4+6+9+8+12];  A43 = &A33[16+10+15];
  pkn_MultMatrixf ( 4, 4, 4, A33, 1, 1, p3, 1, q3 );
  pkn_MultMatrixf ( 5, 4, 4, A43, 1, 1, p3, 1, q4 );
  switch ( sw ) {
case 0:
    br0bc[(G2H_OMCDEG+1)*i+3]    = q3[0];
    br0bc[(G2H_OMCDEG+1)*i+4]    = q4[0];
    br0cr1bc[G2H_OMCDEG*i+2] = q3[1];
    br0cr1bc[G2H_OMCDEG*i+3] = q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+1] = q3[2];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] = q4[2];
    break;
case 1:
    br0bc[(G2H_OMCDEG+1)*i+3]    += q3[0];
    br0bc[(G2H_OMCDEG+1)*i+4]    += q4[0];
    br0cr1bc[G2H_OMCDEG*i+2] += q3[1];
    br0cr1bc[G2H_OMCDEG*i+3] += q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+1] += q3[2];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] += q4[2];
    break;
case 2:
    br0bc[(G2H_OMCDEG+1)*i+3]    -= q3[0];
    br0bc[(G2H_OMCDEG+1)*i+4]    -= q4[0];
    br0cr1bc[G2H_OMCDEG*i+2] -= q3[1];
    br0cr1bc[G2H_OMCDEG*i+3] -= q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+1] -= q3[2];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] -= q4[2];
    break;
  }
} /*ReparHomoPoly3Der*/

static void ReparHomoPoly4Der ( float *p4, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc, float *br0cr2bc,
                                char sw )
{
  float q4[5];
  float *A44;

  A44 = &AuxiMat[i*125+4+6+9+8+12+16+10+15+20];
  pkn_MultMatrixf ( 5, 5, 5, A44, 1, 1, p4, 1, q4 );
  switch ( sw ) {
case 0:
    br0bc[(G2H_OMCDEG+1)*i+4]    = q4[0];
    br0cr1bc[G2H_OMCDEG*i+3] = q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] = q4[2];
    break;
case 1:
    br0bc[(G2H_OMCDEG+1)*i+4]    += q4[0];
    br0cr1bc[G2H_OMCDEG*i+3] += q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] += q4[2];
    break;
case 2:
    br0bc[(G2H_OMCDEG+1)*i+4]    -= q4[0];
    br0cr1bc[G2H_OMCDEG*i+3] -= q4[1];
    br0cr2bc[(G2H_OMCDEG-1)*i+2] -= q4[2];
    break;
  }
} /*ReparHomoPoly4Der*/

boolean _g2h_GetABasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *br0, float *br0cr1, float *br0cr2 )
{
#define TOL 1.0e-6
  void     *sp;
  G2HolePrivateRecf *privateG2;
  int      hole_k;
  unsigned char *spdimen;
  float    *partition;
  GHoleSgnPartf *spartition;
  float    *br0bc, *br0cr1bc, *br0cr2bc;
  float    p1[2], p2[3], p3[4], p4[5];
  float    alpha0, sa0, ca0, sa, sb, ca, cb, a, b;
  float    u, n, m, normf;
  float    ddc[6];
  int      sfn;
  int      i, j, k;


  sp = pkv_GetScratchMemTop ();

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  spdimen = privateG2->spdimen;

  br0bc = pkv_GetScratchMemf ( 3*G2H_OMCDEG*hole_k );
  if ( !br0bc )
    goto failure;

  br0cr1bc = &br0bc[hole_k*(G2H_OMCDEG+1)];   br0cr2bc = &br0cr1bc[hole_k*G2H_OMCDEG];

  memset ( br0bc, 0, hole_k*3*G2_CROSS00DEG*sizeof(float) );
  switch ( fn ) {
case 0:    /* f(x,y) = 1 */
    for ( i = 0; i < hole_k; i++ )
      br0bc[(G2H_OMCDEG+1)*i] = 1.0;
    break;

case 1:    /* f(x,y) = x */
    p1[0] = 1.0;  p1[1] = 0.0;
    goto continue_deg1;

case 2:    /* f(x,y) = y */
    p1[0] = 0.0;  p1[1] = 1.0;
continue_deg1:
    for ( i = 0; i < hole_k;  i++ )
      ReparHomoPoly1Der ( p1, i, privateG2->AuxiMat, br0bc, br0cr1bc, br0cr2bc );
    break;

case 3:    /* x^2 */
    p2[0] = 2.0;  p2[1] = p2[2] = 0.0;
    goto continue_deg2;

case 4:    /* xy */
    p2[0] = 0.0;  p2[1] = 1.0;  p2[2] = 0.0;
    goto continue_deg2;

case 5:    /* y^2 */
    p2[0] = p2[1] = 0.0;  p2[2] = 2.0;
continue_deg2:
    for ( i = 0;  i < hole_k;  i++ )
      ReparHomoPoly2Der ( p2, i, privateG2->AuxiMat, br0bc, br0cr1bc, br0cr2bc );
    break;

case 6:    /* x^3 */
    p3[0] = 6.0;  p3[1] = p3[2] = p3[3] = 0.0;
    goto continue_deg3;

case 7:    /* x^2y */
    p3[0] = 0.0;  p3[1] = 2.0;  p3[2] = p3[3] = 0.0;
    goto continue_deg3;

case 8:    /* xy^2 */
    p3[0] = p3[1] = 0.0;  p3[2] = 2.0;  p3[3] = 0.0;
    goto continue_deg3;

case 9:    /* y^3 */
    p3[0] = p3[1] = p3[2] = 0.0;  p3[3] = 6.0;
continue_deg3:
    for ( i = 0;  i < hole_k;  i++ )
      ReparHomoPoly3Der ( p3, i, privateG2->AuxiMat,
                          br0bc, br0cr1bc, br0cr2bc, 0 );
    break;

case 10:   /* x^4 */
    p4[0] = 24.0;  p4[1] = p4[2] = p4[3] = p4[4] = 0.0;
    goto continue_deg4;

case 11:   /* x^3y */
    p4[0] = 0.0;  p4[1] = 6.0;  p4[2] = p4[3] = p4[4] = 0.0;
    goto continue_deg4;

case 12:   /* x^2y^2 */
    p4[0] = p4[1] = 0.0;  p4[2] = 4.0;  p4[3] = p4[4] = 0.0;
    goto continue_deg4;

case 13:   /* xy^3 */
    p4[0] = p4[1] = p4[2] = 0.0;  p4[3] = 6.0;  p4[4] = 0.0;
    goto continue_deg4;

case 14:   /* y^4 */
    p4[0] = p4[1] = p4[2] = p4[3] = 0.0;  p4[4] = 24.0;
continue_deg4:
    for ( i = 0;  i < hole_k;  i++ )
      ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                          br0bc, br0cr1bc, br0cr2bc, 0 );
    break;

default:   /* splines */
    partition = privateG2->partition;
    spartition = privateG2->spartition;
    alpha0 = privateG2->spart_alpha0;
    sfn = fn-15;        /* which one after the polynomials? */

    if ( sfn < spdimen[0] ) {                 /* a cubic half-polynomial */
      for ( j = k = -1;  k < sfn;  j++ )
        if ( spartition[j+1].both ) k++;
      sa = (float)sin ( spartition[j].malpha );
      ca = (float)cos ( spartition[j].malpha );
      p3[0] = (float)(-6.0*sa*sa*sa);  p3[1] = (float)(6.0*sa*sa*ca);
      p3[2] = (float)(-6.0*sa*ca*ca);  p3[3] = (float)(6.0*ca*ca*ca);
      for ( i = 0; i < hole_k; i++ ) {
        sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
        if ( sa*cb-ca*sb > TOL )
          ReparHomoPoly3Der ( p3, i, privateG2->AuxiMat,
                              br0bc, br0cr1bc, br0cr2bc, 0 );
      }
    }
    else if ( sfn-spdimen[0] < spdimen[1] ) { /* a cubic B-spline */
      sfn -= spdimen[0];
      sa0 = (float)sin ( alpha0 );  ca0 = (float)cos ( alpha0 );
      FindDivDiff4Coeff ( spartition[sfn+4].knot, spartition[sfn+3].knot,
                          spartition[sfn+2].knot, spartition[sfn+1].knot,
                          spartition[sfn].knot, ddc );
      normf = spartition[sfn].knot-spartition[sfn+4].knot;
      for ( j = 0; j < 5; j++ )
        ddc[j] *= normf;

      for ( j = 0; j < 5; j++ ) {
        sa = (float)sin ( spartition[sfn+4-j].malpha );
        ca = (float)cos ( spartition[sfn+4-j].malpha );
        u = spartition[sfn+4-j].knot;  a = sa0-u*ca0;  b = -(ca0+u*sa0);
            /* compute the partial derivatives of the third order */
            /* of the polynomial (X-uY)^3 = (ax+by)^3, times ddc[j] */
        p3[0] = (float)(6.0*ddc[j]*a*a*a);
        p3[1] = (float)(6.0*ddc[j]*a*a*b);
        p3[2] = (float)(6.0*ddc[j]*a*b*b);
        p3[3] = (float)(6.0*ddc[j]*b*b*b);

        for ( i = 0; i < hole_k; i++ ) {
          sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
          n = sa*cb-ca*sb;
          if ( n > TOL ) {
            m = ca0*cb+sa0*sb;
            if ( m > TOL ) {
              if ( !spartition[sfn+4-j].sgn || spartition[sfn+4-j].both )
                ReparHomoPoly3Der ( p3, i, privateG2->AuxiMat,
                                    br0bc, br0cr1bc, br0cr2bc, 1 );
            }
            else if ( m < -TOL ) {
              if ( spartition[sfn+4-j].sgn && !spartition[sfn+4-j].both )
                ReparHomoPoly3Der ( p3, i, privateG2->AuxiMat,
                                    br0bc, br0cr1bc, br0cr2bc, 2 );
            }
          }
        }
      }
    }
    else if ( sfn-spdimen[0]-spdimen[1] < spdimen[2] ) {
                                              /* a quartic half-polynomial */
      sfn -= spdimen[0]+spdimen[1];
      for ( j = k = -1; k < sfn/2; j++ )
        if ( spartition[j+1].both ) k++;
      sa = (float)sin ( spartition[j].malpha );
      ca = (float)cos ( spartition[j].malpha );
      if ( !(sfn & 0x0001) ) {
        p4[0] = (float)(24.0*sa*sa*sa*sa);  p4[1] = (float)(-24.0*ca*sa*sa*sa);
        p4[2] = (float)(24.0*ca*ca*sa*sa);  p4[3] = (float)(-24.0*ca*ca*ca*sa);
        p4[4] = (float)(24.0*ca*ca*ca*ca);
      }
      else {
        p4[0] = (float)(-24.0*ca*sa*sa*sa);
        p4[1] = (float)( 6.0*sa*sa*(3.0*ca*ca-sa*sa));
        p4[2] = (float)( 12.0*ca*sa*(sa*sa-ca*ca));
        p4[3] = (float)( 6.0*ca*ca*(ca*ca-3.0*sa*sa));
        p4[4] = (float)( 24.0*ca*ca*ca*sa);
      }
      for ( i = 0; i < hole_k; i++ ) {
        sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
        if ( sa*cb-ca*sb > 0.0 )
          ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                              br0bc, br0cr1bc, br0cr2bc, 0 );
      }
    }
    else {                                    /* a quartic B-spline */
      sfn -= spdimen[0]+spdimen[1]+spdimen[2];
      sa0 = (float)sin ( alpha0 );  ca0 = (float)cos ( alpha0 );
      if ( !(sfn & 0x0001) ) {
        sfn = sfn >> 1;
        FindDivDiff5Coeff1 ( spartition[sfn+2].knot, spartition[sfn+1].knot,
                             spartition[sfn].knot, ddc );
        normf = spartition[sfn+2].knot-spartition[sfn].knot;
        for ( j = 0; j < 6; j++ )
          ddc[j] *= normf;

        for ( j = 0; j < 3; j++ ) {
          sa = (float)sin ( spartition[sfn+2-j].malpha );
          ca = (float)cos ( spartition[sfn+2-j].malpha );
          u = spartition[sfn+2-j].knot;  a = sa0-u*ca0;  b = -(ca0+u*sa0);
              /* compute the partial derivatives of the fourth order */
              /* of the polynomial (X-uY)^4 = (ax+by)^4, times ddc[2j] */
              /* and of -4Y(X-uY)^3 = (cx+sy)(ax+by)^3, times ddc[2j+1] */
          p4[0] = (float)(24.0*ddc[j+j]*a*a*a*a - 96.0*ddc[j+j+1]*a*a*a*ca0);
          p4[1] = (float)(24.0*ddc[j+j]*a*a*a*b - 24.0*ddc[j+j+1]*a*a*(a*sa0+3.0*b*ca0));
          p4[2] = (float)(24.0*ddc[j+j]*a*a*b*b - 48.0*ddc[j+j+1]*a*b*(a*sa0+b*ca0));
          p4[3] = (float)(24.0*ddc[j+j]*a*b*b*b - 24.0*ddc[j+j+1]*b*b*(3.0*a*sa0+b*ca0));
          p4[4] = (float)(24.0*ddc[j+j]*b*b*b*b - 96.0*ddc[j+j+1]*b*b*b*sa0);

          for ( i = 0; i < hole_k; i++ ) {
            sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
            n = sa*cb-ca*sb;
            if ( n > TOL ) {
              m = ca0*cb+sa0*sb;
              if ( m > TOL ) {
                if ( !spartition[sfn+2-j].sgn || spartition[sfn+2-j].both )
                  ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                                      br0bc, br0cr1bc, br0cr2bc, 1 );
              }
              else if ( m < -TOL ) {
                if ( spartition[sfn+2-j].sgn && !spartition[sfn+2-j].both )
                  ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                                      br0bc, br0cr1bc, br0cr2bc, 2 );
              }
            }
          }
        }
      }
      else {
        sfn = sfn >> 1;
        FindDivDiff5Coeff2 ( spartition[sfn+3].knot, spartition[sfn+2].knot,
                             spartition[sfn+1].knot, spartition[sfn].knot, ddc );
        normf = spartition[sfn+3].knot-spartition[sfn].knot;
        for ( j = 0; j < 6; j++ )
          ddc[j] *= normf;

        for ( j = 0; j < 4; j++ ) {
          sa = (float)sin ( spartition[sfn+3-j].malpha );
          ca = (float)cos ( spartition[sfn+3-j].malpha );
          u = spartition[sfn+3-j].knot;  a = sa0-u*ca0;  b = -(ca0+u*sa0);
          switch ( j ) {
        case 0: k = 0;  goto SINGLE;
        case 1: k = 1;  goto DOUBLE;
        case 2: k = 3;
DOUBLE:     p4[0] = (float)(24.0*ddc[k]*a*a*a*a - 96.0*ddc[k+1]*a*a*a*ca0);
            p4[1] = (float)(24.0*ddc[k]*a*a*a*b - 24.0*ddc[k+1]*a*a*(a*sa0+3.0*b*ca0));
            p4[2] = (float)(24.0*ddc[k]*a*a*b*b - 48.0*ddc[k+1]*a*b*(a*sa0+b*ca0));
            p4[3] = (float)(24.0*ddc[k]*a*b*b*b - 24.0*ddc[k+1]*b*b*(3.0*a*sa0+b*ca0));
            p4[4] = (float)(24.0*ddc[k]*b*b*b*b - 96.0*ddc[k+1]*b*b*b*sa0);
            break;
        case 3: k = 5;
SINGLE:     p4[0] = (float)(24.0*ddc[k]*a*a*a*a);
            p4[1] = (float)(24.0*ddc[k]*a*a*a*b);
            p4[2] = (float)(24.0*ddc[k]*a*a*b*b);
            p4[3] = (float)(24.0*ddc[k]*a*b*b*b);
            p4[4] = (float)(24.0*ddc[k]*b*b*b*b);
            break;
          } /* switch */
          for ( i = 0; i < hole_k; i++ ) {
            sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
            n = sa*cb-ca*sb;
            if ( n > TOL ) {
              m = ca0*cb+sa0*sb;
              if ( m > TOL ) {
                if ( !spartition[sfn+3-j].sgn || spartition[sfn+3-j].both )
                  ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                                      br0bc, br0cr1bc, br0cr2bc, 1 );
              }
              else if ( m < -TOL ) {
                if ( spartition[sfn+3-j].sgn && !spartition[sfn+3-j].both )
                  ReparHomoPoly4Der ( p4, i, privateG2->AuxiMat,
                                      br0bc, br0cr1bc, br0cr2bc, 2 );
              }
            }
          }
        }
      }
    }
    break;
  }

    /* find basis functions "auxiliary patches" */
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG,
                    5, (G2H_OMCDEG+1), br0bc, 3, (G2H_OMCDEG+1),
                    &br0bc[5], G2_CROSS00DEG+1, br0 ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG-1,
                    4, G2H_OMCDEG, br0cr1bc, 3, G2H_OMCDEG,
                    &br0cr1bc[4], G2_CROSS00DEG, br0cr1 ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG-2,
                    3, (G2H_OMCDEG-1), br0cr2bc, 3, (G2H_OMCDEG-1),
                    &br0cr2bc[3], G2_CROSS00DEG-1, br0cr2 ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef TOL
} /*_g2h_GetABasisAuxpf*/

static boolean GetABasisFunctionf ( GHoleDomainf *domain, int fn,
                                    float *bbr0,
                                    float *bbr0cr1, float *bbq0cr1,
                                    float *bbr0cr2, float *bbq0cr2 )
{
#define TOL 1.0e-6
  void     *sp;
  G2HolePrivateRecf *privateG2;
  int      hole_k;
  float    *br0cr1, *br0cr2;
  float    *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
           *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
           *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01;
  int      i, j;


  sp = pkv_GetScratchMemTop ();

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;

  br0cr1 = pkv_GetScratchMemf ( (2*G2H_OMCDEG-1)*hole_k );
  if ( !br0cr1 )
    goto failure;
  br0cr2 = &br0cr1[hole_k*G2H_OMCDEG];

  if ( !_g2h_GetABasisAuxpf ( domain, fn, bbr0, br0cr1, br0cr2 ) )
    goto failure;

    /* now find the proper cross derivatives of the basis function */
  G2GetPolynomialAddresses0 ( privateG2->jfunc,
      b01, c01, f01, g01, b11, c11, f11, g11,
      b02, c02, f02, g02, b12, c12, f12, g12, b01b01, twob01c01, c01c01,
      f01f01, twof01g01, g01g01 );

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    FindDiCrossDeraf ( 1,
        &b01[i*(G2_BF01DEG+1)], &c01[i*(G2_CG01DEG+1)],
        &b02[i*(G2_BF02DEG+1)], &c02[i*(G2_CG02DEG+1)],
        &b01b01[i*(2*G2_BF01DEG+1)], &twob01c01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &c01c01[i*(2*G2_CG01DEG+1)],
        &bbr0[i*(G2H_OMCDEG+1)], &br0cr1[i*G2H_OMCDEG], &br0cr2[i*(G2H_OMCDEG-1)],
        &bbr0cr1[i*(G2_CROSS01DEG+1)], &bbr0cr2[i*(G2_CROSS02DEG+1)] );
    FindDiCrossDeraf ( 1,
        &f01[i*(G2_BF01DEG+1)], &g01[i*(G2_CG01DEG+1)],
        &f02[i*(G2_BF02DEG+1)], &g02[i*(G2_CG02DEG+1)],
        &f01f01[i*(2*G2_BF01DEG+1)], &twof01g01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &g01g01[i*(2*G2_CG01DEG+1)],
        &bbr0[j*(G2H_OMCDEG+1)], &br0cr1[j*G2H_OMCDEG], &br0cr2[j*(G2H_OMCDEG-1)],
        &bbq0cr1[i*(G2_CROSS01DEG+1)], &bbq0cr2[i*(G2_CROSS02DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef TOL
} /*GetABasisFunctionf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2h_GetBBasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *bezfc, float *fcomc,
                              float *fcomcd, float *fcomcdd )
{
                /* this procedure is similar to FindAuxDPatchesf */
  void     *sp;
  int      hole_k;
  float    *bez, *bezbc;
  float    *fcomcbc, *fcomcbcd, *fcomcbcdd;
  int      i, j;
  float    *hole_knots;
  float    g11, g11p2, h11, h11p2;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  hole_knots = domain->hole_knots;

  fcomcbc = pkv_GetScratchMemf ( 3*G2H_OMCDEG*hole_k );
  bez     = pkv_GetScratchMemf ( 24*hole_k );
  bezbc   = pkv_GetScratchMemf ( 36*hole_k );
  if ( !fcomcbc || !bez || !bezbc )
    goto failure;
  fcomcbcd = &fcomcbc[(G2H_OMCDEG+1)*hole_k];
  fcomcbcdd = &fcomcbcd[G2H_OMCDEG*hole_k];

        /* initialise the boundary conditions to zero */
  memset ( fcomcbc, 0, 3*G2H_OMCDEG*hole_k*sizeof(float) );

        /* setup the nonzero boundary conditions */
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      gh_GetDomSurrndBFuncf ( domain, fn, i, j+1, &bezfc[(2*i+j)*16] );
      if ( !mbs_multiBCHornerDer2f ( 3, 1, 4, 0, &bezfc[(2*i+j)*16], 0.0,
              &bez[(2*i+j)*12], &bez[(2*i+j)*12+4], &bez[(2*i+j)*12]+8 ) )
        goto failure;
      memcpy ( &bezfc[(2*i+j)*16+4], &bez[(2*i+j)*12+4], 8*sizeof(float) );
      if ( !mbs_multiBCHornerDer2f ( 3, 3, 1, 4, &bez[(2*i+j)*12], 0.0,
              &bezbc[(2*i+j)*18], &bezbc[(2*i+j)*18+3], &bezbc[(2*i+j)*18+6] ) )
        goto failure;
      if ( !mbs_multiBCHornerDer2f ( 3, 3, 1, 4, &bez[(2*i+j)*12], 1.0,
              &bezbc[(2*i+j)*18+9], &bezbc[(2*i+j)*18+12], &bezbc[(2*i+j)*18+15] ) )
        goto failure;
    }
  }

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    g11 = (hole_knots[11*j+5]-hole_knots[11*j+4])/
          (hole_knots[11*j+4]-hole_knots[11*j+3]);
    g11p2 = g11*g11;
    h11 = (float)((hole_knots[11*i+6]-hole_knots[11*i+4])/
                  (2.0*(hole_knots[11*i+6]-hole_knots[11*i+5])));
    h11p2 = h11*h11;
    fcomcbc[i*(G2H_OMCDEG+1)+5] = bezbc[36*i];
    fcomcbc[i*(G2H_OMCDEG+1)+6] = g11*bezbc[(2*i+1)*18+1];
    fcomcbc[i*(G2H_OMCDEG+1)+7] = g11p2*bezbc[(2*i+1)*18+2];
    fcomcbcd[i*G2H_OMCDEG+4] = h11*bezbc[(2*i+1)*18+3];
    fcomcbcd[i*G2H_OMCDEG+5] = g11*h11*bezbc[(2*i+1)*18+4];
    fcomcbcd[i*G2H_OMCDEG+6] = g11p2*h11*bezbc[(2*i+1)*18+5];
    fcomcbcdd[i*(G2H_OMCDEG-1)+3] = h11p2*bezbc[(2*i+1)*18+6];
    fcomcbcdd[i*(G2H_OMCDEG-1)+4] = g11*h11p2*bezbc[(2*i+1)*18+7];
    fcomcbcdd[i*(G2H_OMCDEG-1)+5] = g11p2*h11p2*bezbc[(2*i+1)*18+8];
  }
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG,
           5, (G2H_OMCDEG+1), (float*)fcomcbc, 3, (G2H_OMCDEG+1),
           (float*)&fcomcbc[5], (G2H_OMCDEG+1), (float*)fcomc ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG-1,
           4, G2H_OMCDEG, (float*)fcomcbcd, 3, G2H_OMCDEG,
           (float*)&fcomcbcd[4], G2H_OMCDEG, (float*)fcomcd ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G2_CROSS00DEG-2,
           3, (G2H_OMCDEG-1), (float*)fcomcbcdd, 3, (G2H_OMCDEG-1),
           (float*)&fcomcbcdd[3], (G2H_OMCDEG-1), (float*)fcomcdd ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_GetBBasisAuxpf*/

static boolean GetBBasisFunctionf ( GHoleDomainf *domain, int fn,
                                    float *bbr0, float *bbr1, float *bbq1,
                                    float *bbr0cr1, float *bbq0cr1,
                                    float *bbr1cr1, float *bbq1cr1,
                                    float *bbr0cr2, float *bbq0cr2,
                                    float *bbr1cr2, float *bbq1cr2 )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  int   hole_k;
  int   i, j;
  float *bezfc, *fcomc, *fcomcd, *fcomcdd;
  float *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11,
        *b02, *c02, *f02, *g02, *b12, *c12, *f12, *g12,
        *b01b01, *twob01c01, *c01c01, *f01f01, *twof01g01, *g01g01,
        *b11b11, *twob11c11, *c11c11, *f11f11, *twof11g11, *g11g11;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  bezfc = pkv_GetScratchMemf ( 32*hole_k );
  fcomc = pkv_GetScratchMemf ( 3*G2H_OMCDEG*hole_k );
  if ( !bezfc || !fcomc )
    goto failure;
  fcomcd = &fcomc[(G2H_OMCDEG+1)*hole_k];
  fcomcdd = &fcomcd[G2H_OMCDEG*hole_k];

  if ( !_g2h_GetBBasisAuxpf ( domain, fn, bezfc, fcomc, fcomcd, fcomcdd ) )
    goto failure;

  G2GetPolynomialAddresses ( privateG2->jfunc,
      b01, c01, f01, g01, b11, c11, f11, g11, 
      b02, c02, f02, g02, b12, c12, f12, g12, b01b01, twob01c01, c01c01,
      f01f01, twof01g01, g01g01, b11b11, twob11c11, c11c11,
      f11f11, twof11g11, g11g11 );

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    memcpy ( &bbr0[i*(G2_CROSS00DEG+1)], &fcomc[i*(G2_CROSS00DEG+1)],
             (G2_CROSS00DEG+1)*sizeof(float) );
    memcpy ( &bbr1[i*(G2_CROSS10DEG+1)], &bezfc[32*j],
             (G2_CROSS10DEG+1)*sizeof(float) );
    memcpy ( &bbq1[i*(G2_CROSS10DEG+1)], &bezfc[(2*i+1)*16],
             (G2_CROSS10DEG+1)*sizeof(float) );
    FindDiCrossDeraf ( 1,
        &b01[i*(G2_BF01DEG+1)], &c01[i*(G2_CG01DEG+1)],
        &b02[i*(G2_BF02DEG+1)], &c02[i*(G2_CG02DEG+1)],
        &b01b01[i*(2*G2_BF01DEG+1)], &twob01c01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &c01c01[i*(2*G2_CG01DEG+1)],
        &fcomc[i*(G2H_OMCDEG+1)], &fcomcd[i*G2H_OMCDEG], &fcomcdd[i*(G2H_OMCDEG-1)],
        &bbr0cr1[i*(G2_CROSS01DEG+1)], &bbr0cr2[i*(G2_CROSS02DEG+1)] );
    FindDiCrossDeraf ( 1,
        &f01[i*(G2_BF01DEG+1)], &g01[i*(G2_CG01DEG+1)],
        &f02[i*(G2_BF02DEG+1)], &g02[i*(G2_CG02DEG+1)],
        &f01f01[i*(2*G2_BF01DEG+1)], &twof01g01[i*(G2_BF01DEG+G2_CG01DEG+1)],
        &g01g01[i*(2*G2_CG01DEG+1)],
        &fcomc[j*(G2H_OMCDEG+1)], &fcomcd[j*G2H_OMCDEG], &fcomcdd[j*(G2H_OMCDEG-1)],
        &bbq0cr1[i*(G2_CROSS01DEG+1)], &bbq0cr2[i*(G2_CROSS02DEG+1)] );
    FindDiCrossDerbf ( 1,
        &b11[i*(G2_BF11DEG+1)], &c11[i*(G2_CG11DEG+1)],
        &b12[i*(G2_BF12DEG+1)], &c12[i*(G2_CG12DEG+1)],
        &b11b11[i*(2*G2_BF11DEG+1)], &twob11c11[i*(G2_BF11DEG+G2_CG11DEG+1)],
        &c11c11[i*(2*G2_CG11DEG+1)],
        &bezfc[32*j], &bezfc[32*j+4], &bezfc[32*j+8],
        &bbr1cr1[i*(G2_CROSS11DEG+1)], &bbr1cr2[i*(G2_CROSS12DEG+1)] );
    FindDiCrossDerbf ( 1,
        &f11[i*(G2_BF11DEG+1)], &g11[i*(G2_CG11DEG+1)],
        &f12[i*(G2_BF12DEG+1)], &g12[i*(G2_CG12DEG+1)],
        &f11f11[i*(2*G2_BF11DEG+1)], &twof11g11[i*(G2_BF11DEG+G2_CG11DEG+1)],
        &g11g11[i*(2*G2_CG11DEG+1)],
        &bezfc[16*(2*i+1)], &bezfc[16*(2*i+1)+4], &bezfc[16*(2*i+1)+8],
        &bbq1cr1[i*(G2_CROSS11DEG+1)], &bbq1cr2[i*(G2_CROSS12DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GetBBasisFunctionf*/

static boolean FindBasisFunctionsf ( GHoleDomainf *domain )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  int   hole_k;
  int   i, k, n, nfunc_a, nfunc_b;
  float *bbr0,              *bbr1,    *bbq1,
        *bbr0cr1, *bbq0cr1, *bbr1cr1, *bbq1cr1,
        *bbr0cr2, *bbq0cr2, *bbr1cr2, *bbq1cr2;
  int   option, ndata, *idata;
  float *fdata;

  sp = pkv_GetScratchMemTop ();

  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;

    /* find the Coons representation of the "inner" basis functions;  */
    /* they are represented only by the values at the "inner" curves, */ 
    /* as they vanish together with their derivatives at the boundary */
    /* of the domain */
  if ( !AnalyzePartitionf ( domain ) )
    goto failure;

  option = privateG2->GetOption ( domain, G2HQUERY_BASIS, 0, &ndata, &idata, &fdata );
  privateG2->opt_basis = (char)option;
  switch ( option ) {
case G2H_DEFAULT:
    break;

case G2H_USE_RESTRICTED_BASIS:
    privateG2->nfunc_a = 15;  /* use the subspace spanned by the functions */
    break;                  /* determined by polynomials only */

default:
    domain->error_code = G2H_ERROR_INVALID_OPTION;
    goto failure;
  }
  privateG2->opt_quad = (char)privateG2->GetOption ( domain, G2HQUERY_QUADRATURE, 0,
                                               &ndata, &idata, &fdata );
  switch ( privateG2->opt_quad ) {
case G2H_DEFAULT:
case G2H_QUADRATURE_GAUSS_LEGENDRE:
    break;

default:
    domain->error_code = G2H_ERROR_INVALID_OPTION;
    goto failure;
  }

  nfunc_a = privateG2->nfunc_a;
  n = G2_CROSS00DEG+1+2*(G2_CROSS01DEG+G2_CROSS02DEG+2);
  privateG2->basis_a = malloc ( nfunc_a*hole_k*n*sizeof(float) );
  privateG2->AuxiMat = malloc ( hole_k*125*sizeof(float) );
  if ( !privateG2->basis_a || !privateG2->AuxiMat ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    goto failure;
  }

  ComputeAuxiMat ( domain );

  G2GetBFuncACrossAddresses ();
  for ( i = k = 0;  i < nfunc_a;  i++, k += hole_k )
    if ( !GetABasisFunctionf ( domain, i,
                  &bbr0[k*(G2_CROSS00DEG+1)],
                  &bbr0cr1[k*(G2_CROSS01DEG+1)], &bbq0cr1[k*(G2_CROSS01DEG+1)],
                  &bbr0cr2[k*(G2_CROSS02DEG+1)], &bbq0cr2[k*(G2_CROSS02DEG+1)] ) )
      goto failure;

    /* find the Coons representation of the "outer" basis functions;  */
  privateG2->nfunc_b = nfunc_b = 6*hole_k+1;
  n = G2_CROSS00DEG+1+2*(G2_CROSS01DEG+G2_CROSS02DEG+G2_CROSS10DEG+G2_CROSS11DEG+G2_CROSS12DEG+5);
  privateG2->basis_b = malloc ( i = nfunc_b*hole_k*n*sizeof(float) );
  if ( !privateG2->basis_b ) {
    domain->error_code = G2H_ERROR_CANNOT_MALLOC;
    goto failure;
  }
  G2GetBFuncBCrossAddresses ();
  memset ( bbr0, 0, i );

  for ( i = k = 0;  i < nfunc_b;  i++, k += hole_k )
    if ( !GetBBasisFunctionf ( domain, i, &bbr0[k*(G2_CROSS00DEG+1)],
                  &bbr1[k*(G2_CROSS10DEG+1)], &bbq1[k*(G2_CROSS10DEG+1)],
                  &bbr0cr1[k*(G2_CROSS01DEG+1)], &bbq0cr1[k*(G2_CROSS01DEG+1)],
                  &bbr1cr1[k*(G2_CROSS11DEG+1)], &bbq1cr1[k*(G2_CROSS11DEG+1)],
                  &bbr0cr2[k*(G2_CROSS02DEG+1)], &bbq0cr2[k*(G2_CROSS02DEG+1)],
                  &bbr1cr2[k*(G2_CROSS12DEG+1)], &bbq1cr2[k*(G2_CROSS12DEG+1)] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindBasisFunctionsf*/

