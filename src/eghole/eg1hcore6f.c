
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
  G1HolePrivateRecf *privateG1;
  int            auxp_k;
  unsigned char  *spdimen;

  privateG1 = domain->privateG1;
  spdimen = &privateG1->spdimen[0];

  if ( _gh_AnalyzePartitionf ( domain, G1H_OMCDEG, privateG1->omcbc,
                      &privateG1->hole_m, &auxp_k, &privateG1->partition,
                      &privateG1->spartition, &privateG1->spart_alpha0 ) ) {
    spdimen[0] = (unsigned char)privateG1->hole_m;    /* quadratic half-polynomials */
    spdimen[1] = (unsigned char)max ( 0, auxp_k-3 );  /* quadratic B-splines */
    privateG1->nfunc_a = 6 + privateG1->hole_m + spdimen[1];
/*
printf ( "spdimen: %d %d %d\n", spdimen[0], spdimen[1], private->nfunc_a );
*/
    return true;
  }
  else
    return false;
} /*AnalyzePartitionf*/

static void ComputeAuxiMat ( GHoleDomainf *domain )
{
  G1HolePrivateRecf *privateG1;
  int      i, hole_k;
  vector2f diu, div, diuu, diuv;
  vector2f *omcbc, *omcbcd;
  float    *A11, *A21, *A22;

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  omcbc = privateG1->omcbc;
  omcbcd = &omcbc[(G1H_OMCDEG+1)*hole_k];

  for ( i = 0; i < hole_k; i++ ) {

        /* setup the addresses of the transformation matrices */
    A11 = &privateG1->AuxiMat[i*19];
    A21 = &A11[4];   A22 = &A21[6];

        /* extract the derivatives of the domain patches */
    diu = omcbc[(G1H_OMCDEG+1)*i+1];     div = omcbcd[G1H_OMCDEG*i];
    diuu = omcbc[(G1H_OMCDEG+1)*i+2];    diuv = omcbcd[G1H_OMCDEG*i+1];

        /* compute the derivative transformation matrices */
    pkn_Setup2DerA11Matrixf ( diu.x, diu.y, div.x, div.y, A11 );
    pkn_Setup2DerA21Matrixf ( diuu.x, diuu.y, diuv.x,
                              diuv.y, 0.0, 0.0, A21 );
    pkn_Setup2DerA22Matrixf ( diu.x, diu.y, div.x, div.y, A22 );
  }
} /*ComputeAuxiMat*/

/* ///////////////////////////////////////////////////////////////////////// */
static void FindDivDiff3Coeff ( float a, float b, float c, float d, float *dd )
{
    /* compute the coefficients of the divided difference of 3rd order */
    /* for single knots a,b,c,d, i.e. the numbers such that            */
    /* f[a,b,c,d] = dd[0]f(a)+dd[1]f(b)+dd[2]f(c)+dd[3]f(d)            */
  dd[0] = (float)(1.0/((a-b)*(a-c)*(a-d)));
  dd[1] = (float)(1.0/((b-a)*(b-c)*(b-d)));
  dd[2] = (float)(1.0/((c-a)*(c-b)*(c-d)));
  dd[3] = (float)(1.0/((d-a)*(d-b)*(d-c)));
} /*FindDivDiff4Coeff*/

static void ReparHomoPoly1Der ( float *p1, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc )
{
  float q1[2], q2[3];
  float *A11, *A21;

  A11 = &AuxiMat[i*19];  A21 = &A11[4];
  pkn_MultMatrixf ( 2, 2, 2, A11, 1, 1, p1, 1, q1 );
  pkn_MultMatrixf ( 3, 2, 2, A21, 1, 1, p1, 1, q2 );
  br0bc[(G1H_OMCDEG+1)*i+1]    = q1[0];
  br0bc[(G1H_OMCDEG+1)*i+2]    = q2[0];
  br0cr1bc[G1H_OMCDEG*i]   = q1[1];
  br0cr1bc[G1H_OMCDEG*i+1] = q2[1];
} /*ReparHomoPoly1Der*/

static void ReparHomoPoly2Der ( float *p2, int i, float *AuxiMat,
                                float *br0bc, float *br0cr1bc,
                                char sw )
{
  float q2[4];
  float *A22;

  A22 = &AuxiMat[i*19+4+6];
  pkn_MultMatrixf ( 3, 3, 3, A22, 1, 1, p2, 1, q2 );
  switch ( sw ) {
case 0:
    br0bc[(G1H_OMCDEG+1)*i+2]    = q2[0];
    br0cr1bc[G1H_OMCDEG*i+1] = q2[1];
    break;
case 1:
    br0bc[(G1H_OMCDEG+1)*i+2]    += q2[0];
    br0cr1bc[G1H_OMCDEG*i+1] += q2[1];
    break;
case 2:
    br0bc[(G1H_OMCDEG+1)*i+2]    -= q2[0];
    br0cr1bc[G1H_OMCDEG*i+1] -= q2[1];
    break;
  }
} /*ReparHomoPoly2Der*/

boolean _g1h_GetABasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *br0, float *br0cr1 )
{
#define TOL 1.0e-6
  void     *sp;
  G1HolePrivateRecf *privateG1;
  int      hole_k;
  unsigned char *spdimen;
  float    *partition;
  GHoleSgnPartf *spartition;
  float    *br0bc, *br0cr1bc;
  float    p1[2], p2[3];
  float    alpha0, sa0, ca0, sa, sb, ca, cb, a, b;
  float    u, n, m, normf;
  float    ddc[6];
  int      sfn;
  int      i, j, k;


  sp = pkv_GetScratchMemTop ();

  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  spdimen = privateG1->spdimen;

  br0bc = pkv_GetScratchMemf ( (2*G1H_OMCDEG+1)*hole_k );
  if ( !br0bc )
    goto failure;

  br0cr1bc = &br0bc[hole_k*(G1H_OMCDEG+1)];

  memset ( br0bc, 0, hole_k*3*G1_CROSS00DEG*sizeof(float) );
  switch ( fn ) {
case 0:    /* f(x,y) = 1 */
    for ( i = 0; i < hole_k; i++ )
      br0bc[(G1H_OMCDEG+1)*i] = 1.0;
    break;

case 1:    /* f(x,y) = x */
    p1[0] = 1.0;  p1[1] = 0.0;
    goto continue_deg1;

case 2:    /* f(x,y) = y */
    p1[0] = 0.0;  p1[1] = 1.0;
continue_deg1:
    for ( i = 0; i < hole_k;  i++ )
      ReparHomoPoly1Der ( p1, i, privateG1->AuxiMat, br0bc, br0cr1bc );
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
      ReparHomoPoly2Der ( p2, i, privateG1->AuxiMat, br0bc, br0cr1bc, 0 );
    break;

default:   /* splines */
    partition = privateG1->partition;
    spartition = privateG1->spartition;
    alpha0 = privateG1->spart_alpha0;
    sfn = fn-6;        /* which one after the polynomials? */

    if ( sfn < spdimen[0] ) {                 /* a quadratic half-polynomial */
      for ( j = k = -1;  k < sfn;  j++ )
        if ( spartition[j+1].both ) k++;
      sa = (float)sin ( spartition[j].malpha );
      ca = (float)cos ( spartition[j].malpha );
      p2[0] = (float)(-2.0*sa*sa);  p2[1] = (float)(2.0*sa*ca);
      p2[2] = (float)(-2.0*ca*ca);
      for ( i = 0; i < hole_k; i++ ) {
        sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
        if ( sa*cb-ca*sb > TOL )
          ReparHomoPoly2Der ( p2, i, privateG1->AuxiMat, br0bc, br0cr1bc, 0 );
      }
    }
    else if ( sfn-spdimen[0] < spdimen[1] ) { /* a quadratic B-spline */
      sfn -= spdimen[0];
      sa0 = (float)sin ( alpha0 );  ca0 = (float)cos ( alpha0 );
      FindDivDiff3Coeff ( spartition[sfn+3].knot, spartition[sfn+2].knot,
                          spartition[sfn+1].knot, spartition[sfn].knot, ddc );
      normf = spartition[sfn].knot-spartition[sfn+3].knot;
      for ( j = 0; j < 4; j++ )
        ddc[j] *= normf;

      for ( j = 0; j < 4; j++ ) {
        sa = (float)sin ( spartition[sfn+3-j].malpha );
        ca = (float)cos ( spartition[sfn+3-j].malpha );
        u = spartition[sfn+3-j].knot;  a = sa0-u*ca0;  b = -(ca0+u*sa0);
            /* compute the partial derivatives of the third order */
            /* of the polynomial (X-uY)^2 = (ax+by)^2, times ddc[j] */
        p2[0] = (float)(2.0*ddc[j]*a*a);
        p2[1] = (float)(2.0*ddc[j]*a*b);
        p2[2] = (float)(2.0*ddc[j]*b*b);

        for ( i = 0; i < hole_k; i++ ) {
          sb = (float)sin ( partition[i] );  cb = (float)cos ( partition[i] );
          n = sa*cb-ca*sb;
          if ( n > TOL ) {
            m = ca0*cb+sa0*sb;
            if ( m > TOL ) {
              if ( !spartition[sfn+3-j].sgn || spartition[sfn+3-j].both )
                ReparHomoPoly2Der ( p2, i, privateG1->AuxiMat, br0bc, br0cr1bc, 1 );
            }
            else if ( m < -TOL ) {
              if ( spartition[sfn+3-j].sgn && !spartition[sfn+3-j].both )
                ReparHomoPoly2Der ( p2, i, privateG1->AuxiMat, br0bc, br0cr1bc, 2 );
            }
          }
        }
      }
    }
    break;
  }

    /* find basis functions "auxiliary patches" */
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G1_CROSS00DEG,
                    3, 5, br0bc, 2, 5, &br0bc[3], G1_CROSS00DEG+1, br0 ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G1_CROSS00DEG-1,
                    2, 4, br0cr1bc, 2, 4, &br0cr1bc[2], G1_CROSS00DEG, br0cr1 ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef TOL
} /*_g1h_GetABasisAuxpf*/

static boolean _g1h_GetABasisFunctionf ( GHoleDomainf *domain, int fn,
                                         float *bbr0,
                                         float *bbr0cr1, float *bbq0cr1 )
{
#define TOL 1.0e-6
  void     *sp;
  G1HolePrivateRecf *privateG1;
  int      hole_k;
  float    *br0cr1;
  float    *b01, *c01, *f01, *g01;
  int      i, j;


  sp = pkv_GetScratchMemTop ();

  privateG1 = domain->privateG1;
  hole_k    = domain->hole_k;

  br0cr1 = pkv_GetScratchMemf ( hole_k*G1H_OMCDEG );
  if ( !br0cr1 )
    goto failure;

  if ( !_g1h_GetABasisAuxpf ( domain, fn, bbr0, br0cr1 ) )
    goto failure;

    /* now find the proper cross derivatives of the basis function */
  G1GetPolyAddr0 ( privateG1->jfunc, b01, c01, f01, g01 );

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    FindDiCrossDeraf ( 1, &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
        &bbr0[i*(G1H_OMCDEG+1)], &br0cr1[i*G1H_OMCDEG], &bbr0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDeraf ( 1, &f01[i*(G1_BF01DEG+1)], &g01[i*(G1_CG01DEG+1)],
        &bbr0[j*(G1H_OMCDEG+1)], &br0cr1[j*G1H_OMCDEG], &bbq0cr1[i*(G1_CROSS01DEG+1)] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef TOL
} /*_g1h_GetABasisFunctionf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_GetBBasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *bezfc, float *fcomc, float *fcomcd )
{
                /* this procedure is similar to FindAuxDPatchesf */
  void     *sp;
  int      hole_k;
  float    *bez, *bezbc;
  float    *fcomcbc, *fcomcbcd;
  int      i, j;
  float    *hole_knots;
  float    g11, h11;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  hole_knots = domain->hole_knots;

  fcomcbc = pkv_GetScratchMemf ( (2*G1H_OMCDEG+1)*hole_k );
  bez     = pkv_GetScratchMemf ( 24*hole_k );
  bezbc   = pkv_GetScratchMemf ( 36*hole_k );
  if ( !fcomcbc || !bez || !bezbc )
    goto failure;
  fcomcbcd = &fcomcbc[(G1H_OMCDEG+1)*hole_k];

        /* initialise the boundary conditions to zero */
  memset ( fcomcbc, 0, (2*G1H_OMCDEG+1)*hole_k*sizeof(float) );

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
    h11 = (float)((hole_knots[11*i+6]-hole_knots[11*i+4])/
                  (2.0*(hole_knots[11*i+6]-hole_knots[11*i+5])));
    fcomcbc[i*(G1H_OMCDEG+1)+3] = bezbc[36*i];
    fcomcbc[i*(G1H_OMCDEG+1)+4] = g11*bezbc[(2*i+1)*18+1];
    fcomcbcd[i*G1H_OMCDEG+2] = h11*bezbc[(2*i+1)*18+3];
    fcomcbcd[i*G1H_OMCDEG+3] = g11*h11*bezbc[(2*i+1)*18+4];
  }
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G1_CROSS00DEG,
           3, (G1H_OMCDEG+1), (float*)fcomcbc, 2, (G1H_OMCDEG+1),
           (float*)&fcomcbc[3], (G1H_OMCDEG+1), (float*)fcomc ) )
    goto failure;
  if ( !mbs_multiInterp2knHermiteBezf ( hole_k, 1, G1_CROSS00DEG-1,
           2, G1H_OMCDEG, (float*)fcomcbcd, 2, G1H_OMCDEG,
           (float*)&fcomcbcd[2], G1H_OMCDEG, (float*)fcomcd ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_GetBBasisAuxpf*/

static boolean GetBBasisFunctionf ( GHoleDomainf *domain, int fn,
                                    float *bbr0, float *bbr1, float *bbq1,
                                    float *bbr0cr1, float *bbq0cr1,
                                    float *bbr1cr1, float *bbq1cr1 )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   hole_k;
  int   i, j;
  float *bezfc, *fcomc, *fcomcd;
  float *b01, *c01, *f01, *g01, *b11, *c11, *f11, *g11;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  bezfc = pkv_GetScratchMemf ( 32*hole_k );
  fcomc = pkv_GetScratchMemf ( (2*G1H_OMCDEG+1)*hole_k );
  if ( !bezfc || !fcomc )
    goto failure;
  fcomcd = &fcomc[(G1H_OMCDEG+1)*hole_k];

  if ( !_g1h_GetBBasisAuxpf ( domain, fn, bezfc, fcomc, fcomcd ) )
    goto failure;

  G1GetPolyAddr ( privateG1->jfunc, b01, c01, f01, g01, b11, c11, f11, g11 );

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    memcpy ( &bbr0[i*(G1_CROSS00DEG+1)], &fcomc[i*(G1_CROSS00DEG+1)],
             (G1_CROSS00DEG+1)*sizeof(float) );
    memcpy ( &bbr1[i*(G1_CROSS10DEG+1)], &bezfc[32*j],
             (G1_CROSS10DEG+1)*sizeof(float) );
    memcpy ( &bbq1[i*(G1_CROSS10DEG+1)], &bezfc[(2*i+1)*16],
             (G1_CROSS10DEG+1)*sizeof(float) );
    FindDiCrossDeraf ( 1, &b01[i*(G1_BF01DEG+1)], &c01[i*(G1_CG01DEG+1)],
        &fcomc[i*(G1H_OMCDEG+1)], &fcomcd[i*G1H_OMCDEG], &bbr0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDeraf ( 1, &f01[i*(G1_BF01DEG+1)], &g01[i*(G1_CG01DEG+1)],
        &fcomc[j*(G1H_OMCDEG+1)], &fcomcd[j*G1H_OMCDEG], &bbq0cr1[i*(G1_CROSS01DEG+1)] );
    FindDiCrossDerbf ( 1, &b11[i*(G1_BF11DEG+1)], &c11[i*(G1_CG11DEG+1)],
        &bezfc[32*j], &bezfc[32*j+4], &bbr1cr1[i*(G1_CROSS11DEG+1)] );
    FindDiCrossDerbf ( 1, &f11[i*(G1_BF11DEG+1)], &g11[i*(G1_CG11DEG+1)],
        &bezfc[16*(2*i+1)], &bezfc[16*(2*i+1)+4], &bbq1cr1[i*(G1_CROSS11DEG+1)] );
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
  G1HolePrivateRecf *privateG1;
  int   hole_k;
  int   i, k, n, nfunc_a, nfunc_b;
  float *bbr0,              *bbr1,    *bbq1,
        *bbr0cr1, *bbq0cr1, *bbr1cr1, *bbq1cr1;
  int   option, ndata, *idata;
  float *fdata;

  sp = pkv_GetScratchMemTop ();

  hole_k    = domain->hole_k;
  privateG1 = domain->privateG1;

    /* find the Coons representation of the "inner" basis functions;  */
    /* they are represented only by the values at the "inner" curves, */ 
    /* as they vanish together with their derivatives at the boundary */
    /* of the domain */
  if ( !AnalyzePartitionf ( domain ) )
    goto failure;

  option = privateG1->GetOption ( domain, G1HQUERY_BASIS, 0, &ndata, &idata, &fdata );
  privateG1->opt_basis = (char)option;
  switch ( option ) {
case G1H_DEFAULT:
    break;

case G1H_USE_RESTRICTED_BASIS:
    privateG1->nfunc_a = 6;   /* use the subspace spanned by the functions */
    break;                  /* determined by polynomials only */

default:
    domain->error_code = G1H_ERROR_INVALID_OPTION;
    goto failure;
  }
  privateG1->opt_quad = (char)privateG1->GetOption ( domain, G1HQUERY_QUADRATURE, 0,
                                                     &ndata, &idata, &fdata );
  switch ( privateG1->opt_quad ) {
case G1H_DEFAULT:
case G1H_QUADRATURE_GAUSS_LEGENDRE:
    break;

default:
    domain->error_code = G1H_ERROR_INVALID_OPTION;
    goto failure;
  }

  nfunc_a = privateG1->nfunc_a;
  n = G1_CROSS00DEG+1+2*(G1_CROSS01DEG+1);
  privateG1->basis_a = malloc ( nfunc_a*hole_k*n*sizeof(float) );
  privateG1->AuxiMat = malloc ( hole_k*19*sizeof(float) );
  if ( !privateG1->basis_a || !privateG1->AuxiMat ) {
    domain->error_code = G1H_ERROR_CANNOT_MALLOC;
    goto failure;
  }

  ComputeAuxiMat ( domain );

  G1GetBFuncACrossAddresses ();
  for ( i = k = 0;  i < nfunc_a;  i++, k += hole_k )
    if ( !_g1h_GetABasisFunctionf ( domain, i,
                  &bbr0[k*(G1_CROSS00DEG+1)],
                  &bbr0cr1[k*(G1_CROSS01DEG+1)], &bbq0cr1[k*(G1_CROSS01DEG+1)] ) )
      goto failure;

    /* find the Coons representation of the "outer" basis functions;  */
  privateG1->nfunc_b = nfunc_b = 6*hole_k+1;
  n = G1_CROSS00DEG+1+2*(G1_CROSS01DEG+G1_CROSS10DEG+G1_CROSS11DEG+3);
  privateG1->basis_b = malloc ( i = nfunc_b*hole_k*n*sizeof(float) );
  if ( !privateG1->basis_b ) {
    domain->error_code = G1H_ERROR_CANNOT_MALLOC;
    goto failure;
  }
  G1GetBFuncBCrossAddresses ();
  memset ( bbr0, 0, i );

  for ( i = k = 0;  i < nfunc_b;  i++, k += hole_k )
    if ( !GetBBasisFunctionf ( domain, i,
                  &bbr0[k*(G1_CROSS00DEG+1)],
                  &bbr1[k*(G1_CROSS10DEG+1)], &bbq1[k*(G1_CROSS10DEG+1)],
                  &bbr0cr1[k*(G1_CROSS01DEG+1)], &bbq0cr1[k*(G1_CROSS01DEG+1)],
                  &bbr1cr1[k*(G1_CROSS11DEG+1)], &bbq1cr1[k*(G1_CROSS11DEG+1)] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindBasisFunctionsf*/

