
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void _g1h_DiJacobian2f ( const vector2f *du, const vector2f *dv,
                         const vector2f *duu, const vector2f *duv,
                         const vector2f *dvv,
                         float *jac, float *trd )
{
  vector2f gx, gy, gxx, gxy, gyy;
  float    A21[6], A22[9];

  *jac = (float)fabs ( du->x*dv->y - du->y*dv->x );

  pkn_f2iDerivatives2f ( du->x, du->y, dv->x, dv->y,
      duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
      (float*)&gx, (float*)&gy, (float*)&gxx, (float*)&gxy, (float*)&gyy );

  pkn_Setup2DerA21Matrixf ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
  pkn_Setup2DerA22Matrixf ( gx.x, gx.y, gy.x, gy.y, A22 );

  trd[0]  = A21[0]+A21[4];   trd[1]  = A21[1]+A21[5];
  trd[2]  = A22[0]+A22[6];   trd[3]  = A22[1]+A22[7];   trd[4] = A22[2]+A22[8];
} /*_g1h_DiJacobian2f*/

boolean _g1h_TabDiPatchJac2f ( int nkn, const float *kn, const float *hfunc,
                               const float *dhfunc, const float *ddhfunc,
                               const vector2f *c00, const vector2f *c01,
                               const vector2f *c10, const vector2f *c11,
                               const vector2f *d00, const vector2f *d01,
                               const vector2f *d10, const vector2f *d11,
                               float *jac, float *trd )
{
  void     *sp;
  int      i, k;
  vector2f *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2f*)pkv_GetScratchMem ( 5*k*sizeof(vector2f) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];  tabpuu = &tabpv[k];  tabpuv = &tabpuu[k];  tabpvv = &tabpuv[k];

  if ( !mbs_TabBezC1CoonsDer2f ( 2, nkn, kn, hfunc, dhfunc, ddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc,
           G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
           G1_CROSS10DEG, (float*)c10, G1_CROSS11DEG, (float*)c11,
           G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
           G1_CROSS10DEG, (float*)d10, G1_CROSS11DEG, (float*)d11,
           NULL, (float*)tabpu, (float*)tabpv, (float*)tabpuu,
           (float*)tabpuv, (float*)tabpvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    _g1h_DiJacobian2f ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                        &jac[i], &trd[5*i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabDiPatchJac2f*/

boolean _g1h_TabLaplacianf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *fc00, const float *fc01,
        const float *fc10, const float *fc11,
        const float *fd00, const float *fd01,
        const float *fd10, const float *fd11,
        const float *trd,
        float *lap )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];

  if ( !mbs_TabBezC1CoonsDer2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 5 )
    lap[i] = trd[k+0]*tabfu[i] + trd[k+1]*tabfv[i] +
             trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] + trd[k+4]*tabfvv[i];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabLaplacianf*/

boolean _g1h_TabLaplacian0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *fc00, const float *fc01,
        const float *fd00, const float *fd01,
        const float *trd, float *lap )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];

  if ( !mbs_TabBezC1Coons0Der2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 5 )
    lap[i] = trd[k+0]*tabfu[i] + trd[k+1]*tabfv[i] +
             trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] + trd[k+4]*tabfvv[i];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabLaplacian0f*/

float _g1h_Integralf ( int hole_k, int nknsq, float *jac,
                       unsigned short supp1, float *func1,
                       unsigned short supp2, float *func2 )
{
  int            i, k;
  double         s;
  unsigned short supp;
  float          *f1, *f2, *jc;

  supp = (unsigned short)(supp1 & supp2);
  s = 0.0;
  for ( k = 0; k < hole_k; k++ )
    if ( supp & (0x0001 << k) ) {
      jc = &jac[k*nknsq];
      f1 = &func1[k*nknsq];
      f2 = &func2[k*nknsq];
      for ( i = 0; i < nknsq; i++ )
        s += f1[i]*f2[i]*jc[i];
    }
  s /= (double)(nknsq);
  return (float)s;
} /*_g1h_Integralf*/

boolean g1h_ComputeFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float    *trd, *jac, *lap, *Amat, *Bmat;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  support_b = privateG->support_b;
  Amat = privateG1->Amat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
  Bmat = privateG1->Bmat = malloc ( nfunc_a*nfunc_b*sizeof(float) );
  if ( !Amat || !Bmat )
    goto failure;

  n = hole_k*G1_NQUADSQ;

  tkn = pkv_GetScratchMemf ( 25*G1_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 5*n );
  lap = pkv_GetScratchMemf ( (nfunc_a+1)*n );
  if ( !tkn || !jac || !trd || !lap )
    goto failure;
  hfunc = &tkn[G1_NQUAD];  dhfunc = &hfunc[4*G1_NQUAD];  ddhfunc = &dhfunc[4*G1_NQUAD];
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
  mbs_TabCubicHFuncDer2f ( 0.0, 1.0, G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc );

  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesf ( domain, i, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    _g1h_TabDiPatchJac2f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                           c00, c01, c10, c11, d00, d01, d10, d11,
                           &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*5] );
  }

  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, j, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabLaplacian0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            fc00, fc01, fd00, fd01, &trd[i*G1_NQUADSQ*5],
                            &lap[(j*hole_k+i)*G1_NQUADSQ] );
    }

  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      Amat[l] = _g1h_Integralf ( hole_k, G1_NQUADSQ, jac,
                                 0xFFFF, &lap[i*n], 0xFFFF, &lap[j*n] );
  }
  
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, i, &fc00, &fc01, &fc10, &fc11,
                                  &fd00, &fd01, &fd10, &fd11 );
        _g1h_TabLaplacianf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                             fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                             &trd[i*G1_NQUADSQ*5],
                             &lap[(nfunc_a*hole_k+i)*G1_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      Bmat[i*nfunc_b+j] = _g1h_Integralf ( hole_k, G1_NQUADSQ, jac,
                            0xFFFF, &lap[i*n], support_b[j], &lap[nfunc_a*n] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeFormMatrixf*/

