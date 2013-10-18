
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean _g1h_DiJacobian2d ( const vector2d *du, const vector2d *dv,
                            const vector2d *duu, const vector2d *duv,
                            const vector2d *dvv,
                            double *jac, double *trd )
{
  vector2d gx, gy, gxx, gxy, gyy;
  double   A21[6], A22[9];

  *jac = (double)fabs ( du->x*dv->y - du->y*dv->x );

  if ( !pkn_f2iDerivatives2d ( du->x, du->y, dv->x, dv->y,
          duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
          (double*)&gx, (double*)&gy, (double*)&gxx, (double*)&gxy, (double*)&gyy ) )
    return false;

  pkn_Setup2DerA21Matrixd ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
  pkn_Setup2DerA22Matrixd ( gx.x, gx.y, gy.x, gy.y, A22 );

  trd[0]  = A21[0]+A21[4];   trd[1]  = A21[1]+A21[5];
  trd[2]  = A22[0]+A22[6];   trd[3]  = A22[1]+A22[7];   trd[4] = A22[2]+A22[8];
  return true;
} /*_g1h_DiJacobian2d*/

boolean _g1h_TabDiPatchJac2d ( int nkn, const double *kn, const double *hfunc,
                               const double *dhfunc, const double *ddhfunc,
                               const vector2d *c00, const vector2d *c01,
                               const vector2d *c10, const vector2d *c11,
                               const vector2d *d00, const vector2d *d01,
                               const vector2d *d10, const vector2d *d11,
                               double *jac, double *trd )
{
  void     *sp;
  int      i, k;
  vector2d *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2d*)pkv_GetScratchMem ( 5*k*sizeof(vector2d) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];  tabpuu = &tabpv[k];  tabpuv = &tabpuu[k];  tabpvv = &tabpuv[k];

  if ( !mbs_TabBezC1CoonsDer2d ( 2, nkn, kn, hfunc, dhfunc, ddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc,
           G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
           G1_CROSS10DEG, (double*)c10, G1_CROSS11DEG, (double*)c11,
           G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
           G1_CROSS10DEG, (double*)d10, G1_CROSS11DEG, (double*)d11,
           NULL, (double*)tabpu, (double*)tabpv, (double*)tabpuu,
           (double*)tabpuv, (double*)tabpvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    _g1h_DiJacobian2d ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                        &jac[i], &trd[5*i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabDiPatchJac2d*/

boolean _g1h_TabLaplaciand ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *fc00, const double *fc01,
        const double *fc10, const double *fc11,
        const double *fd00, const double *fd01,
        const double *fd10, const double *fd11,
        const double *trd,
        double *lap )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];

  if ( !mbs_TabBezC1CoonsDer2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
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
} /*_g1h_TabLaplaciand*/

boolean _g1h_TabLaplacian0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *fc00, const double *fc01,
        const double *fd00, const double *fd01,
        const double *trd, double *lap )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];  tabfuu = &tabfv[k];  tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];

  if ( !mbs_TabBezC1Coons0Der2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
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
} /*_g1h_TabLaplacian0d*/

double _g1h_Integrald ( int hole_k, int nknsq, double *jac,
                        unsigned short supp1, double *func1,
                        unsigned short supp2, double *func2 )
{
  int            i, k;
  long double    s;
  unsigned short supp;
  double         *f1, *f2, *jc;

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
  return (double)s;
} /*_g1h_Integrald*/

boolean g1h_ComputeFormMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *trd, *jac, *lap, *Amat, *Bmat;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisd ( domain ) )
      goto failure;
  }
  hole_k  = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  support_b = privateG->support_b;
  Amat = privateG1->Amat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
  Bmat = privateG1->Bmat = malloc ( nfunc_a*nfunc_b*sizeof(double) );
  if ( !Amat || !Bmat )
    goto failure;

  n = hole_k*G1_NQUADSQ;

  tkn = pkv_GetScratchMemd ( 25*G1_NQUAD );
  jac = pkv_GetScratchMemd ( n );
  trd = pkv_GetScratchMemd ( 5*n );
  lap = pkv_GetScratchMemd ( (nfunc_a+1)*n );
  if ( !tkn || !jac || !trd || !lap )
    goto failure;
  hfunc = &tkn[G1_NQUAD];  dhfunc = &hfunc[4*G1_NQUAD];  ddhfunc = &dhfunc[4*G1_NQUAD];
  _gh_PrepareTabKnotsd ( G1_NQUAD, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesd ( domain, i, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    _g1h_TabDiPatchJac2d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                           c00, c01, c10, c11, d00, d01, d10, d11,
                           &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*5] );
  }

  for ( j = 0; j < nfunc_a; j++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabLaplacian0d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            fc00, fc01, fd00, fd01, &trd[i*G1_NQUADSQ*5],
                            &lap[(j*hole_k+i)*G1_NQUADSQ] );
    }

  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      Amat[l] = _g1h_Integrald ( hole_k, G1_NQUADSQ, jac,
                                 0xFFFF, &lap[i*n], 0xFFFF, &lap[j*n] );
  }
  
  for ( j = 0; j < nfunc_b; j++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[j] & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesd ( domain, j, i, &fc00, &fc01, &fc10, &fc11,
                                  &fd00, &fd01, &fd10, &fd11 );
        _g1h_TabLaplaciand ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                             fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                             &trd[i*G1_NQUADSQ*5],
                             &lap[(nfunc_a*hole_k+i)*G1_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      Bmat[i*nfunc_b+j] = _g1h_Integrald ( hole_k, G1_NQUADSQ, jac,
                            0xFFFF, &lap[i*n], support_b[j], &lap[nfunc_a*n] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeFormMatrixd*/

