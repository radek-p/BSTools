/* /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ */
double g2mbl_MLFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                       int nv, point3d *mvcp, int ndomel, int *domelind,
                       meshdom_elem *domelem, int *domelcpind,
                       boolean recalc, double *ftab )
{
  double f;
  int    i, k, el, ind;

        /* compute the integrals */
  if ( recalc )
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind  = domelem[el].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind], mvcp,
                          domelem[el].C, &ftab[el] );
      else
        g2mbl_UFuncSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind], mvcp,
                          domelem[el].C, &ftab[el] );
    }
        /* sum the integrals over the elements */
  f = 0.0;
  for ( i = 0; i < ndomel; i++ )
    f += ftab[domelind[i]];
  return f;
} /*g2mbl_MLFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                            int nv, point3d *mvcp, int ndomel, int *domelind,
                            meshdom_elem *domelem, int *domelcpind,
                            int nvcp, int *vncpi,
                            boolean recalc, double *ftab, double *gtab,
                            double *func, double *grad )
{
  void   *sp;
  int    *nncpi;
  double f;
  int    i, j, k, el, fp, ncp, ind;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
  if ( recalc ) {
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind = domelem[el].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncGradRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind],
                              nncpi, mvcp, domelem[el].C,
                              &ftab[el], &gtab[3*ind] );
      else {
        if ( !g2mbl_UFuncGradSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                                    &domelcpind[ind],
                                    nncpi, mvcp, domelem[el].C,
                                    &ftab[el], &gtab[3*ind] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, 3*nvcp*sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      k = 3*nncpi[domelcpind[fp+j]];
      if ( k >= 0 )
        AddVector3d ( (point3d*)&grad[k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[k] );
    }
  }
  *func = f;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLFuncGradd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLFuncGradHessianAd ( int nkn, double *qcoeff,
             double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
             int nv, point3d *mvcp, int ndomel, int *domelind,
             meshdom_elem *domelem, int *domelcpind,
             int nvcp, int *vncpi,
             boolean recalc, double *ftab, double *gtab, double *htab,
             double *func, double *grad,
             int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl )
{
  void   *sp;
  int    i, j, k, l, m, b, s, t, hi, el, ind, fp, ncp, hti, *nncpi;
  double f;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
  if ( recalc ) {
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind = domelem[el].firstcpi;
      hti = domelem[el].hti;
      if ( k == 4 )
        g2mbl_UFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp,
                                  domelem[el].C,
                                  &ftab[el], &gtab[3*ind], &htab[hti] );
      else {
        if ( !g2mbl_UFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp,
                          domelem[el].C,
                          &ftab[el], &gtab[3*ind], &htab[hti] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, 3*nvcp*sizeof(double) );
  for ( i = 0; i < nHbl; i++ )
    memset ( &Hbl[tHbl[i]], 0, 9*sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    hti = domelem[el].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = nncpi[domelcpind[fp+j]];
      if ( k >= 0 ) {
        AddVector3d ( (point3d*)&grad[3*k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[3*k] );
        for ( l = 0; l <= j; l++ ) {
          m = nncpi[domelcpind[fp+l]];
          if ( m >= 0 ) {
            hi = hti + 9*pkn_SymMatIndex ( j, l );
            if ( k > m ) {        /* (k,m) is below the diagonal */
              s = iHbl[k].firsthbl;
              t = s + iHbl[k].nhbl-1;
              do {
                b = (s+t)/2;
                if ( cHbl[b] > m )
                  t = b;
                else
                  s = b;
              } while ( t-s > 1 );
              b = tHbl[s];
            }
            else if ( k < m ) {   /* (m,k) is below the diagonal */
              s = iHbl[m].firsthbl;
              t = s + iHbl[m].nhbl-1;
              do {
                b = (s+t) / 2;
                if ( cHbl[b] > k )
                  t = b;
                else
                  s = b;
              } while ( t-s > 1 );
              b = tHbl[s];
            }
            else                  /* a diagonal block */
              b = tHbl[iHbl[k].firsthbl+iHbl[k].nhbl-1];
            pkn_AddMatrixd ( 1, 9, 0, &Hbl[b], 0, &htab[hi], 0, &Hbl[b] );
          }
        }
      }
    }
  }
  *func = f;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLFuncGradHessianAd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLFuncGradHessianBd ( int nkn, double *qcoeff,
             double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
             int nv, point3d *mvcp, int ndomel, int *domelind,
             meshdom_elem *domelem, int *domelcpind,
             int nvcp, int *vncpi,
             boolean recalc, double *ftab, double *gtab, double *htab,
             double *func, double *grad,
             int hsize, int *hprof, double **hrows )
{
  void   *sp;
  int    i, j, k, l, m, b, s, t, el, ind, fp, ncp, hti, *nncpi;
  double f;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
  if ( recalc ) {
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind = domelem[el].firstcpi;
      hti = domelem[el].hti;
      if ( k == 4 )
        g2mbl_UFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp,
                                  domelem[el].C,
                                  &ftab[el], &gtab[3*ind], &htab[hti] );
      else {
        if ( !g2mbl_UFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp,
                          domelem[el].C,
                          &ftab[el], &gtab[3*ind], &htab[hti] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, 3*nvcp*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    hti = domelem[el].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = 3*nncpi[domelcpind[fp+j]];
      if ( k >= 0 ) {
        AddVector3d ( (point3d*)&grad[k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[k] );
        for ( l = 0; l <= j; l++ ) {
          m = 3*nncpi[domelcpind[fp+l]];
          if ( m >= 0 ) {
            b = hti + 9*pkn_SymMatIndex ( j, l );
            if ( k > m ) {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t < 3; t++ )
                  hrows[k+s][m+t] += htab[b+3*s+t];
            }
            else if ( k < m ) {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t < 3; t++ )
                  hrows[m+t][k+s] += htab[b+3*t+s];
            }
            else {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t <= s; t++ )
                  hrows[k+s][k+t] += htab[b+3*s+t];
            }
          }
        }
      }
    }
  }
  *func = f;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLFuncGradHessianBd*/

