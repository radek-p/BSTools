/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLSFuncGradd ( int nkn, double *qcoeff,
              double **Nitabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab,
              double *func, double *grad )
{
  void   *sp;
  int    *nncpi;
  double f;
  int    i, j, k, n, el, fp, ncp, ind;

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
        g2mbl_SFuncGradRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind],
                              nncpi, mvcp, mvcpn, &ftab[el], &gtab[ind] );
      else {
        if ( !g2mbl_SFuncGradSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                                    &domelcpind[ind],
                                    nncpi, mvcp, mvcpn, &ftab[el], &gtab[ind] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvcp*sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      k = domelcpind[fp+j];
      n = nncpi[k];
      if ( n >= 0 )
        grad[n] += gtab[fp+j];
    }
  }
  *func = f;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLSFuncGradd*/

boolean g2mbl_MLSFuncGradHessianAd ( int nkn, double *qcoeff,
              double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab, double *htab,
              double *func, double *grad,
              int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl )
{
  void     *sp;
  int      i, j, k, l, m, n, p, b, s, t, hi, el, ind, fp, ncp, hti, *nncpi;
  double   f;

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
        g2mbl_SFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp, mvcpn,
                                  &ftab[el], &gtab[ind], &htab[hti] );
      else {
        if ( !g2mbl_SFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp, mvcpn,
                          &ftab[el], &gtab[ind], &htab[hti] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvcp*sizeof(double) );
  for ( i = 0; i < nHbl; i++ )
    memset ( &Hbl[tHbl[i]], 0, sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    hti = domelem[el].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = domelcpind[fp+j];
      n = nncpi[k];
      if ( n >= 0 ) {
        grad[n] += gtab[fp+j];
        for ( l = 0; l <= j; l++ ) {
          m = domelcpind[fp+l];
          p = nncpi[m];
          if ( p >= 0 ) {
            hi = hti + pkn_SymMatIndex ( j, l );
            if ( n > p ) {  /* (n,p) is below the diagonal */
              s = iHbl[n].firsthbl;
              t = s + iHbl[n].nhbl-1;
              do {
                b = (s+t)/2;
                if ( cHbl[b] > p )
                  t = b;
                else
                  s = b;
              } while ( t-s > 1 );
              b = tHbl[s];
            }
            else {
              if ( n < p ) {  /* (p,n) is below the diagonal */
                s = iHbl[p].firsthbl;
                t = s + iHbl[p].nhbl-1;
                do {
                  b = (s+t) / 2;
                  if ( cHbl[b] > n )
                    t = b;
                  else
                    s = b;
                } while ( t-s > 1 );
                b = tHbl[s];
              }
              else
                b = tHbl[iHbl[p].firsthbl+iHbl[p].nhbl-1];
            }
            Hbl[b] += htab[hi];
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
} /*g2mbl_MLSFuncGradHessianAd*/

boolean g2mbl_MLSFuncGradHessianBd ( int nkn, double *qcoeff,
              double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab, double *htab,
              double *func, double *grad,
              int hsize, int *hprof, double **hrows )
{
  void     *sp;
  int      i, j, k, l, m, n, p, hi, el, ind, fp, ncp, hti, *nncpi;
  double   f;

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
        g2mbl_SFuncGradHessRSQd ( nkn, qcoeff,
                                   Nitabs[1], Nijtabs[1], Mijtabs[1],
                                   &domelcpind[ind], nncpi, mvcp, mvcpn,
                                   &ftab[el], &gtab[ind], &htab[hti] );
      else {
        if ( !g2mbl_SFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp, mvcpn,
                          &ftab[el], &gtab[ind], &htab[hti] ) )
          goto failure;
      }
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvcp*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    hti = domelem[el].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = domelcpind[fp+j];
      n = nncpi[k];
      if ( n >= 0 ) {
        grad[n] += gtab[fp+j];
        for ( l = 0; l <= j; l++ ) {
          m = domelcpind[fp+l];
          p = nncpi[m];
          if ( p >= 0 ) {
            hi = hti + pkn_SymMatIndex ( j, l );
            if ( n >= p )      /* (n,p) is below or on the diagonal */
              hrows[n][p] += htab[hi];
            else if ( n < p )  /* (p,n) is below the diagonal */
              hrows[p][n] += htab[hi];
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
} /*g2mbl_MLSFuncGradHessianBd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLSGetHessianRowsd ( int nv, int nvcp1, int *vncpi1,
                                   int nHbl, nzHbl_rowdesc *iHbl,
                                   int *cHbl, int *tHbl, double *Hbl,
                                   int nvcp, int *vncpi,
                                   int hsize, int *hprof, double **hrows )
{
  void *sp;
  int  *nncpi;
  int  i, j, k, l, fhbl, nhbl;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;

        /* copy the coefficients */
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0;  i < nvcp1; i++ )
    if ( (k = nncpi[vncpi1[i]]) >= 0 ) {
      fhbl = iHbl[i].firsthbl;
      nhbl = iHbl[i].nhbl;
      for ( j = 0; j < nhbl; j++ )
        if ( (l = nncpi[vncpi1[cHbl[fhbl+j]]]) >= 0 ) {
          if ( k >= l )
            hrows[k][l] = Hbl[tHbl[fhbl+j]];
          else
            hrows[l][k] = Hbl[tHbl[fhbl+j]];
        }
    }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLSGetHessianRowsd*/

