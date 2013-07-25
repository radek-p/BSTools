
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */   
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"
#include "g2mblmlprivated.h"

/* /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ */
typedef struct {
    int          nkn;
    double       *qcoeff;
    double       **Nitabs, **Nijtabs, **Mijtabs, **Jac;
    int          nv;
    int          *nncpi;
    point3d      *mvcp;
    vector3d     *mvcpn;
    int          ndomel;
    int          *domelind;
    meshdom_elem *domelem;
    int          *domelcpind;
    double       *ftab, *gtab, *htab;
  } g2mbl_sjob_desc;

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2mbl_TMLSFuncGradd ( void *usrdata, int3 *jobnum )
{
  void            *sp;
  g2mbl_sjob_desc *data;
  int             i, k, el, ind;

  sp = pkv_GetScratchMemTop ();
  data = (g2mbl_sjob_desc*)usrdata;
  i = jobnum->x;

  el = data->domelind[i];
  k = data->domelem[el].type;
  ind = data->domelem[el].firstcpi;
  if ( k == 4 )
    g2mbl_SFuncGradRSQd ( data->nkn, data->qcoeff, data->Nitabs[1],
                          &data->domelcpind[ind],
                          data->nncpi, data->mvcp, data->mvcpn,
                          &data->ftab[el], &data->gtab[ind] );
  else {
    if ( !g2mbl_SFuncGradSSQd ( data->nkn, data->qcoeff, k, data->Nitabs[k-3],
                                data->Jac[k-3], &data->domelcpind[ind],
                                data->nncpi, data->mvcp, data->mvcpn,
                                &data->ftab[el], &data->gtab[ind] ) )
      goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_TMLSFuncGradd*/

boolean g2mbl_MLSFuncGradd ( int nkn, double *qcoeff,
              double **Nitabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab,
              double *func, double *grad )
{
  void            *sp;
  g2mbl_sjob_desc data;
  int3            size;
  int             *nncpi;
  double          f;
  int             i, j, k, n, el, fp, ncp;
  boolean         success;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* compute the integrals */
  if ( recalc ) {
          /* copy the parameters */
    data.nkn        = nkn;
    data.qcoeff     = qcoeff;
    data.Nitabs     = Nitabs;
    data.Jac        = Jac;
    data.nv         = nv;
    data.nncpi      = nncpi;
    data.mvcp       = mvcp;
    data.mvcpn      = mvcpn;
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    size.y = size.z = 1;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( &size, _g2mbl_npthreads, 1048576, 16*1048576,
                             (void*)&data, _g2mbl_TMLSFuncGradd,
                             NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLSFuncGradd ( (void*)&data, &size );
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

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2mbl_TMLSFuncGradHessiand ( void *usrdata, int3 *jobnum )
{
  void            *sp;
  g2mbl_sjob_desc *data;
  int             i, k, el, ind, hti;

  sp = pkv_GetScratchMemTop ();
  data = (g2mbl_sjob_desc*)usrdata;
  i = jobnum->x;

  el = data->domelind[i];
  k = data->domelem[el].type;
  ind = data->domelem[el].firstcpi;
  hti = data->domelem[el].hti;
  if ( k == 4 )
    g2mbl_SFuncGradHessRSQd ( data->nkn, data->qcoeff,
                              data->Nitabs[1], data->Nijtabs[1], data->Mijtabs[1],
                              &data->domelcpind[ind], data->nncpi,
                              data->mvcp, data->mvcpn,
                              &data->ftab[el], &data->gtab[ind], &data->htab[hti] );
  else {
    if ( !g2mbl_SFuncGradHessSSQd ( data->nkn, data->qcoeff, k,
                      data->Nitabs[k-3], data->Nijtabs[k-3],
                      data->Mijtabs[k-3], data->Jac[k-3],
                      &data->domelcpind[ind], data->nncpi, data->mvcp, data->mvcpn,
                      &data->ftab[el], &data->gtab[ind], &data->htab[hti] ) )
      goto failure;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_TMLSFuncGradHessiand*/

boolean g2mbl_MLSFuncGradHessianAd ( int nkn, double *qcoeff,
              double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab, double *htab,
              double *func, double *grad,
              int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl )
{
  void            *sp;
  g2mbl_sjob_desc data;
  int3            size;
  int             i, j, k, l, m, n, p, b, s, t, hi, el, fp, ncp, hti, *nncpi;
  double          f;
  boolean         success;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* compute the integrals */
  if ( recalc ) {
          /* copy the parameters */
    data.nkn        = nkn;
    data.qcoeff     = qcoeff;
    data.Nitabs     = Nitabs;
    data.Nijtabs    = Nijtabs;
    data.Mijtabs    = Mijtabs;
    data.Jac        = Jac;
    data.nv         = nv;
    data.nncpi      = nncpi;
    data.mvcp       = mvcp;
    data.mvcpn      = mvcpn;
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    data.htab       = htab;
    size.y = size.z = 1;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( &size, _g2mbl_npthreads, 1048576, 16*1048576,
                              (void*)&data, _g2mbl_TMLSFuncGradHessiand,
                              NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLSFuncGradHessiand ( (void*)&data, &size );
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
  void            *sp;
  g2mbl_sjob_desc data;
  int3            size;
  int             i, j, k, l, m, n, p, hi, el, fp, ncp, hti, *nncpi;
  double          f;
  boolean         success;

  sp = pkv_GetScratchMemTop ();
  nncpi = pkv_GetScratchMemi ( nv );
  if ( !nncpi )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* compute the integrals */
  if ( recalc ) {
          /* copy the parameters */
    data.nkn        = nkn;
    data.qcoeff     = qcoeff;
    data.Nitabs     = Nitabs;
    data.Nijtabs    = Nijtabs;
    data.Mijtabs    = Mijtabs;
    data.Jac        = Jac;
    data.nv         = nv;
    data.nncpi      = nncpi;
    data.mvcp       = mvcp;
    data.mvcpn      = mvcpn;
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    data.htab       = htab;
    size.y = size.z = 1;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( &size, _g2mbl_npthreads, 1048576, 16*1048576,
                              (void*)&data, _g2mbl_TMLSFuncGradHessiand,
                              NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLSFuncGradHessiand ( (void*)&data, &size );
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

