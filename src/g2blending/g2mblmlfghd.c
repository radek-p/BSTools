
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2015                            */
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
    int          ndomel;
    int          *domelind;
    meshdom_elem *domelem;
    int          *domelcpind;
    double       *ftab, *gtab, *htab;
  } g2mbl_job_desc;

int _g2mbl_npthreads = 1;

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2mbl_TMLFuncd ( void *usrdata, int4 *jobnum )
{
  g2mbl_job_desc *data;
  int            i, k, el, ind;

  data = (g2mbl_job_desc*)usrdata;
  i = jobnum->x;

  el = data->domelind[i];
  k = data->domelem[el].type;
  ind  = data->domelem[el].firstcpi;
  if ( k == 4 )
    g2mbl_UFuncRSQd ( data->nkn, data->qcoeff, data->Nitabs[1],
                      &data->domelcpind[ind], data->mvcp,
                      data->domelem[el].C, &data->ftab[el] );
  else
    g2mbl_UFuncSSQd ( data->nkn, data->qcoeff, k, data->Nitabs[k-3],
                      data->Jac[k-3], &data->domelcpind[ind], data->mvcp,
                      data->domelem[el].C, &data->ftab[el] );
  return true;
} /*_g2mbl_TMLFuncd*/

double g2mbl_MLFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                       int nv, point3d *mvcp, int ndomel, int *domelind,
                       meshdom_elem *domelem, int *domelcpind,
                       boolean recalc, double *ftab )
{
  g2mbl_job_desc data;
  int4           size;
  double         f;
  int            i;
  boolean        success;

        /* compute the integrals */
  if ( recalc ) {
         /* copy the parameters */
    data.nkn        = nkn;
    data.qcoeff     = qcoeff;
    data.Nitabs     = Nitabs;
    data.Jac        = Jac;
    data.nv         = nv;
    data.mvcp       = mvcp;
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    if ( _g2mbl_npthreads > 1 ) {
        /* setup the number of jobs */
      size.x = ndomel;
        /* set threads to work */
      pkv_SetPThreadsToWork ( 1, &size, _g2mbl_npthreads, 1048576, 16*1048576,
                              (void*)&data, _g2mbl_TMLFuncd,
                              NULL, NULL, &success );
    }
    else {
        /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLFuncd ( (void*)&data, &size );
    }
  }
        /* sum the integrals over the elements */
  f = 0.0;
  for ( i = 0; i < ndomel; i++ )
    f += ftab[domelind[i]];
  return f;
} /*g2mbl_MLFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2mbl_TMLFuncGradd ( void *usrdata, int4 *jobnum )
{
  void           *sp;
  g2mbl_job_desc *data;
  int            i, k, el, ind;

  sp = pkv_GetScratchMemTop ();
  data = (g2mbl_job_desc*)usrdata;
  i = jobnum->x;

  el = data->domelind[i];
  k = data->domelem[el].type;
  ind = data->domelem[el].firstcpi;
  if ( k == 4 )
    g2mbl_UFuncGradRSQd ( data->nkn, data->qcoeff, data->Nitabs[1],
                          &data->domelcpind[ind],
                          data->nncpi, data->mvcp, data->domelem[el].C,
                          &data->ftab[el], &data->gtab[3*ind] );
  else {
    if ( !g2mbl_UFuncGradSSQd ( data->nkn, data->qcoeff, k, data->Nitabs[k-3],
                                data->Jac[k-3], &data->domelcpind[ind],
                                data->nncpi, data->mvcp, data->domelem[el].C,
                                &data->ftab[el], &data->gtab[3*ind] ) )
      goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_TMLFuncGradd*/

boolean g2mbl_MLFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                            int nv, point3d *mvcp, int ndomel, int *domelind,
                            meshdom_elem *domelem, int *domelcpind,
                            int nvcp, int *vncpi,
                            boolean recalc, double *ftab, double *gtab,
                            double *func, double *grad )
{
  void           *sp;
  g2mbl_job_desc data;
  int4           size;
  int            *nncpi;
  double         f;
  int            i, j, k, el, fp, ncp;
  boolean        success;

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
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( 1, &size, _g2mbl_npthreads, 1048576, 16*1048576,
                             (void*)&data, _g2mbl_TMLFuncGradd,
                             NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLFuncGradd ( (void*)&data, &size );
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
static boolean _g2mbl_TMLFuncGradHessiand ( void *usrdata, int4 *jobnum )
{
  void           *sp;
  g2mbl_job_desc *data;
  int            i, k, el, ind, hti;

  sp = pkv_GetScratchMemTop ();
  data = (g2mbl_job_desc*)usrdata;
  i = jobnum->x;

  el = data->domelind[i];
  k = data->domelem[el].type;
  ind = data->domelem[el].firstcpi;
  hti = data->domelem[el].hti;
  if ( k == 4 )
    g2mbl_UFuncGradHessRSQd ( data->nkn, data->qcoeff, data->Nitabs[1],
                              data->Nijtabs[1], data->Mijtabs[1],
                              &data->domelcpind[ind], data->nncpi,
                              data->mvcp, data->domelem[el].C, &data->ftab[el],
                              &data->gtab[3*ind], &data->htab[hti] );
  else {
    if ( !g2mbl_UFuncGradHessSSQd ( data->nkn, data->qcoeff, k,
                      data->Nitabs[k-3], data->Nijtabs[k-3],
                      data->Mijtabs[k-3], data->Jac[k-3],
                      &data->domelcpind[ind], data->nncpi, data->mvcp,
                      data->domelem[el].C,
                      &data->ftab[el], &data->gtab[3*ind], &data->htab[hti] ) )
      goto failure;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_TMLFuncGradHessiand*/

boolean g2mbl_MLFuncGradHessianAd ( int nkn, double *qcoeff,
             double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
             int nv, point3d *mvcp, int ndomel, int *domelind,
             meshdom_elem *domelem, int *domelcpind,
             int nvcp, int *vncpi,
             boolean recalc, double *ftab, double *gtab, double *htab,
             double *func, double *grad,
             int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl )
{
  void           *sp;
  g2mbl_job_desc data;
  int4           size;
  int            i, j, k, l, m, b, s, t, hi, el, fp, ncp, hti, *nncpi;
  double         f;
  boolean        success;

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
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    data.htab       = htab;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( 1, &size, _g2mbl_npthreads, 1048576, 16*1048576,
                             (void*)&data, _g2mbl_TMLFuncGradHessiand,
                             NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLFuncGradHessiand ( (void*)&data, &size );
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
  g2mbl_job_desc data;
  int4           size;
  int            i, j, k, l, m, b, s, t, el, fp, ncp, hti, *nncpi;
  double         f;
  boolean        success;

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
    data.ndomel     = ndomel;
    data.domelind   = domelind;
    data.domelem    = domelem;
    data.domelcpind = domelcpind;
    data.ftab       = ftab;
    data.gtab       = gtab;
    data.htab       = htab;
    if ( _g2mbl_npthreads > 1 ) {
          /* setup the number of jobs */
      size.x = ndomel;
          /* set threads to work */
      pkv_SetPThreadsToWork ( 1, &size, _g2mbl_npthreads, 1048576, 16*1048576,
                             (void*)&data, _g2mbl_TMLFuncGradHessiand,
                             NULL, NULL, &success );
    }
    else {
          /* no concurrent threads - do it yourself */
      for ( size.x = 0; size.x < ndomel; size.x++ )
        _g2mbl_TMLFuncGradHessiand ( (void*)&data, &size );
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

/* /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ */
boolean g2mbl_MLGetHessianRowsd ( int nv, int nvcp1, int *vncpi1,
                                  int nHbl, nzHbl_rowdesc *iHbl,
                                  int *cHbl, int *tHbl, double *Hbl,
                                  int nvcp, int *vncpi,
                                  int hsize, int *hprof, double **hrows )
{
  void *sp;
  int  *nncpi;
  int  i, j, k, l, b, s, t, fhbl, nhbl;

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
  for ( i = 0; i < nvcp1; i++ )
    if ( (k = nncpi[vncpi1[i]]) >= 0 ) {
      fhbl = iHbl[i].firsthbl;
      nhbl = iHbl[i].nhbl;
      for ( j = 0; j < nhbl; j++ )
        if ( (l = nncpi[vncpi1[cHbl[fhbl+j]]]) >= 0 ) {
          b = tHbl[fhbl+j];
          if ( k > l ) {
            for ( s = 0; s < 3; s++ )
              for ( t = 0; t < 3; t++ )
                hrows[3*k+s][3*l+t] = Hbl[b+3*s+t];
          }
          else if ( k < l ) {
            for ( s = 0; s < 3; s++ )
              for ( t = 0; t < 3; t++ )
                hrows[3*l+t][3*k+s] = Hbl[b+3*s+t];
          }
          else {
            for ( s = 0; s < 3; s++ )
              for ( t = 0; t <= s; t++ )
                hrows[3*k+s][3*k+t] = Hbl[b+3*s+t];
          }
        }
    }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLGetHessianRowsd*/

