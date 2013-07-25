
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#include <sys/times.h>
#include <unistd.h>

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
#include "msgpool.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_CMPSSetupCoarseHessiand ( mesh_ml_optdata *d, int bl, double nu )
{
  void         *sp, *sp1;
  mlblock_desc *bd;
  double       *Hbl, *rmnzc;
  nzHbl_rowdesc *iHbl;
  int          brmnnz, *cHbl, *tHbl, *rrows, *hrcols;
  unsigned int *rpermut, *hrpermut, *auxpermut;
  index3       *brmnzi, *bhnzi;
  index2       *hri, *aikbkj, hind;
  int          hsize;
  double       **hrows, *hrc, *hnu, s;
  int          bnvcp, bnwcp, nnz1, bnnzH;
  int          wsps, wspn, *wsp;
  int          i, j, j0, j1, k, k0, k1, l, l0, l1, m, n, cnt;

  sp = pkv_GetScratchMemTop ();
  bd = &d->bd[bl];
  bhnzi = NULL;
  hnu = NULL;
  hri = NULL;
  hrc = NULL;
  rpermut = NULL;
  if ( bd->rmnzi ) {
    Hbl   = d->Hbl;
    rmnzc = d->rmnzc;
        /* get block description data */
    brmnnz = bd->rmnnz;
    brmnzi = bd->rmnzi;
    hsize  = bd->hsize;
    hrows  = bd->hrows;
    iHbl   = bd->iHbl;
    cHbl   = bd->cHbl;
    tHbl   = bd->tHbl;
    bnvcp  = bd->nvcp;
    bnwcp  = bd->nwcp;
    nnz1   = bd->nnz1;
        /* find the sparse Hessian representation */
    bnnzH = 2*bd->nHbl - bnvcp;
    PKV_MALLOC ( bhnzi, bnnzH*sizeof(index3) );
    PKV_MALLOC ( hnu, bnnzH*sizeof(double) );
    if ( !bhnzi || !hnu )
      goto failure;
    for ( i = l = 0; i < bnvcp; i++ ) {
      j0 = iHbl[i].firsthbl;
      j1 = j0 + iHbl[i].nhbl - 1;
      for ( j = j0; j < j1; j++ ) {
        k = cHbl[j];
        bhnzi[l].i = i;
        bhnzi[l].j = k;
        bhnzi[l].k = l;
        hnu[l++] = Hbl[tHbl[j]];
        bhnzi[l].i = k;
        bhnzi[l].j = i;
        bhnzi[l].k = l;
        hnu[l++] = Hbl[tHbl[j]];
      }
      bhnzi[l].i = bhnzi[l].j = i;
      bhnzi[l].k = l;
      hnu[l++] = Hbl[tHbl[j]] + nu;
    }
        /* multiply the Hessian by the refinement matrix */
    PKV_MALLOC ( hri, nnz1*sizeof(index2) );
    PKV_MALLOC ( hrc, nnz1*sizeof(double) );
    if ( !hri || !hrc )
      goto failure;
    if ( !pkn_SPsubMmultMMCd ( bnvcp, bnvcp, bnwcp,
                               bnnzH, bhnzi, hnu, NULL, NULL, false,
                               brmnnz, brmnzi, rmnzc, NULL, NULL, false,
                               hri, hrc ) )
      goto failure;

        /* multiply the transposed refinement matrix by the product */
        /* obtained above */
    PKV_MALLOC ( rpermut, (nnz1+brmnnz+bnvcp+bnwcp+2)*sizeof(int) );
    if ( !rpermut )
      goto failure;
    rrows = (int*)&rpermut[brmnnz];
    hrpermut = (unsigned int*)&rrows[bnvcp+1];
    hrcols = (int*)&hrpermut[nnz1];
    if ( !pkn_SPsubMFindRows ( bnvcp, bnwcp, brmnnz, brmnzi,
                               rpermut, false, rrows ) )
      goto failure;
    if ( !pkn_SPMFindCols ( bnvcp, bnwcp, nnz1, hri,
                            hrpermut, false, hrcols ) )
      goto failure;

    memset ( hrows[0], 0, hsize*sizeof(double) );
    sp1 = pkv_GetScratchMemTop ();
    wsp = NULL;
    aikbkj = NULL;
    wsps = 0;
    for ( i = 0; i < bnwcp; i++ ) {
      k0 = hrcols[i];
      k1 = hrcols[i+1];
      wspn = 0;
      for ( k = k0; k < k1; k++ ) {
        j = hri[hrpermut[k]].i;
        wspn += rrows[j+1]-rrows[j];
      }
      if ( wspn > wsps ) {
        pkv_SetScratchMemTop ( sp1 );
        wsp = pkv_GetScratchMemi ( 4*wspn );
        if ( !wsp )
          goto failure;
        aikbkj = (index2*)&wsp[wspn];
        wsps = wspn;
      }
      cnt = 0;
      for ( k = k0; k < k1; k++ ) {
        m = hrpermut[k];
        j = hri[m].i;
        l0 = rrows[j];
        l1 = rrows[j+1];
        for ( l = l0; l < l1; l++ ) {
          n = rpermut[l];
          if ( brmnzi[n].j >= i ) {
            aikbkj[cnt].i = brmnzi[n].k;
            aikbkj[cnt].j = m;
            wsp[cnt++] = brmnzi[n].j;
          }
        }
      }
      if ( cnt ) {
        auxpermut = (unsigned int*)&aikbkj[wspn];
        for ( l = 0; l < cnt; l++ )
          auxpermut[l] = l;
        if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(int),
                              0, cnt, wsp, auxpermut ) != SORT_OK )
          goto failure;
        k = auxpermut[0];
        hind.i = l0 = wsp[k];
        hind.j = i;
        s = rmnzc[aikbkj[k].i]*hrc[aikbkj[k].j];
        for ( l = 1; l < cnt; l++ ) {
          k = auxpermut[l];
          if ( wsp[k] != l0 ) {
#ifdef DEBUG
if ( hind.i < hind.j || hind.j < hprof[hind.i] ) {
  printf ( "blad profilu!" );
  exit ( 1 );
}
#endif
            hrows[hind.i][hind.j] = s;
            s = rmnzc[aikbkj[k].i]*hrc[aikbkj[k].j];
            hind.i = l0 = wsp[k];
          }
          else
            s += rmnzc[aikbkj[k].i]*hrc[aikbkj[k].j];
        }
#ifdef DEBUG
if ( hind.i < hind.j || hind.j < hprof[hind.i] ) {
  printf ( "blad profilu!" );
  exit ( 1 );
}
#endif
        hrows[hind.i][hind.j] = s;
      }
    }
  }
  PKV_FREE ( bhnzi );
  PKV_FREE ( hnu );
  PKV_FREE ( hri );
  PKV_FREE ( hrc );
  PKV_FREE ( rpermut );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( bhnzi )   PKV_FREE ( bhnzi );
  if ( hnu )     PKV_FREE ( hnu );
  if ( hri )     PKV_FREE ( hri );
  if ( hrc )     PKV_FREE ( hrc );
  if ( rpermut ) PKV_FREE ( rpermut );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPSSetupCoarseHessiand*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLSmultAxd ( int nvars, void *usrdata, double *x, double *Ax )
{
  mesh_ml_cg_data *cgd;
  mesh_ml_optdata *d;
  mlblock_desc    *bd;
  int             bl;
  nzHbl_rowdesc   *iHbl;
  int             *cHbl, *tHbl, i, j, j0, j1;
  double          *Hbl, *bH, *bx1, *bAx1, *bx2, *bAx2, nu;
  
  cgd = (mesh_ml_cg_data*)usrdata;
  d = cgd->d;
  Hbl  = d->Hbl;
  bl   = cgd->bl;
  bd   = &d->bd[bl];
  /*nvcp = bd->nvcp;*/  /* must be equal to nvars */
  iHbl = bd->iHbl;
  cHbl = bd->cHbl;
  tHbl = bd->tHbl;
  memset ( Ax, 0, nvars*sizeof(double) );
  for ( i = 0; i < nvars; i++ ) {
    bAx1 = &Ax[i];
    bx2 = &x[i];
    j0 = iHbl[i].firsthbl;
    j1 = j0 + iHbl[i].nhbl-1;
    for ( j = j0; j < j1; j++ ) {
      bH = &Hbl[tHbl[j]];
      bAx2 = &Ax[cHbl[j]];
      bx1 = &x[cHbl[j]];
        /* block below the diagonal */
      bAx1[0] += bH[0]*bx1[0];
        /* block above the diagonal */
      bAx2[0] += bH[0]*bx2[0];
    }
        /* diagonal block */
    bH = &Hbl[tHbl[j1]];
    bAx1[0] += bH[0]*bx2[0];
  }
  if ( (nu = cgd->nu) > 0.0 ) {
    for ( i = 0; i < nvars; i++ )
      Ax[i] += nu*x[i];
  }
  return true;
} /*_g2mbl_MLSmultAxd*/

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
    mesh_ml_cg_data *cgd;
    int             nvars, *nncpi;
    double          *x, *Qix;
    short int       bl, bl0, bl1;
    pthread_mutex_t mutex;
  } g2mbl_multqixdata;

static boolean _g2mbl_MLSmultQIxd4 ( void *usrdata, int3 *jobnum )
{
  boolean           *sp;
  g2mbl_multqixdata *qdata;
  mesh_ml_optdata   *d;
  mesh_ml_cg_data   *cgd;
  short int         bl;
  mlblock_desc      *bd;
  int               nvcp, *nncpi, *vncpi, *avncpi;
  int               i, *hprof;
  double            *x, *ax, **lhrows;

  sp = pkv_GetScratchMemTop ();
  qdata = (g2mbl_multqixdata*)usrdata;
  cgd = qdata->cgd;
  d = cgd->d;
  bl = qdata->bl0 + jobnum->x;
  bd = &d->bd[bl];
        /* get the right selection/permutation of the vertices */
  nvcp = bd->nvcp;
  avncpi = pkv_GetScratchMemi ( nvcp );
  if ( !avncpi )
    goto failure;
  nncpi = qdata->nncpi;
  vncpi = bd->vncpi;
  for ( i = 0; i < nvcp; i++ )
    avncpi[i] = nncpi[vncpi[i]];
        /* solve the system using the Cholesky's decomposition factor */
  ax = pkv_GetScratchMemd ( nvcp );
  if ( !ax )
    goto failure;
  x = qdata->x;
  for ( i = 0; i < nvcp; i++ )
    ax[i] = x[avncpi[i]];
  hprof = bd->hprof;
  lhrows = bd->lhrows;
  pkn_NRBLowerTrSolved ( nvcp, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
  pkn_NRBUpperTrSolved ( nvcp, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
        /* add the multiplication result */
  x = qdata->Qix;
  if ( _g2mbl_npthreads > 1 ) {
          /* with pthreads in use, mutual exclusion is needed  */
    pthread_mutex_lock ( &qdata->mutex );
    for ( i = 0; i < nvcp; i++ )
      x[avncpi[i]] += ax[i];
    pthread_mutex_unlock ( &qdata->mutex );
  }
  else {
    for ( i = 0; i < nvcp; i++ )
      x[avncpi[i]] += ax[i];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSmultQIxd4*/

static boolean _g2mbl_MLSmultQIxd5 ( void *usrdata, int3 *jobnum )
{
  boolean           *sp;
  g2mbl_multqixdata *qdata;
  mesh_ml_optdata   *d;
  mesh_ml_cg_data   *cgd;
  short int         bl;
  mlblock_desc      *bd;
  int               *hprof, nvcp, nwcp;
  double            **lhrows, *rmnzc, *x, *ax, *ay;
  int               rmnnz;
  index3            *rmnzi;

  sp = pkv_GetScratchMemTop ();
  qdata = (g2mbl_multqixdata*)usrdata;
  cgd = qdata->cgd;
  d = cgd->d;
  bl = cgd->cbl;
  bd = &d->bd[bl];
          /* compute the coarse mesh preconditioner term, if present */
  if ( bd->rmnzi ) {
    hprof  = bd->hprof;
    lhrows = bd->lhrows;
    nvcp   = bd->nvcp;
    nwcp   = bd->nwcp;
    rmnnz = bd->rmnnz;
    rmnzi = bd->rmnzi;
    rmnzc = d->rmnzc;
    ax = pkv_GetScratchMemd ( nwcp+nvcp );
    if ( !ax )
      goto failure;
    ay = &ax[nwcp];
    x = qdata->x;
    if ( !pkn_MultSPsubMTVectord ( nvcp, nwcp, rmnnz, rmnzi, rmnzc,
                                   1, x, ax ) )
      goto failure;
    pkn_NRBLowerTrSolved ( nwcp, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
    pkn_NRBUpperTrSolved ( nwcp, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
    if ( !pkn_MultSPsubMVectord ( nvcp, nwcp, rmnnz, rmnzi, rmnzc,
                                  1, ax, ay ) )
      goto failure;
    x = qdata->Qix;
    if ( _g2mbl_npthreads > 1 ) {
      pthread_mutex_lock ( &qdata->mutex );
      pkn_AddMatrixd ( 1, nvcp, 0, x, 0, ay, 0, x );
      pthread_mutex_unlock ( &qdata->mutex );
    }
    else
      pkn_AddMatrixd ( 1, nvcp, 0, x, 0, ay, 0, x );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSmultQIxd5*/

boolean _g2mbl_MLSmultQIxd ( int nvars, void *usrdata, double *x, double *Qix )
{
  void                *sp;
  g2mbl_multqixdata   qdata;
  mlblock_desc        *bd;
  mesh_ml_optdata     *d;
  int3                jobsize;
  boolean             result;
  short int           bl, nsmbl;
  int                 i, nvcp, *nncpi, *vncpi;
  pthread_mutexattr_t attr;

  sp = pkv_GetScratchMemTop ();
  qdata.cgd   = (mesh_ml_cg_data*)usrdata;
  d           = qdata.cgd->d;
  qdata.nvars = nvars;
  qdata.bl    = bl = qdata.cgd->cbl;
  qdata.x     = x;
  qdata.Qix   = Qix;
        /* find out, how many terms to compute */
  nsmbl = 0x0002;
  qdata.bl0 = 2*bl+1;
  bd = &d->bd[qdata.bl0];
  while ( bd->iHbl ) {
    nsmbl <<= 1;
    qdata.bl0 = 2*qdata.bl0+1;
    bd = &d->bd[qdata.bl0];
  }
  qdata.bl1 = qdata.bl0 + nsmbl - 1;
  jobsize.y = jobsize.z = 1;
        /* setup the array of numbers of vertices in the current block; */
        /* this will be used to select the right variables for the subblocks */
  nncpi = pkv_GetScratchMemi ( d->nv );
  if ( !nncpi )
    goto failure;
  qdata.nncpi = nncpi;
  bd = &d->bd[qdata.bl];
  nvcp = bd->nvcp;
  vncpi = bd->vncpi;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* do the multiplication */
  memset ( Qix, 0, nvars*sizeof(double) );
  if ( _g2mbl_npthreads > 1 ) {
          /* parallel computations */
    if ( pthread_mutexattr_init ( &attr ) )
      goto failure;
    if ( pthread_mutex_init ( &qdata.mutex, &attr ) )
      goto failure;
    pthread_mutexattr_destroy ( &attr );
    jobsize.x = nsmbl;
    result = pkv_SetPThreadsToWork ( &jobsize, _g2mbl_npthreads,
                                  4*1048576, 16*1048576,
                                  (void*)&qdata, _g2mbl_MLSmultQIxd4,
                                  (void*)&qdata, _g2mbl_MLSmultQIxd5,
                                  &result );
    pthread_mutex_destroy ( &qdata.mutex );
    if ( !result )
      goto failure;
  }
  else {
          /* sequential computations */
    result = true;
    for ( i = 0; i < nsmbl; i++ ) {
      jobsize.x = i;
      if ( !_g2mbl_MLSmultQIxd4 ( (void*)&qdata, &jobsize ) ) {
        result = false;
        break;
      }
    }
    bd = &d->bd[bl];
    if ( result && bd->rmnzi )
      result = _g2mbl_MLSmultQIxd5 ( (void*)&qdata, NULL );
  }
  pkv_SetScratchMemTop ( sp );
  return result;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSmultQIxd*/

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
    mesh_ml_optdata *d;
    short int       bl, bl0, bl1;  /* main, first small and last small block */
    double          nu;
    boolean         positive, abort;
  } g2mbl_blockdecompdata;

static boolean _g2mbl_MLSDecompSmallBlockd ( void *usrdata, int3 *jobnum )
{
  g2mbl_blockdecompdata *qdata;
  mesh_ml_optdata       *d;
  short int             bl;
  mlblock_desc          *bd0, *bd;
  double                nu;
  int                   nv, nvcp, *vncpi, i;
  int                   hsize, *hprof;
  double                **hrows, **lhrows;

        /* decompose the matrix for a small block */
  qdata = (g2mbl_blockdecompdata*)usrdata;
  d = qdata->d;
  nv = d->nv;
  nu = qdata->nu;
  bl = qdata->bl0 + jobnum->x;
  bd = &d->bd[bl];
  bd0 = &d->bd[qdata->bl];

  nvcp  = bd->nvcp;
  vncpi = bd->vncpi;
  hsize = bd->hsize;
  hprof = bd->hprof;
  hrows = bd->hrows;
  lhrows = bd->lhrows;
  if ( !(bd->fghflag & FLAG_H) ) {  /* invalid Hessian, recompute */
    if ( !g2mbl_MLSGetHessianRowsd ( nv, bd0->nvcp, bd0->vncpi,
                          bd0->nHbl, bd0->iHbl, bd0->cHbl, bd0->tHbl, d->Hbl,
                          nvcp, vncpi, hsize, hprof, hrows ) ) {
printf ( "%s\n", ERRMSG_17 );
      goto failure;
    }
    bd->fghflag |= FLAG_H;
    bd->fghflag &= ~FLAG_LH;
  }
  if ( !(bd->fghflag & FLAG_LH) ) {  /* invalid triangular factor, recompute */
    memcpy ( lhrows[0], hrows[0], hsize*sizeof(double) );
    if ( nu > 0.0 )
      for ( i = 0; i < nvcp; i++ )
        lhrows[i][i] += nu;
    if ( !pkn_NRBSymCholeskyDecompd ( nvcp, hprof, lhrows[0], lhrows,
                                      &qdata->abort ) ) {
if ( !qdata->abort )
  printf ( "%d-", bl );
      qdata->positive = false;
      goto failure;
    }
printf ( "%d+", bl );
    bd->fghflag |= FLAG_LH;
  }
  return true;

failure:
  qdata->abort = true;
  return false;
} /*_g2mbl_MLSDecompSmallBlockd*/

static boolean _g2mbl_MLSDecompCoarseBlockd ( void *usrdata, int3 *jobnum )
{
  g2mbl_blockdecompdata *qdata;
  mesh_ml_optdata       *d;
  short int             bl;
  mlblock_desc          *bd;
  double                nu;
  int                   nwcp;
  int                   hsize, *hprof;
  double                **hrows, **lhrows;

        /* decompose the matrix for a coarse mesh preconditioner */
  qdata = (g2mbl_blockdecompdata*)usrdata;
  d = qdata->d;
  nu = qdata->nu;
  bl = qdata->bl;
  bd = &d->bd[bl];
  if ( !bd->rmnzi )
    return true;    /* nothing to do */

  hsize  = bd->hsize;
  hprof  = bd->hprof;
  lhrows = bd->lhrows;
  nwcp   = bd->nwcp;
  if ( !(bd->fghflag & FLAG_CMH) ) {
    if ( !_g2mbl_CMPSSetupCoarseHessiand ( d, bl, nu ) ) {
printf ( "qq\n" );
      goto failure;
    }
    bd->fghflag |= FLAG_CMH;
    bd->fghflag &= ~FLAG_CMLH;
  }
  if ( (bd->fghflag & FLAG_CMH) && !(bd->fghflag & FLAG_CMLH) ) {
printf ( "(" );
    hrows = bd->hrows;
    memcpy ( lhrows[0], hrows[0], hsize*sizeof(double) );
    if ( pkn_NRBSymCholeskyDecompd ( nwcp, hprof, lhrows[0], lhrows, &qdata->abort ) )
      bd->fghflag |= FLAG_CMLH;
    else {
if ( !qdata->abort )
  printf ( "%d-)", bl );
      bd->fghflag &= ~FLAG_CMH;
      qdata->positive = false;
      goto failure;
    }
printf ( "%d+)", bl );
  }
  return true;

failure:
  qdata->abort = true;
  return false;
} /*_g2mbl_MLSDecompCoarseBlockd*/

boolean _g2mbl_MLSDecomposeBlockPrecond ( mesh_ml_optdata *d, short int bl,
                                          double nu, boolean *positive )
{
  void                  *sp;
  g2mbl_blockdecompdata qdata;
  int                   nsmbl;  /* number of small blocks */
  int                   i;
  mlblock_desc          *bd;
  int3                  jobnum;
  boolean               success;

  sp = pkv_GetScratchMemTop ();
  qdata.d = d;
  qdata.bl = bl;
  qdata.nu = nu;
        /* find out, how many small matrices there are to decompose */
  nsmbl = 0x0002;
  qdata.bl0 = 2*bl+1;
  bd = &d->bd[qdata.bl0];
  while ( bd->iHbl ) {
    nsmbl <<= 1;
    qdata.bl0 = 2*qdata.bl0+1;
    bd = &d->bd[qdata.bl0];
  }
  qdata.bl1 = qdata.bl0 + nsmbl - 1;
  bd = &d->bd[bl];
  *positive = qdata.positive = true;
  jobnum.y = jobnum.z = 1;

  qdata.abort = false;
  if ( _g2mbl_npthreads > 1 ) {
        /* parallel computations */
    jobnum.x = nsmbl;
    if ( !pkv_SetPThreadsToWork ( &jobnum, _g2mbl_npthreads,
                                  4*1048576, 32*1048576,
                                  (void*)&qdata, _g2mbl_MLSDecompSmallBlockd,
                                  (void*)&qdata, _g2mbl_MLSDecompCoarseBlockd,
                                  &success ) )
      goto failure;
    if ( !qdata.positive )
      *positive = false;
  }
  else {
        /* sequential computations */
    for ( i = 0; i < nsmbl; i++ ) {
      jobnum.x = i;
      if ( !_g2mbl_MLSDecompSmallBlockd ( (void*)&qdata, &jobnum ) ) {
        *positive = false;
        break;
      }
    }
    if ( *positive && bd->rmnzi )
      *positive = _g2mbl_MLSDecompCoarseBlockd ( (void*)&qdata, NULL );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSDecomposeBlockPrecond*/

