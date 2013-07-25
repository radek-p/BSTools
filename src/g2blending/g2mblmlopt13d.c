
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
#define _CONST
#include "g2mblmlprivated.h"
#include "msgpool.h"

#define _DEBUG

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_CMPMultRTHR3x3d ( int nrowsa, int nnza, index3 *ai, double *ac,
                                 int ncolsb, int nnzb, index3 *bi, double *bc,
                                 int nnz1,
                                 double nu,
                                 int hsize, int *hprof, double **hrows )
{
  void         *sp, *sp1;
  unsigned int *apermut, *bpermut, *abpermut;
  int          *colsa, *colsab, *colsb, *rowsb;
  index2       *abi;
  double       *abnzc, s[9], t[9], b;
  int          i, j, k, k0, k1, l, l0, l1, m, n, cnt, nnz, ii, jj, kk;
  int          wsps, wspn, *wsp;
  unsigned int *auxpermut;
  index2       *aikbkj, hind;

  sp = pkv_GetScratchMemTop ();
  apermut = abpermut = NULL;
  abi = NULL;
  abnzc = NULL;
        /* multiplication M = A*B */
  PKV_MALLOC ( apermut, (nnza+nnzb+nrowsa+ncolsb+2)*sizeof(int) );
  if ( !apermut )
    goto failure;
  colsa = (int*)&apermut[nnza];
  bpermut = (unsigned int*)&colsa[nrowsa+1];
  colsb = (int*)&bpermut[nnzb];
  PKV_MALLOC ( abi, nnz1*sizeof(index2) );
  PKV_MALLOC ( abnzc, 9*nnz1*sizeof(double) );
  if ( !abi || !abnzc )
    goto failure;
  if ( !pkn_SPsubMFindCols ( nrowsa, nrowsa, nnza, ai, apermut, false, colsa ) )
    goto failure;
  if ( !pkn_SPsubMFindCols ( nrowsa, ncolsb, nnzb, bi, bpermut, false, colsb ) )
    goto failure;
  sp1 = pkv_GetScratchMemTop ();
  wsp = NULL;
  aikbkj = NULL;
  wsps = 0;
  nnz = 0;
  for ( i = 0; i < ncolsb; i++ ) {
    k0 = colsb[i];
    k1 = colsb[i+1];
        /* find the workspace size and allocate the workspace */
    wspn = 0;
    for ( k = k0; k < k1; k++ ) {
      j = bi[bpermut[k]].i;
      wspn += colsa[j+1]-colsa[j];
    }
    if ( wspn > wsps ) {
      pkv_SetScratchMemTop ( sp1 );
      wsp = pkv_GetScratchMemi ( 4*wspn );
      if ( !wsp )
        goto failure;
      aikbkj = (index2*)&wsp[wspn];
      wsps = wspn;
    }
        /* find the positions (row numbers) of nonzero product coefficients */
    cnt = 0;
    for ( k = k0; k < k1; k++ ) {
      m = bpermut[k];
      j = bi[m].i;
      l0 = colsa[j];
      l1 = colsa[j+1];
      for ( l = l0; l < l1; l++ ) {
        n = apermut[l];
        aikbkj[cnt].i = (ai[n].i >= ai[n].j) ? ai[n].k : -1-ai[n].k;
        aikbkj[cnt].j = bi[m].k;
        wsp[cnt++] = ai[n].i;
      }
    }
    if ( cnt ) {
          /* set up the identity permutation */
      auxpermut = (unsigned int*)&aikbkj[wspn];
      for ( l = 0; l < cnt; l++ )
        auxpermut[l] = l;
          /* sort by columns */
      if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(int),
                            0, cnt, wsp, auxpermut ) != SORT_OK )
        goto failure;
      k = auxpermut[0];
      l0 = wsp[k];
      abi[nnz].j = i;
      abi[nnz].i = l0;
      b = bc[aikbkj[k].j];
      if ( aikbkj[k].i >= 0 )
        pkn_MultMatrixNumd ( 1, 9, 0, &ac[aikbkj[k].i], b, 0, s );
      else {
        pkv_TransposeMatrixd ( 3, 3, 3, &ac[-1-aikbkj[k].i], 3, t );
        pkn_MultMatrixNumd ( 1, 9, 0, t, b, 0, s );
      }
      for ( l = 1; l < cnt; l++ ) {
        k = auxpermut[l];
        b = bc[aikbkj[k].j];
        if ( wsp[k] != l0 ) {
          memcpy ( &abnzc[9*nnz], s, 9*sizeof(double) );
          nnz ++;
          if ( aikbkj[k].i >= 0 )
            pkn_MultMatrixNumd ( 1, 9, 0, &ac[aikbkj[k].i], b, 0, s );
          else {
            pkv_TransposeMatrixd ( 3, 3, 3, &ac[-1-aikbkj[k].i], 3, t );
            pkn_MultMatrixNumd ( 1, 9, 0, t, b, 0, s );
          }
          l0 = wsp[k];
          abi[nnz].j = i;
          abi[nnz].i = l0;
        }
        else {
          if ( aikbkj[k].i >= 0 )
            pkn_AddMatrixMd ( 1, 9, 0, s, 0, &ac[aikbkj[k].i], b, 0, s );
          else {
            pkv_TransposeMatrixd ( 3, 3, 3, &ac[-1-aikbkj[k].i], 3, t );
            pkn_AddMatrixMd ( 1, 9, 0, s, 0, t, b, 0, s );
          }
        }
      }
      memcpy ( &abnzc[9*nnz], s, 9*sizeof(double) );
      nnz ++;
    }
  }
  pkv_SetScratchMemTop ( sp1 );

  if ( nu != 0.0 ) {
        /* add nu*B to the product A*B; both representations are */
        /* ordered columnwise, B via the bpermut array, A*B as is */
    for ( i = j = 0; i < nnzb; i++ ) {
      k = bpermut[i];
      while ( abi[j].i != bi[k].i || abi[j].j != bi[k].j )
        j++;
      l = 9*j;
      b = nu*bc[bi[k].k];
      abnzc[l]   += b;
      abnzc[l+4] += b;
      abnzc[l+8] += b;
    }
  }

        /* multiplication B^T*M */
  PKV_MALLOC ( abpermut, (nnz1+nrowsa+ncolsb+2)*sizeof(int) );
  if ( !abpermut )
    goto failure;
  colsab = (int*)&abpermut[nnz1];
  rowsb = &colsab[ncolsb+1];
  if ( !pkn_SPsubMFindRows ( nrowsa, ncolsb, nnzb, bi, bpermut, false, rowsb ) )
    goto failure;
  if ( !pkn_SPMFindCols ( nrowsa, ncolsb, nnz1, abi, abpermut, false, colsab ) )
    goto failure;
  
  memset ( hrows[0], 0, hsize*sizeof(double) );
  sp1 = pkv_GetScratchMemTop ();
  wsps = 0;
  for ( i = 0; i < ncolsb; i++ ) {
    k0 = colsab[i];
    k1 = colsab[i+1];
    wspn = 0;
    for ( k = k0; k < k1; k++ ) {
      j = abi[abpermut[k]].i;
      wspn += rowsb[j+1]-rowsb[j];
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
      m = abpermut[k];
      j = abi[m].i;
      l0 = rowsb[j];
      l1 = rowsb[j+1];
      for ( l = l0; l < l1; l++ ) {
        n = bpermut[l];
        if ( bi[n].j >= i ) {
          aikbkj[cnt].i = bi[n].k;
          aikbkj[cnt].j = 9*m;
          wsp[cnt++] = bi[n].j;
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
      pkn_MultMatrixNumd ( 1, 9, 0, &abnzc[aikbkj[k].j], bc[aikbkj[k].i], 0, s );
      for ( l = 1; l < cnt; l++ ) {
        k = auxpermut[l];
        if ( wsp[k] != l0 ) {
#define DEBUG
#ifdef DEBUG
if ( hind.i < hind.j || 3*hind.j < hprof[3*hind.i] ) {
  printf ( "blad profilu!" );
  exit ( 1 );
}
#endif
          if ( hind.i > hind.j ) {
            for ( ii = kk = 0;  ii < 3;  ii++ )
              for ( jj = 0;  jj < 3;  jj++, kk++ )
                hrows[3*hind.i+ii][3*hind.j+jj] = s[kk];
          }
          else {
            for ( ii = 0;  ii < 3;  ii++ )
              for ( jj = 0;  jj <= ii;  jj++ )
                hrows[3*hind.i+ii][3*hind.j+jj] = s[3*ii+jj];
          }
          pkn_MultMatrixNumd ( 1, 9, 0, &abnzc[aikbkj[k].j], bc[aikbkj[k].i], 0, s );
          hind.i = l0 = wsp[k];
        }
        else
          pkn_AddMatrixMd ( 1, 9, 0, s, 0, &abnzc[aikbkj[k].j],
                            bc[aikbkj[k].i], 0, s );
      }
#ifdef DEBUG
if ( hind.i < hind.j || 3*hind.j < hprof[3*hind.i] ) {
  printf ( "blad profilu!" );
  exit ( 1 );
}
#endif
      if ( hind.i > hind.j ) {
        for ( ii = kk = 0;  ii < 3;  ii++ )
          for ( jj = 0;  jj < 3;  jj++, kk++ )
            hrows[3*hind.i+ii][3*hind.j+jj] = s[kk];
      }
      else {
        for ( ii = 0;  ii < 3;  ii++ )
          for ( jj = 0;  jj <= ii;  jj++ )
            hrows[3*hind.i+ii][3*hind.j+jj] = s[3*ii+jj];
      }
    }
  }

  PKV_FREE ( apermut );
  PKV_FREE ( abi );
  PKV_FREE ( abnzc );
  PKV_FREE ( abpermut );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( apermut )  PKV_FREE ( apermut );
  if ( abi )      PKV_FREE ( abi );
  if ( abnzc )    PKV_FREE ( abnzc );
  if ( abpermut ) PKV_FREE ( abpermut );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPMultRTHR3x3d*/

boolean _g2mbl_CMPSetupCoarseHessiand ( mesh_ml_optdata *d, int bl, double nu )
{
  void          *sp;
  mlblock_desc  *bd;
  double        *Hbl, *rmnzc;
  int           nHbl, *cHbl, *tHbl;
  nzHbl_rowdesc *iHbl;
  int           brmnnz, bHnnz;
  index3        *brmnzi, *Hblnzi;
  int           *hprof, hsize, nvcp, nwcp;
  double        **hrows;
  int           i, j, k, l, j0, j1;

  sp = pkv_GetScratchMemTop ();
  bd = &d->bd[bl];
  Hblnzi = NULL;
  if ( bd->rmnzi ) {
    Hbl   = d->Hbl;
    rmnzc = d->rmnzc;
        /* get block description data */
    nHbl   = bd->nHbl;
    iHbl   = bd->iHbl;
    cHbl   = bd->cHbl;
    tHbl   = bd->tHbl;
    nvcp   = bd->nvcp;
    nwcp   = bd->nwcp;
    brmnnz = bd->rmnnz;
    brmnzi = bd->rmnzi;
    hprof  = bd->hprof;
    hsize  = bd->hsize;
    hrows  = bd->hrows;
        /* find the sparse representation of the Hessian matrix */
    bHnnz = 2*nHbl-nvcp;
    PKV_MALLOC ( Hblnzi, bHnnz*sizeof(index3) );
    if ( !Hblnzi ) {
printf ( "QQ a\n" );
      goto failure;
    }
    for ( i = l = 0;  i < nvcp;  i++ ) {
      j0 = iHbl[i].firsthbl;
      j1 = j0 + iHbl[i].nhbl - 1;
      for ( j = j0; j < j1; j++ ) {
        k = cHbl[j];
        Hblnzi[l  ].i = i;
        Hblnzi[l  ].j = k;
        Hblnzi[l++].k = tHbl[j];
        Hblnzi[l  ].i = k;
        Hblnzi[l  ].j = i;
        Hblnzi[l++].k = tHbl[j];
      }
      Hblnzi[l].i = Hblnzi[l].j = i;
      Hblnzi[l++].k = tHbl[j1];
    }
        /* multiply the block Hessian by the proper refinement submatrix */
    if ( !_g2mbl_CMPMultRTHR3x3d ( nvcp, bHnnz, Hblnzi, Hbl,
                                   nwcp, brmnnz, brmnzi, rmnzc, bd->nnz1,
                                   nu,
                                   hsize, hprof, hrows ) ) {
printf ( "QQ b\n" );
      goto failure;
    }
  }
  pkv_SetScratchMemTop ( sp );
  PKV_FREE ( Hblnzi );
  return true;

failure:
  if ( Hblnzi ) PKV_FREE ( Hblnzi );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPSetupCoarseHessiand*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLmultAxd ( int nvars, void *usrdata,
                           _CONST double *x, double *Ax )
{
  mesh_ml_cg_data *cgd;
  mesh_ml_optdata *d;
  mlblock_desc    *bd;
  int             bl, nvcp;
  nzHbl_rowdesc   *iHbl;
  int             *cHbl, *tHbl, i, j, j0, j1;
  double          *Hbl, *bH, *bx1, *bAx1, *bx2, *bAx2, nu;
  
  cgd = (mesh_ml_cg_data*)usrdata;
  d = cgd->d;
  Hbl = d->Hbl;
  bl  = cgd->bl;
  bd = &d->bd[bl];
  nvcp  = bd->nvcp;
  iHbl  = bd->iHbl;
  cHbl  = bd->cHbl;
  tHbl  = bd->tHbl;
  memset ( Ax, 0, nvars*sizeof(double) );
  for ( i = 0; i < nvcp; i++ ) {
    bAx1 = &Ax[3*i];
    bx2 = &x[3*i];
    j0 = iHbl[i].firsthbl;
    j1 = j0 + iHbl[i].nhbl-1;
    for ( j = j0; j < j1; j++ ) {
      bH = &Hbl[tHbl[j]];
      bAx2 = &Ax[3*cHbl[j]];
      bx1 = &x[3*cHbl[j]];
        /* block below the diagonal */
      bAx1[0] += bH[0]*bx1[0] + bH[1]*bx1[1] + bH[2]*bx1[2];
      bAx1[1] += bH[3]*bx1[0] + bH[4]*bx1[1] + bH[5]*bx1[2];
      bAx1[2] += bH[6]*bx1[0] + bH[7]*bx1[1] + bH[8]*bx1[2];
        /* block above the diagonal */
      bAx2[0] += bH[0]*bx2[0] + bH[3]*bx2[1] + bH[6]*bx2[2];
      bAx2[1] += bH[1]*bx2[0] + bH[4]*bx2[1] + bH[7]*bx2[2];
      bAx2[2] += bH[2]*bx2[0] + bH[5]*bx2[1] + bH[8]*bx2[2];
    }
        /* diagonal block */
    bH = &Hbl[tHbl[j1]];
    bAx1[0] += bH[0]*bx2[0] + bH[1]*bx2[1] + bH[2]*bx2[2];
    bAx1[1] += bH[3]*bx2[0] + bH[4]*bx2[1] + bH[5]*bx2[2];
    bAx1[2] += bH[6]*bx2[0] + bH[7]*bx2[1] + bH[8]*bx2[2];
  }
  if ( (nu = cgd->nu) > 0.0 )
    for ( i = 0; i < nvars; i++ )
      Ax[i] += nu*x[i];
  return true;
} /*_g2mbl_MLmultAxd*/

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
    mesh_ml_cg_data *cgd;
    int             nvars, *nncpi;
    double          *x, *Qix;
    short int       bl, bl0, bl1;
    pthread_mutex_t mutex;
  } g2mbl_multqixdata;

static boolean _g2mbl_MLmultQIxd4 ( void *usrdata, int3 *jobnum )
{
  boolean           *sp;
  g2mbl_multqixdata *qdata;
  mesh_ml_optdata   *d;
  mesh_ml_cg_data   *cgd;
  short int         bl;
  mlblock_desc      *bd;
  int               nvcp, *nncpi, *vncpi, *avncpi;
  int               i, j, k, nvars, *hprof;
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
  nvars = 3*nvcp;
  ax = pkv_GetScratchMemd ( nvars );
  if ( !ax )
    goto failure;
  x = qdata->x;
  for ( i = j = 0;  i < nvcp;  i++, j += 3 )
    memcpy ( &ax[j], &x[3*avncpi[i]], sizeof(point3d) );
  hprof = bd->hprof;
  lhrows = bd->lhrows;
  pkn_NRBLowerTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
  pkn_NRBUpperTrSolved ( nvars, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
        /* add the multiplication result */
  x = qdata->Qix;
  if ( _g2mbl_npthreads > 1 ) {
          /* with pthreads in use, mutual exclusion is needed  */
    pthread_mutex_lock ( &qdata->mutex );
    for ( i = j = 0;  i < nvcp;  i++, j += 3 ) {
      k = 3*avncpi[i];
      AddVector3d ( (vector3d*)&x[k], (vector3d*)&ax[j], (vector3d*)&x[k] );
    }
    pthread_mutex_unlock ( &qdata->mutex );
  }
  else {
    for ( i = j = 0;  i < nvcp;  i++, j += 3 ) {
      k = 3*avncpi[i];
      AddVector3d ( (vector3d*)&x[k], (vector3d*)&ax[j], (vector3d*)&x[k] );
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLmultQIxd4*/

static boolean _g2mbl_MLmultQIxd5 ( void *usrdata, int3 *jobnum )
{
  boolean           *sp;
  g2mbl_multqixdata *qdata;
  mesh_ml_optdata   *d;
  mesh_ml_cg_data   *cgd;
  short int         bl;
  mlblock_desc      *bd;
  int               *hprof, nvcp, nvars, nwcp, nvars3;
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
    nvars  = 3*nvcp;
    nwcp   = bd->nwcp;
    nvars3 = 3*nwcp;
    rmnnz = bd->rmnnz;
    rmnzi = bd->rmnzi;
    rmnzc = d->rmnzc;
    ax = pkv_GetScratchMemd ( nvars3+nvars );
    if ( !ax )
      goto failure;
    ay = &ax[nvars3];
    x = qdata->x;
    if ( !pkn_MultSPsubMTVectord ( nvcp, nwcp, rmnnz, rmnzi, rmnzc,
                                   3, x, ax ) )
      goto failure;
    pkn_NRBLowerTrSolved ( nvars3, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
    pkn_NRBUpperTrSolved ( nvars3, hprof, lhrows[0], lhrows, 1, 1, ax, 1, ax );
    if ( !pkn_MultSPsubMVectord ( nvcp, nwcp, rmnnz, rmnzi, rmnzc,
                                  3, ax, ay ) )
      goto failure;
    x = qdata->Qix;
    if ( _g2mbl_npthreads > 1 ) {
      pthread_mutex_lock ( &qdata->mutex );
      pkn_AddMatrixd ( 1, nvars, 0, x, 0, ay, 0, x );
      pthread_mutex_unlock ( &qdata->mutex );
    }
    else
      pkn_AddMatrixd ( 1, nvars, 0, x, 0, ay, 0, x );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLmultQIxd5*/

boolean _g2mbl_MLmultQIxd ( int nvars, void *usrdata,
                            _CONST double *x, double *Qix )
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
                                  (void*)&qdata, _g2mbl_MLmultQIxd4,
                                  (void*)&qdata, _g2mbl_MLmultQIxd5,
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
      if ( !_g2mbl_MLmultQIxd4 ( (void*)&qdata, &jobsize ) ) {
        result = false;
        break;
      }
    }
    bd = &d->bd[bl];
    if ( result && bd->rmnzi )
      result = _g2mbl_MLmultQIxd5 ( (void*)&qdata, NULL );
  }
  pkv_SetScratchMemTop ( sp );
  return result;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLmultQIxd*/

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
    mesh_ml_optdata *d;
    short int       bl, bl0, bl1;  /* main, first small and last small block */
    double          nu;
    boolean         positive, abort;
  } g2mbl_blockdecompdata;

static boolean _g2mbl_MLDecompSmallBlockd ( void *usrdata, int3 *jobnum )
{
  g2mbl_blockdecompdata *qdata;
  mesh_ml_optdata       *d;
  short int             bl;
  mlblock_desc          *bd0, *bd;
  double                nu;
  int                   nv, nvcp, nvars, *vncpi, i;
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
  nvars = bd->nvars;
  vncpi = bd->vncpi;
  hsize = bd->hsize;
  hprof = bd->hprof;
  hrows = bd->hrows;
  lhrows = bd->lhrows;
  if ( !(bd->fghflag & FLAG_H) ) {  /* invalid Hessian, recompute */
    if ( !g2mbl_MLGetHessianRowsd ( nv, bd0->nvcp, bd0->vncpi,
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
      for ( i = 0; i < nvars; i++ )
        lhrows[i][i] += nu;
    if ( !pkn_NRBSymCholeskyDecompd ( nvars, hprof, lhrows[0], lhrows,
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
} /*_g2mbl_MLDecompSmallBlockd*/

static boolean _g2mbl_MLDecompCoarseBlockd ( void *usrdata, int3 *jobnum )
{
  g2mbl_blockdecompdata *qdata;
  mesh_ml_optdata       *d;
  short int             bl;
  mlblock_desc          *bd;
  double                nu;
  int                   nwcp, nvars3;
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
  nvars3 = 3*nwcp;
  if ( !(bd->fghflag & FLAG_CMH) ) {
    if ( !_g2mbl_CMPSetupCoarseHessiand ( d, bl, nu ) ) {
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
    if ( pkn_NRBSymCholeskyDecompd ( nvars3, hprof, lhrows[0], lhrows,
                                     &qdata->abort ) )
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
} /*_g2mbl_MLDecompCoarseBlockd*/

boolean _g2mbl_MLDecomposeBlockPrecond ( mesh_ml_optdata *d, short int bl,
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
                                  4*1048576, 16*1048576,
                                  (void*)&qdata, _g2mbl_MLDecompSmallBlockd,
                                  (void*)&qdata, _g2mbl_MLDecompCoarseBlockd,
                                  &success ) )
      goto failure;
    if ( !qdata.positive )
      *positive = false;
  }
  else {
        /* sequential computations */
    for ( i = 0; i < nsmbl; i++ ) {
      jobnum.x = i;
      if ( !_g2mbl_MLDecompSmallBlockd ( (void*)&qdata, &jobnum ) ) {
        *positive = false;
        break;
      }
    }
    if ( *positive && bd->rmnzi )
      *positive = _g2mbl_MLDecompCoarseBlockd ( (void*)&qdata, NULL );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLDecomposeBlockPrecond*/

