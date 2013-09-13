
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

#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
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

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLSSetupBlockCGHessiand ( mesh_ml_optdata *d, int bn )
{
  void          *sp;
  int           Hblsize, i, j, fnz, nnz, fnz1, nnz1;
  mlblock_desc  *bd, *bd1;
  int           nbcp1, *vncpi1, *ind1, *cHbl1, *tHbl1;
  int           nbcp, *vncpi, *ind, *cHbl, *tHbl;
  nzHbl_rowdesc *iHbl, *iHbl1;

  sp = pkv_GetScratchMemTop ();
  bd = &d->bd[bn];
  if ( bn == 0 ) {
    if ( !_g2mbl_MLFindVCPNeighboursd ( d, &bd->iHbl, &bd->cHbl ) ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_18, ERRMSG_18 );
      goto failure;
    }
    Hblsize = d->Hblsize;
    PKV_MALLOC ( d->Hbl, Hblsize*sizeof(double) );
    if ( !d->Hbl ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
      goto failure;
    }
    bd->nHbl = Hblsize;
    tHbl = bd->tHbl = &bd->cHbl[Hblsize];
    for ( i = 0; i < Hblsize; i++ )
      tHbl[i] = i;
  }
  else {
    bd1 = &d->bd[(bn-1)/2];  /* number of the block 1 up in the hierarchy */
    nbcp   = bd->nvcp;
    vncpi  = bd->vncpi;
    nbcp1  = bd1->nvcp;
    vncpi1 = bd1->vncpi;
    iHbl1  = bd1->iHbl;
    cHbl1  = bd1->cHbl;
    tHbl1  = bd1->tHbl;
    ind  = pkv_GetScratchMemi ( nbcp );
    ind1 = pkv_GetScratchMemi ( nbcp1 );
    if ( !ind || !ind1 ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
        /* vncpi1 contains an increasing sequence of integers (indices of */
        /* control points of the block above) and vncpi contains a subsequence */
        /* the first task is to find the indices cross-referencing the arrays */
    for ( i = j = 0;  i < nbcp1;  i++ )
      if ( vncpi[j] == vncpi1[i] ) {
        ind[j] = i;
        ind1[i] = j ++;
      }
      else
        ind1[i] = -1;
        /* after that the job is to select the rows and vertices of the Hessian */
        /* first, count the nonzero coefficients in the selected rows and columns */
      Hblsize = 0;
      for ( i = 0; i < nbcp1; i++ )
        if ( ind1[i] >= 0 ) {
          fnz1 = iHbl1[i].firsthbl;
          nnz1 = iHbl1[i].nhbl;
          for ( j = 0; j < nnz1; j++ )
            if ( ind1[cHbl1[fnz1+j]] >= 0 )
              Hblsize ++;
        }
    PKV_MALLOC ( bd->iHbl, nbcp*sizeof(nzHbl_rowdesc) );
    PKV_MALLOC ( bd->cHbl, 2*Hblsize*sizeof(int) );
    if ( !bd->iHbl || !bd->cHbl ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
      goto failure;
    }
    iHbl = bd->iHbl;
    cHbl = bd->cHbl;
    tHbl = bd->tHbl = &cHbl[Hblsize];
        /* now select the coefficients, i.e. store their indices in the arrays */
    for ( i = fnz = 0;  i < nbcp1;  i++ )
      if ( ind1[i] >= 0 ) {
        fnz1 = iHbl1[i].firsthbl;
        nnz1 = iHbl1[i].nhbl;
        for ( j = nnz = 0;  j < nnz1;  j++ )
          if ( ind1[cHbl1[fnz1+j]] >= 0 ) {
            cHbl[fnz+nnz] = ind1[cHbl1[fnz1+j]];
            tHbl[fnz+nnz] = tHbl1[fnz1+j];
            nnz ++;
          }
        iHbl[ind1[i]].firsthbl = fnz;
        iHbl[ind1[i]].nhbl = nnz;
        fnz += nnz;
      }
  }
  bd->nHbl = Hblsize;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSSetupBlockCGHessiand*/

boolean _g2mbl_MLSSetupBlockCholHessiand ( mesh_ml_optdata *d, int bn )
{
  void         *sp;
  meshdom_elem *domelem;
  mlblock_desc *bd;
  int          nv, nvcp, *vncpi, ndomelems,
               *domelcpind, ndomel, *domelind, *nncpi;
  int          hsize;
  int          nzcdsize;
  byte         *nzcdistr;
  int          i, j, k, vj, vk, fcp, ncp, el;
  vertex_desc  *vertd;

  sp = pkv_GetScratchMemTop ();
  ndomelems  = d->ndomelems;   
  domelem    = d->domelem;     
  domelcpind = d->domelcpind;  
  nv         = d->nv;
  bd         = &d->bd[bn];
  nvcp       = bd->nvcp;  
  vncpi      = bd->vncpi;
  ndomel     = bd->ndomel;
  domelind   = bd->domelind;
  nzcdsize   = pkn_TMBSize ( nvcp );
  nzcdistr   = pkv_GetScratchMem ( nzcdsize );
  nncpi      = pkv_GetScratchMemi ( nv );
  vertd      = pkv_GetScratchMem ( nvcp*sizeof(vertex_desc) );
  if ( !nzcdistr || !nncpi || !vertd ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
        /* find the number of variable vertices in the block */
        /* for each vertex, -1 for the vertices not from the block */
  memset ( nncpi, 0xFF, nv*sizeof(int) );
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* find the distribution of the nonzero Hessian coefficients */
        /* and the neighbours for each variable vertex */
  for ( i = 0; i < nvcp; i++ ) {
    vertd[i].vnum = vncpi[i];
    vertd[i].firstsel = nvcp;
    vertd[i].nneigh = 0;
  }
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < nvcp; i++ )
    pkn_TMBElemSet ( nzcdistr, i, i );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    fcp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      vj = nncpi[domelcpind[fcp+j]];
      if ( vj >= 0 ) {
        for ( k = 0; k < ncp; k++ ) {
          vk = nncpi[domelcpind[fcp+k]];
          if ( vk > vj )
            if ( !pkn_TMBTestAndSet ( nzcdistr, vj, vk ) ) {
              vertd[vj].nneigh ++;
              vertd[vk].nneigh ++;
            }
        }
      }
    }
  }
  if ( g2mbl_outputnzdistr )
    g2mbl_outputnzdistr ( d->nblocks, bn, false, d->nvcp, nvcp, nzcdistr );

  PKV_MALLOC ( bd->hprof, nvcp*sizeof(int) );
  if ( !bd->hprof ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;
  }
  if ( !_g2mbl_OrderCPoints ( nv, nvcp, 0, nzcdsize, nzcdistr, nncpi, vncpi,
                              vertd, ndomelems, domelem, domelcpind,
                              NULL, 1, &bd->hsize, bd->hprof, NULL ) ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_19, ERRMSG_19 );
    goto failure;
  }
  hsize = bd->hsize;

  if ( g2mbl_outputnzdistr ) {
    memset ( nzcdistr, 0, nzcdsize );
    for ( i = 0; i < nvcp; i++ )
      pkn_TMBElemSet ( nzcdistr, i, i );
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      fcp = domelem[el].firstcpi;
      ncp = domelem[el].ncp;
      for ( j = 0; j < ncp; j++ ) {
        vj = nncpi[domelcpind[fcp+j]];
        if ( vj >= 0 )
          for ( k = 0; k < ncp; k++ ) {
            vk = nncpi[domelcpind[fcp+k]];
            if ( vk > vj )
              pkn_TMBElemSet ( nzcdistr, vk, vj );
          }
      }
    }
    g2mbl_outputnzdistr ( d->nblocks, bn, true, d->nvcp, nvcp, nzcdistr );
  }
  PKV_MALLOC ( bd->hrows, 2*(nvcp*sizeof(double*)+hsize*sizeof(double)) );
  if ( !bd->hrows ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_9, ERRMSG_9 );
    goto failure;  
  }
  bd->lhrows = &bd->hrows[nvcp];
  pkn_NRBFindRowsd ( nvcp, bd->hprof, (double*)&bd->lhrows[nvcp], bd->hrows );
  pkn_NRBFindRowsd ( nvcp, bd->hprof, &bd->hrows[0][hsize], bd->lhrows );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSSetupBlockCholHessiand*/

boolean _g2mbl_MLSSetupBlockHessiansd ( mesh_ml_optdata *d )
{
  int          i, nblocks;
  mlblock_desc *bd;

  nblocks = d->nblocks;
  for ( i = 0; 2*i+1 < nblocks; i++ ) {
    bd = &d->bd[i];
    if ( bd->nvcp <= MAX_NVCPS ) {
      if ( !_g2mbl_MLSSetupBlockCholHessiand ( d, i ) ) {
        PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_20, ERRMSG_20 );
        return false;
      }
    }
    else {
      if ( !_g2mbl_MLSSetupBlockCGHessiand (d, i ) ) {
        PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_21, ERRMSG_21 );
        return false;
      }
    }
  }
  for ( ; i < nblocks; i++ )
    if ( !_g2mbl_MLSSetupBlockCholHessiand ( d, i ) ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_20, ERRMSG_20 );
      return false;
    }
  return true;
} /*_g2mbl_MLSSetupBlockHessiansd*/

