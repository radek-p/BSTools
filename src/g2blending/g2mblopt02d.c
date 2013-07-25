
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

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

#define RR 4
#define SR 0.6

/* ////////////////////////////////////////////////////////////////////////// */
void g2mbl_OptLMTDeallocated ( void **data )
{
  mesh_lmt_optdata *d;

  if ( *data ) {
    d = (mesh_lmt_optdata*)*data;
    if ( d->nncpi )      PKV_FREE ( d->nncpi );
    if ( d->domelem )    PKV_FREE ( d->domelem );
    if ( d->domelcpind ) PKV_FREE ( d->domelcpind );
    if ( d->hprof )      PKV_FREE ( d->hprof );
    if ( d->ftab )       PKV_FREE ( d->ftab );
    if ( d->Hessian )    PKV_FREE ( d->Hessian );
    if ( d->aqcoeff )    PKV_FREE ( d->aqcoeff );
    if ( d->belind )     PKV_FREE ( d->belind );
    if ( d->BHessian )   PKV_FREE ( d->BHessian );
    if ( d->iHbl )       PKV_FREE ( d->iHbl );
    if ( d->Hbl )        PKV_FREE ( d->Hbl );
    if ( d->bltag )      PKV_FREE ( d->bltag );
    if ( d->mvtag )      PKV_FREE ( d->mvtag );

    if ( d->pmnzi )      PKV_FREE ( d->pmnzi );
    if ( d->wncpi )      PKV_FREE ( d->wncpi );
    if ( d->phprof )     PKV_FREE ( d->phprof );
    if ( d->phrows )     PKV_FREE ( d->phrows );

    PKV_FREE ( *data );
  }
} /*g2mbl_OptLMTDeallocated*/

/* This procedure uses the DFS algorithm to tag the vertices */
/* at the distance 1 and 2 from the mesh boundary vertices.  */
/* These are the vertices, which determine the G2 boundary   */
/* condition for the surface. */
void _g2mbl_TagBoundaryCondVert ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe, int cvn,
                                  char *vtag )
{
  int dist, d, fhe, i, j;

  dist = vtag[cvn]+1;
  if ( dist <= 3 ) {
    d = mv[cvn].degree;
    fhe = mv[cvn].firsthalfedge;
    for ( i = fhe; i < fhe+d; i++ ) {
      j = mhe[mvhei[i]].v1;
      if ( vtag[j] == 0 || vtag[j] > dist ) {
        vtag[j] = dist;
        _g2mbl_TagBoundaryCondVert ( nv, mv, mvhei, nhe, mhe, j, vtag );
      }
    }
  }
} /*_g2mbl_TagBoundaryCondVert*/

/* The two procedures below are used to count the mesh facets making */
/* the domain of the surface */
void _g2mbl_CountRegularFacets ( int d, int *vertnum, int *mtab, void *usrptr )
{
  mesh_lmt_optdata *data;

  data = (mesh_lmt_optdata*)usrptr;
  data->ndomelems ++;
  data->ndomelcpind += 16;
  data->eltypes[1] = true;
} /*_g2mbl_CountRegularFacets*/

void _g2mbl_CountSpecialVertices ( int d, int k, int *vertnum, int *mtab,
                                   void *usrptr )
{
  mesh_lmt_optdata *data;

  data = (mesh_lmt_optdata*)usrptr;
  data->ndomelems ++;
  data->ndomelcpind += 6*k+1;
  data->eltypes[k-3] = true;
} /*_g2mbl_CountSpecialVertices*/

/* The two procedures below store the numbers of vertices (i.e. of basis */
/* functions) relevant (i.e. nonzero) in a regular square or in the */
/* squares around the common special vertex */
void _g2mbl_GetRegularFacetVertNum ( int d, int *vertnum, int *mtab, void *usrptr )
{
  mesh_lmt_optdata *data;
  int              nel, fcpi;
  meshdom_elem     *domelem;

  data = (mesh_lmt_optdata*)usrptr;
  nel = data->ndomelems ++;
  fcpi = data->ndomelcpind;
  data->ndomelcpind += 16;
  domelem = &data->domelem[nel];
  domelem->type = 4;  /* regular */
  domelem->ncp = 16;
  domelem->firstcpi = fcpi;
  domelem->hti = data->hti;
  data->hti += 9*136;
  memcpy ( &data->domelcpind[fcpi], vertnum, 16*sizeof(int) );
  domelem->cfn = mtab[24];
} /*_g2mbl_GetRegularFacetVertNum*/

void _g2mbl_GetSpecialElemVertNum ( int d, int k, int *vertnum, int *mtab,
                                    void *usrptr )
{
  mesh_lmt_optdata *data;
  int              nel, fcpi, ncpi;
  meshdom_elem     *domelem;

  data = (mesh_lmt_optdata*)usrptr;
  nel = data->ndomelems ++;
  ncpi = 6*k+1;
  fcpi = data->ndomelcpind;
  data->ndomelcpind += ncpi;
  domelem = &data->domelem[nel];
  domelem->type = k;    /* 3 or 5,...,16 */
  domelem->ncp = ncpi;
  domelem->firstcpi = fcpi;
  domelem->hti = data->hti;
  data->hti += 9*(ncpi*(ncpi+1))/2;
  memcpy ( &data->domelcpind[fcpi], vertnum, ncpi*sizeof(int) );
  domelem->cfn = mtab[6];
} /*_g2mbl_GetSpecialElemVertNum*/

boolean _g2mbl_AssignMeshd ( mesh_lmt_optdata *d,
                             int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                             int nhe, BSMhalfedge *mhe,
                             int nfac, BSMfacet *mfac, int *mfhei,
                             byte *mkcp )
{
  void *sp;
  int  vi, vb, ei, eb, nvcp, i;
  char *mvtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  if ( !bsm_CheckMeshIntegrity ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei ) )
    goto failure;
        /* attach the mesh to the data structure */
  d->nv    = nv;
  d->mv    = mv;
  d->mvhei = mvhei;
  d->mvcp  = mvcp;
  d->nhe   = nhe;
  d->mhe   = mhe;
  d->nfac  = nfac;
  d->mfac  = mfac;
  d->mfhei = mfhei;
  PKV_MALLOC ( d->mvtag, nv );
  if ( !d->mvtag )
    goto failure;
  mvtag = d->mvtag;

  ftag = pkv_GetScratchMem ( nfac );
  if ( !ftag )
    goto failure;
        /* verify the mesh */
          /* all facets have to be quadrangles */
  for ( i = 0; i < nfac; i++ )
    if ( mfac[i].degree != 4 )
      goto failure;
          /* boundary vertices must be of degree less than 4 */
  bsm_TagMesh ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                mvtag, ftag, &vi, &vb, &ei, &eb );
  for ( i = 0; i < nv; i++ )
    if ( mvtag[i] == 1 ) {
      if ( mv[i].degree >= 4 )
        goto failure;
      else  /* tag the vertices, which determine the boundary condition */
        _g2mbl_TagBoundaryCondVert ( nv, mv, mvhei, nhe, mhe, i, mvtag );
    }
    else if ( mv[i].degree < 3 || mv[i].degree > GH_MAX_K )
      goto failure;

            /* in addition tag the vertices fixed by constraints */
  if ( mkcp ) {
    for ( i = 0; i < nv; i++ )
      if ( mkcp[i] )
        mvtag[i] = 1;
  }

          /* count the vertices to be optimized */
  for ( i = nvcp = 0;  i < nv;  i++ )
    if ( !mvtag[i] )
      nvcp ++;
  if ( !nvcp )  /* nothing to optimize */
    goto failure;
  d->nvcp = nvcp;
  d->nvars = 3*nvcp;
  PKV_MALLOC ( d->nncpi, (nv+nvcp)*sizeof(int) );
  PKV_MALLOC ( d->hprof, 3*nvcp*sizeof(int) );
  if ( !d->nncpi || !d->hprof )
    goto failure;
  d->vncpi = &d->nncpi[nv];
#ifdef _DEBUG
printf ( "nvcp = %d, nvars = %d\n", nvcp, d->nvars );
#endif
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_AssignMeshd*/

void _g2mbl_OutputNZDistr ( mesh_lmt_optdata *d )
{
  void         *sp;
  byte         *nzcdistr;
  int          ndomelems, *domelcpind, *nncpi, ncp, nzcsize;
  meshdom_elem *domelem;
  int          i, j, l, p0, p1, gvni, gvnj;

  if ( !g2mbl_outputnzdistr )
    return;
  sp = pkv_GetScratchMemTop ();
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  nncpi = d->nncpi;
  ncp = d->nvcp;
  nzcsize = pkn_TMBSize ( ncp );
  nzcdistr = pkv_GetScratchMem ( nzcsize );
  if ( nzcdistr ) {
    memset ( nzcdistr, 0, nzcsize );
    for ( i = 0; i < ncp; i++ )
      pkn_TMBElemSet ( nzcdistr, i, i );
    for ( l = 0; l < ndomelems; l++ ) {
      p0 = domelem[l].firstcpi;
      p1 = p0 + domelem[l].ncp;
      for ( i = p0+1; i < p1; i++ ) {
        gvni = nncpi[domelcpind[i]];
        if ( gvni >= 0 && gvni < ncp )
          for ( j = p0; j < i; j++ ) {
            gvnj = nncpi[domelcpind[j]];
            if ( gvnj >= 0 && gvnj < ncp )
              pkn_TMBElemSet ( nzcdistr, gvni, gvnj );
          }
      }
    }
    g2mbl_outputnzdistr ( d->nbl, -1, true, ncp, ncp, nzcdistr );
  }
  pkv_SetScratchMemTop ( sp );
} /*_g2mbl_OutputNZDistr*/

boolean _g2mbl_SetupHessianProfiled ( mesh_lmt_optdata *d, boolean use_blocks )
{
  void         *sp;
  int          nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex    *mv;
  BSMhalfedge  *mhe;
  BSMfacet     *mfac;
  int          nvcp, ndomelems, ndomelcpind, *nncpi, *domelcpind;
  char         *mvtag;
  meshdom_elem *domelem;
  vertex_desc  *nvcpi;
  int          nvars, i, j, k, nspv, vi, vb, ei, eb;
  int          size, nzcdsize, hsize;
  int          *hprof, *vncpi;
  byte         *nzcdistr;

  sp = pkv_GetScratchMemTop ();
  nv = d->nv;
  nhe = d->nhe;
  nfac = d->nfac;
  mv = d->mv;
  mvhei = d->mvhei;
  mhe = d->mhe;
  mfac = d->mfac;
  mfhei = d->mfhei;
  nvcp = d->nvcp;
  nncpi = d->nncpi;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  vncpi = d->vncpi;
  hprof = d->hprof;
  nvars = d->nvars;
  mvtag = d->mvtag;

  nvcpi = pkv_GetScratchMem ( nvcp*sizeof(vertex_desc) );
  if ( !nvcpi )
    goto failure;
  for ( i = j = 0;  i < nv;  i++ )
    if ( mvtag[i] )
      nncpi[i] = -1;
    else {
      nvcpi[j].vnum = vncpi[j] = i;
      nvcpi[j].firstsel = nvcp;
      nvcpi[j].nneigh = 0;
      nncpi[i] = j++;
    }
          /* count the special vertices and their squares */
  for ( i = nspv = 0; i < nv; i++ )
    if ( mvtag[i] != 1 && mv[i].degree != 4 )
      nspv ++;

          /* compute the length of the array domelcpind */
          /* and count the squares */
  d->ndomelems = d->ndomelcpind = 0;
  bsm_FindRegularSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           4, d, _g2mbl_CountRegularFacets );
  bsm_FindSpecialVSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                            2, d, _g2mbl_CountSpecialVertices );
  ndomelems = d->ndomelems;
  ndomelcpind = d->ndomelcpind;
#ifdef _DEBUG
printf ( "ndomelems = %d, ndomelcpind = %d\n", d->ndomelems, d->ndomelcpind );
#endif

  PKV_MALLOC ( d->domelem, ndomelems*sizeof(meshdom_elem) );
  PKV_MALLOC ( d->domelcpind, ndomelcpind*sizeof(int) );
  if ( !d->domelem || !d->domelcpind )
    goto failure;
          /* get the numbers of functions nonzero in the domain elements */
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  d->ndomelems = d->ndomelcpind = d->hti = 0;
  bsm_FindRegularSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           4, d, _g2mbl_GetRegularFacetVertNum );
  bsm_FindSpecialVSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                            2, d, _g2mbl_GetSpecialElemVertNum );
  size = ndomelems + 3*d->ndomelcpind + d->hti;
#ifdef _DEBUG
printf ( "ftabsize = %d\n", size );
#endif
  PKV_MALLOC ( d->ftab, size*sizeof(double) );
  if ( !d->ftab )
    goto failure;
  d->gtab = &d->ftab[ndomelems];
  d->htab = &d->gtab[3*ndomelcpind];
  memset ( d->ftab, 0, size*sizeof(double) );

  if ( !use_blocks || nvars <= CG_THRESHOLD ) {
        /* find the distribution of the Hessian nonzero coefficients */
    nzcdsize = pkn_TMBSize ( nvcp );
    nzcdistr = pkv_GetScratchMem ( nzcdsize );
    if ( !nzcdistr )
      goto failure;
    memset ( nzcdistr, 0, nzcdsize );
    for ( i = 0; i < nvcp; i++ )
      pkn_TMBElemSet ( nzcdistr, i, i );
    for ( i = 0; i < ndomelems; i++ ) {
      vi = domelem[i].ncp;
      vb = domelem[i].firstcpi;
      for ( j = 0; j < vi; j++ ) {
        ei = nncpi[domelcpind[vb+j]];
        if ( ei >= 0 )
          for ( k = 0; k < vi; k++ ) {
            eb = nncpi[domelcpind[vb+k]];
            if ( eb > ei )
              if ( !pkn_TMBTestAndSet ( nzcdistr, eb, ei ) ) {
                nvcpi[eb].nneigh ++;
                nvcpi[ei].nneigh ++;
              }
          }
      }
    }
    if ( g2mbl_outputnzdistr )
      g2mbl_outputnzdistr ( d->nbl, -1, false, nvcp, nvcp, nzcdistr );
    if ( !_g2mbl_OrderCPoints ( nv, nvcp, 0, nzcdsize, nzcdistr,
                                nncpi, vncpi, nvcpi, ndomelems, domelem, domelcpind,
                                NULL, 3, &hsize, hprof, NULL ) )
      goto failure;
    d->hsize = hsize;
#ifdef _DEBUG
printf ( "hsize = %d\n", hsize );
#endif
  }
  if ( use_blocks ) {
    if ( nvars > CG_THRESHOLD ) {
      d->use_cg = true;
      d->hsize = 0;
    }
    else {
      d->use_cg = false;
      size = hsize*sizeof(double) + nvars*sizeof(double*);
      PKV_MALLOC ( d->Hessian, size );
      if ( !d->Hessian )
        goto failure;
      d->hrows = (double**)(&d->Hessian[hsize]);
      pkn_NRBFindRowsd ( nvars, hprof, d->Hessian, d->hrows );
    }
  }
  else {
    d->use_cg = false;
    size = 2*hsize*sizeof(double) + 2*nvars*sizeof(double*);
    PKV_MALLOC ( d->Hessian, size );
    if ( !d->Hessian )
      goto failure;
    d->LHessian = &d->Hessian[hsize];
    d->hrows = (double**)(&d->LHessian[hsize]);
    d->lhrows = &d->hrows[nvars];
    pkn_NRBFindRowsd ( nvars, hprof, d->Hessian, d->hrows );
    pkn_NRBFindRowsd ( nvars, hprof, d->LHessian, d->lhrows );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupHessianProfiled*/

boolean _g2mbl_AllocBFArraysd ( mesh_lmt_optdata *d, int nkn1, int nkn2 )
{
  int    i, j;
  int    size;
  double *aqcoeff, *aqknots, *bqcoeff, *bqknots;
  double *abf, *adbf, *addbf, *adddbf, *bbf, *bdbf, *bddbf, *bdddbf;

        /* compute the size and allocate arrays for basis functions */
  if ( nkn2 != nkn1 ) {
    size = nkn1+nkn2;
    for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
      if ( d->eltypes[j] )
        size += g2mbl_NiSize ( nkn1, i ) + g2mbl_NijSize ( nkn1, i ) +
                g2mbl_MijSize ( nkn1, i ) + g2mbl_NiSize ( nkn2, i ) +
                nkn1*nkn1 + nkn2*nkn2;
  }
  else {
    size = nkn1;
    for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
      if ( d->eltypes[j] )
        size += g2mbl_NiSize ( nkn1, i ) + g2mbl_NijSize ( nkn1, i ) +
                g2mbl_MijSize ( nkn1, i ) + nkn1*nkn1;
  }
#ifdef _DEBUG
printf ( "Nsize = %d\n", size );
#endif
  PKV_MALLOC ( d->aqcoeff, size*sizeof(double) );
  if ( !d->aqcoeff )
    goto failure;
  if ( nkn2 != nkn1 )
    d->bqcoeff = &d->aqcoeff[nkn1];
  else
    d->bqcoeff = d->aqcoeff;
  bqcoeff = &d->bqcoeff[nkn2];
  for ( i = 3, j = 0;  i <= GH_MAX_K;  i++, j++ )
    if ( d->eltypes[j] ) {
      d->aJac[j] = bqcoeff;
      d->aNitabs[j] = d->aJac[j] + nkn1*nkn1;
      d->aNijtabs[j] = d->aNitabs[j] + g2mbl_NiSize ( nkn1, i );
      d->aMijtabs[j] = d->aNijtabs[j] + g2mbl_NijSize ( nkn1, i );
      if ( nkn2 != nkn1 ) {
        d->bJac[j] = d->aMijtabs[j] + g2mbl_MijSize ( nkn1, i );
        d->bNitabs[j] = d->bJac[j] + nkn2*nkn2;
        bqcoeff = d->bNitabs[j] + g2mbl_NiSize ( nkn2, i );
      }
      else {
        d->bJac[j] = d->aJac[j];
        d->bNitabs[j] = d->aNitabs[j];
        bqcoeff = d->aMijtabs[j] + g2mbl_MijSize ( nkn1, i );
      }
    }
        /* evaluate the basis functions for regular squares */
  if ( !_g2bl_TabBasisFuncd ( nkn1, &aqknots, &aqcoeff,
                              &abf, &adbf, &addbf, &adddbf ) )
    goto failure;
  memcpy ( d->aqcoeff, aqcoeff, nkn1*sizeof(double) );
  g2bl_TabNid ( nkn1, abf, adbf, addbf, adddbf, d->aNitabs[1] );
  g2bl_TabNijd ( nkn1, abf, adbf, addbf, d->aNijtabs[1] );
  g2bl_TabMijd ( nkn1, abf, adbf, addbf, adddbf, d->aMijtabs[1] );
  for ( i = 0; i < nkn1*nkn1; i++ )
    d->aJac[1][i] = 1.0;
  if ( nkn2 != nkn1 ) {
    if ( !_g2bl_TabBasisFuncd ( nkn2, &bqknots, &bqcoeff,
                                &bbf, &bdbf, &bddbf, &bdddbf ) )
      goto failure;
    memcpy ( d->bqcoeff, bqcoeff, nkn2*sizeof(double) );
    g2bl_TabNid ( nkn2, bbf, bdbf, bddbf, bdddbf, d->bNitabs[1] );
    for ( i = 0; i < nkn2*nkn2; i++ )
      d->bJac[1][i] = 1.0;
  }
  d->nkn1 = nkn1;
  d->nkn2 = nkn2;
        /* evaluate the basis functions for special elements */
  if ( d->eltypes[0] ) {
    if ( !g2mbl_SetupHolePatchMatrixd ( 3 ) )
      goto failure;
    if ( !g2mbl_TabNid ( 3, nkn1, aqknots, d->aNitabs[0], d->aJac[0], true ) )
      goto failure;
    g2mbl_TabNijd ( 3*6+1, nkn1, d->aNitabs[0], d->aNijtabs[0] );
    g2mbl_TabMijd ( 3*6+1, nkn1, d->aNitabs[0], d->aMijtabs[0] );
    if ( nkn2 != nkn1 ) {
      if ( !g2mbl_TabNid ( 3, nkn2, bqknots, d->bNitabs[0], d->bJac[0], true ) )
        goto failure;;
    }
  }
  for ( i = 5; i <= GH_MAX_K; i++ )
    if ( d->eltypes[i-3] ) {
      if ( !g2mbl_SetupHolePatchMatrixd ( i ) )
        goto failure;
      if ( !g2mbl_TabNid ( i, nkn1, aqknots, d->aNitabs[i-3], d->aJac[i-3], true ) )
        goto failure;
      g2mbl_TabNijd ( i*6+1, nkn1, d->aNitabs[i-3], d->aNijtabs[i-3] );
      g2mbl_TabMijd ( i*6+1, nkn1, d->aNitabs[i-3], d->aMijtabs[i-3] );
      if ( nkn2 != nkn1 ) {
        if ( !g2mbl_TabNid ( i, nkn2, bqknots, d->bNitabs[i-3], d->bJac[i-3], true ) )
          goto failure;
      }
    }
  return true;

failure:
  return false;
} /*_g2mbl_AllocBFArraysd*/

boolean _g2mbl_SetupElemConstd ( mesh_lmt_optdata *d,
                                 double dM, double dO, double C )
{
  void         *sp;
  int          ndomelems, nv;
/*
  int          nhe, nfac;
  BSMvertex    *mv;
  BSMhalfedge  *mhe;
  BSMfacet     *mfac;
  int          *mvhei, *mfhei;
*/
  int          *nncpi;
  point3d      *mvcp;
  meshdom_elem *domelem;
  int          i, j;
  double       r;
/*
  int          *fdist, maxdist;
*/

  sp = pkv_GetScratchMemTop ();
  nv = d->nv;
/*
  nhe = d->nhe;
  nfac = d->nfac;
  mv = d->mv;
  mvhei = d->mvhei;
  mhe = d->mhe;
  mfac = d->mfac;
  mfhei = d->mfhei;
*/
  mvcp = d->mvcp;
  nncpi = d->nncpi;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
        /* set up the constant c for each element */
  if ( dM <= 0.0 ) {
          /* compute the "surface diameter" */
    for ( i = 1; i < nv; i++ )
      if ( nncpi[i] == -1 )
        for ( j = 0; j < i; j++ )
          if ( nncpi[j] == -1 ) {
            r = Point3Distanced ( &mvcp[i], &mvcp[j] );
            dM = max ( dM, r );
          }
    if ( dM <= 0 )
      goto failure;
  }
  dM *= dM;

  if ( dO <= 0.0 ) {
          /* compute the "domain diameter" */
    j = 0;
    for ( i = 0; i < ndomelems; i++ )
      if ( domelem[i].type == 4 )
        j ++;
      else
        j += domelem[i].type;
    dO = (double)j;
  }
  else
    dO *= dO;
  r = C*dO*dO/(dM*dM);
/*
  fdist = pkv_GetScratchMemi ( nfac );
  if ( !fdist )
    goto failure;
  if ( !g2mbl_FindDistances ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                              fdist ) )
    goto failure;
  maxdist = 0;
  for ( i = 0; i < ndomelems; i++ )
    if ( domelem[i].type == 4 )
      maxdist = max ( maxdist, fdist[domelem[i].cfn] );
*/
  for ( i = 0; i < ndomelems; i++ ) {
/*
    if ( domelem[i].type != 4 )
      domelem[i].C = r * pkv_rpower ( SR, RR );
    else if ( fdist[domelem[i].cfn] < RR )
      domelem[i].C = r * pkv_rpower ( SR, RR-fdist[domelem[i].cfn] );
    else
*/
      domelem[i].C = r;
  }
#ifdef _DEBUG
/*  printf ( "maxdist = %d\n", maxdist ); */
  printf ( "C = %8.5f\n", r );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupElemConstd*/

