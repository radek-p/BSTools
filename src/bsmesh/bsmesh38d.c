
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file has been written by Anna Sierhej                                */
/* and modified by Przemyslaw Kiciak                                         */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"
#include "bsfile.h"

/*#define DEBUG*/
#define DEBUG1
/* ////////////////////////////////////////////////////////////////////////// */
boolean bsm_SimplifyMeshd ( int spdimen, int inv,
                            const BSMvertex *imv, const int *imvhei,
                            double *iptc, int inhe, const BSMhalfedge *imhe,
                            int infac, const BSMfacet *imfac, const int *imfhei,
                            int *nboxes,
                            int *onv, BSMvertex **omv, int **omvhei,
                            double **optc, int *onhe, BSMhalfedge **omhe,
                            int *onfac, BSMfacet **omfac, int **omfhei )
{
  double      *xyz, *minxyz, *maxxyz, *eps;
  int         i, j, k, l, vmed, vnext, he, fhe, deg, v, ile;
  int         *nxyz, ind, err;
  int         *Boxes, totalboxes;
  pkv_queue   **Queues;
  pkv_queue   *Bufor;
  int         *newvnum;
  void        *sp;
  boolean     *disabled;
  int         tnv, tnhe, tnfac;
  BSMvertex   *tmv;
  int         *tmvhei, *tmfhei;
  double      *tptc;
  BSMhalfedge *tmhe;
  BSMfacet    *tmfac;
  int         t2nv, t2nhe, t2nfac;
  BSMvertex   *t2mv;
  int         *t2mvhei, *t2mfhei;
  double      *t2ptc;
  BSMhalfedge *t2mhe;
  BSMfacet    *t2mfac;
  boolean     isT;
#ifdef DEBUG
  char        filename[20];
#endif

  sp = pkv_GetScratchMemTop ();
  *omv = NULL;  *omvhei = *omfhei = NULL;  *optc = NULL;
  *omhe = NULL;  *omfac = NULL;

  Bufor = NULL;
  Queues = NULL;
  xyz = pkv_GetScratchMemd ( 4*spdimen );
  nxyz = pkv_GetScratchMemi ( spdimen );
  newvnum = pkv_GetScratchMemi ( inv );
  disabled = pkv_GetScratchMem ( inv*sizeof(boolean) );
  if ( !xyz || !nxyz || !newvnum || !disabled )
    goto failure;
  minxyz = &xyz[spdimen];
  maxxyz = &minxyz[spdimen];
  eps = &maxxyz[spdimen];

    /*zarezerwowanie miejsca na tablice na t*/
  tmfac = pkv_GetScratchMem ( infac*sizeof(BSMfacet) );
  tmfhei = pkv_GetScratchMemi ( inhe );
  tmhe = pkv_GetScratchMem ( inhe*sizeof(BSMhalfedge) );
  tmv = pkv_GetScratchMem ( inv*sizeof(BSMvertex) );
  tmvhei = pkv_GetScratchMemi ( inhe );
  tptc = pkv_GetScratchMemd ( inv*spdimen );
  if ( !tmfac || !tmfhei || !tmhe || !tmv || !tmvhei || !tptc )
    goto failure;

  memcpy ( tmfac, imfac, infac*sizeof(BSMfacet) );
  memcpy ( tmfhei, imfhei, inhe*sizeof(int) );
  memcpy ( tmhe, imhe, inhe*sizeof(BSMhalfedge) );
  memcpy ( tmv, imv, inv*sizeof(BSMvertex) );
  memcpy ( tmvhei, imvhei, inhe*sizeof(int) );
  memcpy ( tptc, iptc, inv*spdimen*sizeof(double) );

  tnv = inv;
  tnhe = inhe;
  tnfac = infac;

  t2mfac = pkv_GetScratchMem ( infac*sizeof(BSMfacet) );
  t2mfhei = pkv_GetScratchMemi ( inhe );
  t2mhe = pkv_GetScratchMem ( inhe*sizeof(BSMhalfedge) );
  t2mv = pkv_GetScratchMem ( inv*sizeof(BSMvertex) );
  t2mvhei = pkv_GetScratchMemi ( inhe );
  t2ptc = pkv_GetScratchMemd ( inv*spdimen );
  if ( !t2mfac || !t2mfhei || !t2mhe || !t2mv || !t2mvhei || !t2ptc )
    goto failure;

  isT = true;

  memcpy ( minxyz, iptc, spdimen*sizeof(double) );
  memcpy ( maxxyz, iptc, spdimen*sizeof(double) );
  for ( i = 1, j = spdimen;  i < inv;  i++, j += spdimen )
    for ( l = 0; l < spdimen; l++ ) {
      if      ( iptc[j+l] < minxyz[l] ) minxyz[l] = iptc[j+l];
      else if ( iptc[j+l] > maxxyz[l] ) maxxyz[l] = iptc[j+l];
     }
  memset ( eps, 0, spdimen*sizeof(double) );
  for ( l = 0; l < spdimen; l++ )
    if ( minxyz[l] < maxxyz[l] )
      eps[l] = (maxxyz[l]-minxyz[l])/nboxes[l];
    else {
      nboxes[l] = 1;
      eps[l] = maxxyz[l]-minxyz[l];
    }

    /* zaalokowanie boksow */
  totalboxes = nboxes[0];
  for ( l = 1; l < spdimen; l++ )
    totalboxes *= nboxes[l];
  Boxes = pkv_GetScratchMemi ( totalboxes );
  if ( !Boxes )
    goto failure;
  Bufor = pkv_InitQueue ( inv, sizeof(int) );
  if ( !Bufor )
    goto failure;
  Queues = pkv_GetScratchMem ( totalboxes*sizeof(pkv_queue*) );
  if ( !Queues )
    goto failure;
  memset ( Queues, 0, totalboxes*sizeof(pkv_queue*) );
        /* obliczenie pojemnosci kolejek dla boksow */
  memset ( Boxes, 0, totalboxes*sizeof(int) );
  for ( v = 0; v < inv; v++ ) {
    memcpy ( xyz, &iptc[spdimen*v], spdimen*sizeof(double) );
    memset ( nxyz, 0, spdimen*sizeof(int) );
    for ( l = 0; l < spdimen; l++ ) {
      nxyz[l] = (xyz[l]-minxyz[l])/eps[l];
      if ( nxyz[l] >= nboxes[l] ) nxyz[l] = nboxes[l]-1;
     }
    ind = nxyz[spdimen-1];
    for ( l = spdimen-2; l >= 0; l-- )
      ind = ind*nboxes[l] + nxyz[l];
    if ( ind >= totalboxes )
      goto failure;
    Boxes[ind] ++;
  }
#ifdef DEBUG1
{
int qcnt = 0;
#endif
  for ( i = 0; i < totalboxes; i++ ) {
    if ( Boxes[i] ) {
      Queues[i] = pkv_InitQueue ( Boxes[i], sizeof(int) );
      if ( !Queues[i] )
        goto failure;
#ifdef DEBUG1
qcnt ++;
#endif
    }
    Boxes[i] = 0;
  }
#ifdef DEBUG1
printf ( "totalboxes = %d, qcnt = %d\n", totalboxes, qcnt );
}
#endif

    /* sprawdzenie ktory wierzcholek jest brzegowy - jak sie dowiem */
    /* jak uzywac to zamienie na bsm_TagBoundaryZoneVertices */
  for ( i = 0; i < inv; i++ ) {
    disabled[i] = false;
    newvnum[i] = i;
    j = imv[i].firsthalfedge+imv[i].degree-1;
    if ( imhe[imvhei[j]].otherhalf < 0 )
      disabled[i] = true;
  }

    /* szukam pierwszego wierzcholka niebrzegowego */
  for ( i = 0; i < inv; i++ )
    if ( !disabled[i] )
      break;
  if ( i >= inv )
    goto failure;

  /* wlozenie do kazdego boksu odpowiednich wierzcholkow */
  pkv_QueueInsert ( Bufor, &i );
  disabled[i] = true;
  do {
    pkv_QueueRemoveFirst ( Bufor, &v );
    memcpy ( xyz, &iptc[spdimen*v], spdimen*sizeof(double) );
    memset ( nxyz, 0, spdimen*sizeof(int) );
    for ( l = 0; l < spdimen; l++ ) {
      nxyz[l] = (xyz[l]-minxyz[l])/eps[l];
      if ( nxyz[l] >= nboxes[l] ) nxyz[l] = nboxes[l]-1;
     }
    ind = nxyz[spdimen-1];
    for ( l = spdimen-2; l >= 0; l-- )
      ind = ind*nboxes[l] + nxyz[l];
    if ( ind >= totalboxes )
      goto failure;
    Boxes[ind] ++;
    pkv_QueueInsert ( Queues[ind], &v );

    for ( i = imv[v].firsthalfedge;
          i < imv[v].degree+imv[v].firsthalfedge;
          i++ ) {
      vnext = imhe[imvhei[i]].v1;
      if ( !disabled[vnext] ) {
        disabled[vnext] = true;
        pkv_QueueInsert ( Bufor, &vnext );
      }
    }
  } while ( !pkv_QueueEmpty ( Bufor ) );

    /* dla kazdego boksu zmniejszamy ilosc wierzcholkow o n-1, */
    /* scian o 2*(n-1), krawedzi o 3*2*(n-1) */
  *onv = inv;
  *onhe = inhe;
  *onfac = infac;
  k = 0;
  for ( i = 0; i < totalboxes; i++ ) {
    if ( Boxes[i] > 1 ) {
        /*przeprowadzam sciskanie krawedzi dla wierzcholkow w Queues[i]*/
      pkv_QueueRemoveFirst ( Queues[i], &vmed );
      memcpy ( xyz, &iptc[vmed*spdimen], spdimen*sizeof(double) );
      ile = 1;
      for ( k = 0; k < Boxes[i]-1; k++ ) {
        pkv_QueueRemoveFirst ( Queues[i], &vnext );
        if ( isT ) {
          fhe = tmv[newvnum[vmed]].firsthalfedge;
          deg = tmv[newvnum[vmed]].degree;
          he = -1;
          while ( he==-1 && k < Boxes[i]-1 ) {
            for ( j = fhe; j < fhe+deg; j++ ) {
              if ( tmhe[tmvhei[j]].v1 == newvnum[vnext] )
                he = tmvhei[j];
            }
                /* vmed nie laczy sie z vnext */
            if ( he == -1 && k < Boxes[i]-1 ) {
              for ( l = 0; l < spdimen; l++ )
                xyz[l] /= (double)ile;
              memcpy ( &tptc[spdimen*newvnum[vmed]], xyz, spdimen*sizeof(double) );
              vmed = vnext;
              ile = 1;
              memcpy ( xyz, &iptc[vmed*spdimen], spdimen*sizeof(double) );
              pkv_QueueRemoveFirst ( Queues[i], &vnext );
              k ++;
              fhe = tmv[newvnum[vmed]].firsthalfedge;
              deg = tmv[newvnum[vmed]].degree;
            }
          }
          if ( he != -1 ) {
#ifdef DEBUG1
printf ( "%i\n", he );
#endif
            bsm_ContractEdgeNum ( tnv, tmv, tmvhei, tnhe, tmhe, tnfac, tmfac, tmfhei,
                                  he, &t2nv, &t2nhe, &t2nfac );
            err = bsm_ContractEdged ( spdimen, tnv, tmv, tmvhei, tptc, tnhe, tmhe,
                                      tnfac, tmfac, tmfhei, he, &t2nv, t2mv, t2mvhei,
                                      t2ptc, &t2nhe, t2mhe, &t2nfac, t2mfac, t2mfhei );
#ifdef DEBUG
                            /*printf("error: %i\n", err); */
sprintf ( filename, "k%d.bs", k );
bsf_OpenOutputFile ( filename, false );
bsf_WriteBSMeshd ( 3, 3, false, 1, t2nv, t2mv, t2mvhei, t2ptc, t2nhe, t2mhe,
                   t2nfac, t2mfac, t2mfhei, NULL, 0, NULL, NULL );
bsf_CloseOutputFile ();
#endif
            if ( err >= 0 ) {
              if ( !bsm_CheckMeshIntegrity ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                                             t2nfac, t2mfac, t2mfhei, NULL, NULL ) )
                err = -1;
            }
            if ( err != -1 ) {
              ile ++;
              for ( l = 0; l < spdimen; l++ )
                xyz[l] += iptc[spdimen*vnext+l];
            }
          }
        }
        else {
          fhe = t2mv[newvnum[vmed]].firsthalfedge;
          deg = t2mv[newvnum[vmed]].degree;
          he = -1;
          while ( he == -1 && k < Boxes[i]-1 ) {
            for ( j = fhe; j < fhe+deg; j++ ) {
              if ( t2mhe[t2mvhei[j]].v1 == newvnum[vnext] )
                he = t2mvhei[j];
            }
                         /* vmed nie laczy sie z vnext */
            if ( he == -1 && k < Boxes[i]-1 ) {
              for ( l = 0; l < spdimen; l++ )
                xyz[l] /= (double)ile;
              memcpy ( &t2ptc[spdimen*newvnum[vmed]], xyz, spdimen*sizeof(double) );
              vmed = vnext;
              ile = 1;
              memcpy ( xyz, &iptc[spdimen*vmed], spdimen*sizeof(double) );
              pkv_QueueRemoveFirst ( Queues[i], &vnext );
              k++;
              fhe = t2mv[newvnum[vmed]].firsthalfedge;
              deg = t2mv[newvnum[vmed]].degree;
            }
          }
          if ( he != -1 ) {
#ifdef DEBUG1
printf ( "%i\n", he );
#endif
            bsm_ContractEdgeNum ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe, t2nfac, t2mfac,
                                  t2mfhei, he, &tnv, &tnhe, &tnfac );
            err = bsm_ContractEdged ( spdimen, t2nv, t2mv, t2mvhei, t2ptc, t2nhe,
                                      t2mhe, t2nfac, t2mfac, t2mfhei, he,
                                      &tnv, tmv, tmvhei, tptc, &tnhe, tmhe, &tnfac,
                                      tmfac, tmfhei );
#ifdef DEBUG
/*printf ( "error: %i\n", err ); */
/*sprintf ( filename, "k%d.bs", k ); */
/*bsf_OpenOutputFile ( filename, false ); */
/*bsf_WriteBSMeshd ( 3, 3, false, 1, t2nv, t2mv, t2mvhei,
  t2ptc, t2nhe, t2mhe, t2nfac, t2mfac, t2mfhei, NULL, NULL, NULL ); */
/*bsf_CloseOutputFile (); */
#endif
            if ( err >= 0 ) {
              if ( !bsm_CheckMeshIntegrity ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                                             t2nfac, t2mfac, t2mfhei, NULL, NULL ) )
                err = -1;
            }
            if ( err != -1 ) {
              ile ++;
              for ( l = 0; l < spdimen; l++ )
                xyz[l] += iptc[spdimen*vnext+l];
            }
          }
        }

        if ( he != -1 && err != -1 ) {
          isT = !isT;
                        /* przenumerowuje wierzcholki */
          newvnum[vnext] = -1;
          for ( j = vnext+1; j < inv; j++ )
            newvnum[j]--;
        }
      }
      for ( l = 0; l < spdimen; l++ )
        xyz[l] /= (double)ile;

      if ( isT ) {
        memcpy ( &tptc[spdimen*newvnum[vmed]], xyz, spdimen*sizeof(double) );
        bsm_CheckMeshIntegrity ( tnv, tmv, tmvhei, tnhe, tmhe,
                                 tnfac, tmfac, tmfhei, NULL, NULL );
      }
      else {
        memcpy ( &t2ptc[spdimen*newvnum[vmed]], xyz, spdimen*sizeof(double) );
        bsm_CheckMeshIntegrity ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                                 t2nfac, t2mfac, t2mfhei, NULL, NULL );
      }
    }
  }

  if ( isT ) {
    if ( !bsm_CheckMeshIntegrity ( tnv, tmv, tmvhei, tnhe, tmhe,
                                   tnfac, tmfac, tmfhei, NULL, NULL ) ) {
      if ( bsm_CheckMeshIntegrity ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                                    t2nfac, t2mfac, t2mfhei, NULL, NULL ) )
        isT=false;
      else
        goto failure;
    }
  }
  else {
    if ( !bsm_CheckMeshIntegrity ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                                   t2nfac, t2mfac, t2mfhei, NULL, NULL ) ) {
      if ( bsm_CheckMeshIntegrity ( tnv, tmv, tmvhei, tnhe, tmhe,
                                    tnfac, tmfac, tmfhei, NULL, NULL ) )
        isT=true;
      else
        goto failure;
    }
  }

  PKV_FREE ( Bufor );
  for ( i = 0; i < totalboxes; i++ )
    if ( Queues[i] ) PKV_FREE ( Queues[i] );

  if ( isT ) { /* najnowsza siatka znajduje sie w t */
    t2nv = tnv;
    t2nhe = tnhe;
    t2nfac = tnfac;
    t2mfac = tmfac;
    t2mfhei = tmfhei;
    t2mhe = tmhe;
    t2mv = tmv;
    t2mvhei = tmvhei;
    t2ptc = tptc;
  }
  PKV_MALLOC ( *omv, t2nv*sizeof(BSMvertex) );
  PKV_MALLOC ( *omvhei, t2nhe*sizeof(int) );
  PKV_MALLOC ( *omfhei, t2nhe*sizeof(int) );
  PKV_MALLOC ( *omhe, t2nhe*sizeof(BSMhalfedge) );
  PKV_MALLOC ( *optc, t2nv*spdimen*sizeof(double) );
  PKV_MALLOC ( *omfac, t2nfac*sizeof(BSMfacet) );
  if ( !*omv || !*omvhei || !*omfhei || !*omhe || !*optc || !*omfac )
    goto failure;
  *onv = t2nv;
  *onhe = t2nhe;
  *onfac = t2nfac;
  memcpy ( *omfac, t2mfac, *onfac*sizeof(BSMfacet) );
  memcpy ( *omfhei, t2mfhei, *onhe*sizeof(int) );
  memcpy ( *omhe, t2mhe, *onhe*sizeof(BSMhalfedge) );
  memcpy ( *omv, t2mv, *onv*sizeof(BSMvertex) );
  memcpy ( *omvhei, t2mvhei, *onhe*sizeof(int) );
  memcpy ( *optc, t2ptc, t2nv*spdimen*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( *omv )    PKV_FREE ( *omv );
  if ( *omvhei ) PKV_FREE ( *omvhei );
  if ( *omfhei ) PKV_FREE ( *omfhei );
  if ( *optc )   PKV_FREE ( *optc );
  if ( *omhe )   PKV_FREE ( *omhe );
  if ( *omfac )  PKV_FREE ( *omfac );
  if ( Bufor )   PKV_FREE ( Bufor );
  if ( Queues ) {
    for ( i = 0; i < totalboxes; i++ )
      if ( Queues[i] ) PKV_FREE ( Queues[i] );
  }
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_SimplifyMeshd*/

