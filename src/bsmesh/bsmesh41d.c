
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Anna Sierhej                                     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "bsmesh.h"

#include "bsmprivate.h"

#define DEBUG
/* ////////////////////////////////////////////////////////////////////////// */
static boolean Compare ( void* v1, void* v2 )
{
  return ((vector2d*)v1)->y < ((vector2d*)v2)->y;
} /*Compare*/

static vector2d getError ( int inv, const BSMvertex *imv, const int *imvhei,
                           double *iptc, int inhe,const BSMhalfedge *imhe,
                           int infac, const BSMfacet *imfac, const int *imfhei,
                           int i )
{
  int      j, fhe, deg, facet, facetFhe, v1, v2, v3;
  vector3d V1, V2, Vn;
  point3d  vc;
  double   x1, y1, z1, x2, y2, z2, x3, y3, z3, area, sumA, distance;
  vector3d sumB, sumN;
  point2d  pom;

  fhe = imv[i].firsthalfedge;
  deg = imv[i].degree;
  sumA = 0;
  SetPoint3d ( &sumB, 0.0, 0.0, 0.0 );
  SetVector3d ( &sumN, 0.0, 0.0, 0.0 );

  for ( j = fhe; j < fhe+deg; j++ ) {
    facet = imhe[imvhei[j]].facetnum;
    facetFhe = imfac[facet].firsthalfedge;
    v1 = imhe[imfhei[facetFhe]].v0;
    v2 = imhe[imfhei[facetFhe+1]].v0;
    v3 = imhe[imfhei[facetFhe+2]].v0;

    x1 = iptc[3*v1];
    y1 = iptc[3*v1+1];
    z1 = iptc[3*v1+2];
    x2 = iptc[3*v2];
    y2 = iptc[3*v2+1];
    z2 = iptc[3*v2+2];
    x3 = iptc[3*v3];
    y3 = iptc[3*v3+1];
    z3 = iptc[3*v3+2];

    SetVector3d ( &V1, x1-x2, y1-y2, z1-z2 );
    SetVector3d ( &V2, x1-x3, y1-y3, z1-z3 );
    SetPoint3d ( &vc, (x1+x2+x3)/3.0, (y1+y2+y3)/3.0, (z1+z2+z3)/3.0 );
    CrossProduct3d ( &V1, &V2, &Vn );

    area = Vn.x*Vn.x + Vn.y*Vn.y +Vn.z*Vn.z;
    NormalizeVector3d ( &Vn );
    sumA += area;
    SetPoint3d ( &sumB, sumB.x + vc.x*area, sumB.y + vc.y*area, sumB.z + vc.z*area );
    SetVector3d ( &sumN, sumB.x + vc.x*area, sumB.y + vc.y*area, sumB.z + vc.z*area );
  }

  SetPoint3d ( &sumB, sumB.x/sumA, sumB.y/sumA, sumB.z/sumA );
  SetVector3d ( &sumN, sumN.x/sumA, sumN.y/sumA, sumN.z/sumA );
  SetVector3d ( &Vn, iptc[3*i]-sumB.x, iptc[3*i+1]-sumB.y, iptc[3*i+2]-sumB.z );
  distance = fabs ( DotProduct3d ( &sumN, &Vn ) );
  SetPoint2d ( &pom, i, distance );
  return pom;
} /*getError*/

boolean bsm_VertexDecimation ( int spdimen, int inv,
                               const BSMvertex *imv, const int *imvhei,
                               double *iptc, int inhe,const BSMhalfedge *imhe,
                               int infac, const BSMfacet *imfac, const int *imfhei,
                               int nIterations,
                               int *onv, BSMvertex **omv, int **omvhei,
                               double **optc, int *onhe, BSMhalfedge **omhe,
                               int *onfac, BSMfacet **omfac, int **omfhei )
{
  int         i, j, v, nHeap;
  point2d     **Error;
  point2d     *pom;

  int         *newvnum;
  void        *sp;
  boolean     *disabled;
  int         tnv;
  BSMvertex   *tmv;
  int         *tmvhei;
  double      *tptc;
  int         tnhe;
  BSMhalfedge *tmhe;
  int         tnfac;
  BSMfacet    *tmfac;
  int         *tmfhei;

  int         t2nv;
  BSMvertex   *t2mv;
  int         *t2mvhei;
  double      *t2ptc;
  int         t2nhe;
  BSMhalfedge *t2mhe;
  int         t2nfac;
  BSMfacet    *t2mfac;
  int         *t2mfhei;

  boolean     isT;

  sp = pkv_GetScratchMemTop ();

        /*workspace allocation*/
  newvnum = pkv_GetScratchMemi ( inv );
  disabled = pkv_GetScratchMem ( inv*sizeof(boolean) );
  tmfac = pkv_GetScratchMem ( infac*sizeof(BSMfacet) );
  tmfhei = pkv_GetScratchMemi ( inhe );
  tmhe = pkv_GetScratchMem ( inhe*sizeof(BSMhalfedge) );
  tmv = pkv_GetScratchMem ( inv*sizeof(BSMvertex) );
  tmvhei = pkv_GetScratchMemi ( inhe );
  tptc = pkv_GetScratchMemd ( inv*spdimen );
  t2mfac = pkv_GetScratchMem ( infac*sizeof(BSMfacet) );
  t2mfhei = pkv_GetScratchMemi ( inhe );
  t2mhe = pkv_GetScratchMem ( inhe*sizeof(BSMhalfedge) );
  t2mv = pkv_GetScratchMem ( inv*sizeof(BSMvertex) );
  t2mvhei = pkv_GetScratchMemi ( inhe );
  t2ptc = pkv_GetScratchMemd ( inv*spdimen );
  Error = pkv_GetScratchMem ( inv*sizeof(point2d*) );
  pom = pkv_GetScratchMem ( inv*sizeof(point2d) );
  if ( !newvnum || !disabled || !tmfac || !tmfhei || !tmhe || !tmv || !tmvhei ||
       !tptc || !t2mfac || !t2mfhei || !t2mhe || !t2mv || !t2mvhei || !t2ptc ||
       !Error || !pom )
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

  isT = true;
  nHeap = -1;

        /* marking boundary vertices */
  for ( i = 0; i < inv; i++ ) {
    newvnum[i] = i;
    j = imv[i].firsthalfedge+imv[i].degree-1;
    disabled[i] = imhe[imvhei[j]].otherhalf < 0;
  }

  for ( i = 0; i < inv; i++ )
    if ( !disabled[i] ) {
      pom[i] = getError ( inv, imv, imvhei, iptc, inhe, imhe,
                          infac, imfac, imfhei, i );
      pkv_HeapInsert ( (void**)Error, &nHeap, &pom[i], Compare );
    }
  pkv_HeapSort ( (void**)Error, nHeap, Compare );

  for ( i = 0; i < nIterations; i++ ) {
    v = newvnum[(int)(*Error[i]).x];
    if ( isT ) {
      bsm_RemoveVertexNum ( tnv, tmv, tmvhei, tnhe, tmhe, tnfac,
                            tmfac, tmfhei, v, &t2nv, &t2nhe, &t2nfac );
      bsm_RemoveVertexd ( spdimen, tnv, tmv, tmvhei, tptc, tnhe, tmhe,
                          tnfac, tmfac, tmfhei, v,
                          &t2nv, t2mv, t2mvhei, t2ptc, &t2nhe, t2mhe,
                          &t2nfac, t2mfac, t2mfhei );
      if ( t2nfac > 0 ) {
#ifdef DEBUG
printf ( "%i \n", v );
#endif
        bsm_TriangulateFacetd ( spdimen, &t2nv, t2mv, t2mvhei, t2ptc,
                                &t2nhe, t2mhe, &t2nfac, t2mfac, t2mfhei, t2nfac-1 );
        for ( j = v+1; j < inv; j++ )
          newvnum[j] --;
        isT = !isT;
      }
      newvnum[(int)(*Error[i]).x] = -1;
    }
    else {
      bsm_RemoveVertexNum ( t2nv, t2mv, t2mvhei, t2nhe, t2mhe,
                            t2nfac, t2mfac, t2mfhei, v, &tnv, &tnhe, &tnfac);
      bsm_RemoveVertexd ( spdimen, t2nv, t2mv, t2mvhei, t2ptc, t2nhe, t2mhe,
                          t2nfac, t2mfac, t2mfhei, v,&tnv, tmv, tmvhei, tptc,
                          &tnhe, tmhe, &tnfac, tmfac, tmfhei );
      if ( tnfac > 0 ) {
#ifdef DEBUG
printf("%i \n", v);
#endif
        bsm_TriangulateFacetd ( spdimen, &tnv, tmv, tmvhei, tptc,
                                &tnhe, tmhe, &tnfac, tmfac, tmfhei, tnfac-1 );
        for ( j = v+1; j < inv; j++ ) {
          newvnum[j] --;
        }
        isT = !isT;
      }
      newvnum[(int)(*Error[i]).x] = -1;
    }
        /* todo - processing of other vertices */
        /* ************* */
  }

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
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_VertexDecimation*/

