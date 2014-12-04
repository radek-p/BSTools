
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef BSMESH_H
#define BSMESH_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

/* error identifiers */
#define BSM_ERROR_NO_SCRATCH_MEMORY           1
#define BSM_ERROR_NOT_ENOUGH_ELEMENTS         2
#define BSM_ERROR_INVALID_MVHEI               3
#define BSM_ERROR_INVALID_MFHEI               4
#define BSM_ERROR_INVALID_PERMUTATION         5
#define BSM_ERROR_INVALID_HE_V0V1             6
#define BSM_ERROR_INVALID_HE_FACETNUM         7
#define BSM_ERROR_INVALID_HE_OTHERHALF        8
#define BSM_ERROR_INVALID_HE_PAIR             9
#define BSM_ERROR_INVALID_HE_PAIR_V0V1       10
#define BSM_ERROR_INVALID_VERTEX_DEGREE      11
#define BSM_ERROR_INVALID_VERTEX_INNERHE     12
#define BSM_ERROR_INVALID_VERTEX_LASTHE      13
#define BSM_ERROR_INVALID_INNER_VERTEX_DEG   14
#define BSM_ERROR_INVALID_FACET_DEGREE       15
#define BSM_ERROR_INVALID_FACET_VERT         16
#define BSM_ERROR_INVALID_FACET_HEDGE        17
#define BSM_ERROR_INVALID_FACET_ORIENTATON   18
#define BSM_ERROR_INVALID_VERTEX_ORIENTATION 19
#define BSM_ERROR_DOUBLE_FACET_VERTEX        20
#define BSM_ERROR_DOUBLE_VERTEX_FACET        21

/* mesh representation data structure */
typedef struct {
    char degree;
    int  firsthalfedge;
  } BSMfacet, BSMvertex;

typedef struct {
    int v0, v1;
    int facetnum;
    int otherhalf;
  } BSMhalfedge;

/* Sabin net list data structure */
typedef struct {
    byte  el_type;               /* mesh special element type, */
                                 /* 0 - vertex, 1 - facet */
    byte  degree;                /* degree of the special vertex or facet */
    byte  snet_rad;              /* Sabin net radius */
    short snet_nvert;            /* number of Sabin net vertices */
    int   first_snet_vertex;     /* index of the first Sabin net vertex */
  } bsm_special_el;

typedef struct {
    int            nspecials;  /* number of special mesh elements */
    int            nspvert;    /* number of vertices in the Sabin nets */
    int            nextravert; /* number of extra vertices */
    bsm_special_el *spel;      /* array of special elements */
    int            *spvert;    /* array of indices of the vertices */
                               /* in the Sabin nets */
  } bsm_special_elem_list;


boolean bsm_CheckMeshIntegrity ( int nv, const BSMvertex *mv, const int *mvhei,
                                 int nhe, const BSMhalfedge *mhe,
                                 int nfac, const BSMfacet *mfac, const int *mfhei,
                                 int *errcode, int *errelem );

void bsm_TagMesh ( int nv, BSMvertex *mv, int *mvhei,
                   int nhe, BSMhalfedge *mhe,
                   int nfac, BSMfacet *mfac, int *mfhei,
                   char *vtag, char *ftag,
                   int *vi, int *vb, int *ei, int *eb );

boolean bsm_DoublingNum ( int inv, BSMvertex *imv, int *imvhei,
                          int inhe, BSMhalfedge *imhe,
                          int infac, BSMfacet *imfac, int *imfhei,
                          int *onv, int *onhe, int *onfac );
boolean bsm_Doublingd ( int spdimen,
                        int inv, BSMvertex *imv, int *imvhei, double *iptc,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int *onv, BSMvertex *omv, int *omvhei, double *optc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei );

int bsm_DoublingMatSize ( int inv, BSMvertex *imv, int *imvhei,
                          int inhe, BSMhalfedge *imhe,
                          int infac, BSMfacet *imfac, int *imfhei );
boolean bsm_DoublingMatd ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int *onv, BSMvertex *omv, int *omvhei,
                           int *onhe, BSMhalfedge *omhe,
                           int *onfac, BSMfacet *omfac, int *omfhei,
                           int *ndmat, index2 *dmi, double *dmc );

boolean bsm_AveragingNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int *onv, int *onhe, int *onfac );
boolean bsm_Averagingd ( int spdimen,
                         int inv, BSMvertex *imv, int *imvhei, double *iptc,
                         int inhe, BSMhalfedge *imhe,
                         int infac, BSMfacet *imfac, int *imfhei,
                         int *onv, BSMvertex *omv, int *omvhei, double *optc,
                         int *onhe, BSMhalfedge *omhe,
                         int *onfac, BSMfacet *omfac, int *omfhei );

int bsm_AveragingMatSize ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei );
boolean bsm_AveragingMatd ( int inv, BSMvertex *imv, int *imvhei,
                            int inhe, BSMhalfedge *imhe,
                            int infac, BSMfacet *imfac, int *imfhei,
                            int *onv, BSMvertex *omv, int *omvhei,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei,
                            int *namat, index2 *ami, double *amc );

boolean bsm_RefineBSMeshd ( int spdimen, int degree,
                    int inv, BSMvertex *imv, int *imvhei, double *iptc,
                    int inhe, BSMhalfedge *imhe,
                    int infac, BSMfacet *imfac, int *imfhei,
                    int *onv, BSMvertex **omv, int **omvhei, double **optc,
                    int *onhe, BSMhalfedge **omhe,
                    int *onfac, BSMfacet **omfac, int **omfhei );

boolean bsm_RefinementMatd ( int degree,
                             int inv, BSMvertex *imv, int *imvhei,
                             int inhe, BSMhalfedge *imhe,
                             int infac, BSMfacet *imfac, int *imfhei,
                             int *onv, BSMvertex **omv, int **omvhei,
                             int *onhe, BSMhalfedge **omhe,
                             int *onfac, BSMfacet **omfac, int **omfhei,
                             int *nrmat, index2 **rmi, double **rmc );

void bsm_MergeMeshesd ( int spdimen,
                        int nv1, BSMvertex *mv1, int *mvhei1, double *vpc1,
                        int nhe1, BSMhalfedge *mhe1,
                        int nfac1, BSMfacet *mfac1, int *mfhei1,
                        int nv2, BSMvertex *mv2, int *mvhei2, double *vpc2,
                        int nhe2, BSMhalfedge *mhe2,
                        int nfac2, BSMfacet *mfac2, int *mfhei2,
                        int *onv, BSMvertex *omv, int *omvhei, double *ovpc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei );

boolean bsm_RemoveFacetNum ( int inv, BSMvertex *imv, int *imvhei,
                             int inhe, BSMhalfedge *imhe,
                             int infac, BSMfacet *imfac, int *imfhei,
                             int nfr,
                             int *onv, int *onhe, int *onfac );
boolean bsm_RemoveFacetd ( int spdimen,
                           int inv, BSMvertex *imv, int *imvhei, double *iptc,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nfr,
                           int *onv, BSMvertex *omv, int *omvhei, double *optc,
                           int *onhe, BSMhalfedge *omhe,
                           int *onfac, BSMfacet *omfac, int *omfhei );

void bsm_FacetEdgeDoublingNum ( int inv, BSMvertex *imv, int *imvhei,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int fn,
                                int *onv, int *onhe, int *onfac );
boolean bsm_FacetEdgeDoublingd ( int spdimen,
                                 int inv, BSMvertex *imv, int *imvhei, double *iptc,
                                 int inhe, BSMhalfedge *imhe,
                                 int infac, BSMfacet *imfac, int *imfhei,
                                 int fn,
                                 int *onv, BSMvertex *omv, int *omvhei, double *optc,
                                 int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei );

void bsm_RemoveVertexNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nvr,
                           int *onv, int *onhe, int *onfac );
boolean bsm_RemoveVertexd ( int spdimen,
                            int inv, BSMvertex *imv, int *imvhei, double *iptc,
                            int inhe, BSMhalfedge *imhe,
                            int infac, BSMfacet *imfac, int *imfhei,
                            int nvr,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei );

void bsm_ContractEdgeNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int nche,
                           int *onv, int *onhe, int *onfac );
int bsm_ContractEdged ( int spdimen,
                        int inv, BSMvertex *imv, int *imvhei, double *iptc,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int nche,
                        int *onv, BSMvertex *omv, int *omvhei, double *optc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei );

/* the procedures with the headers below search regular (i.e. rectangular) */
/* and special subnets of indicated sizes in the mesh */
boolean bsm_FindRegularSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                 int nhe, BSMhalfedge *mhe,
                                 int nfac, BSMfacet *mfac, int *mfhei,
                                 int d, void *usrptr,
                                 void (*output)( int d, int *vertnum, int *mtab,
                                                 void *usrptr ) );
boolean bsm_FindSpecialVSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  int d, void *usrptr,
                                  void (*output)( int d, int k, int *vertnum,
                                                  int *mtab, void *usrptr ) );
boolean bsm_FindSpecialFSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  int d, void *usrptr,
                                  void (*output)( int d, int k, int *vertnum,
                                                  int *mtab, void *usrptr ) );

boolean bsm_CountSpecialVSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte snet_rad,
                                   int *nspecials, int *nspvert );
boolean bsm_FindSpecialVSubnetList ( int nv, BSMvertex *mv, int *mvhei,
                                     int nhe, BSMhalfedge *mhe,
                                     int nfac, BSMfacet *mfac, int *mfhei,
                                     byte snet_rad,
                                     boolean append,
                                     bsm_special_elem_list *list );

boolean bsm_CountSpecialFSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte snet_rad,
                                   int *nspecials, int *nspvert );
boolean bsm_FindSpecialFSubnetLists ( int nv, BSMvertex *mv, int *mvhei,
                                      int nhe, BSMhalfedge *mhe,
                                      int nfac, BSMfacet *mfac, int *mfhei,
                                      boolean append,
                                      byte snet_rad,
                                      bsm_special_elem_list *list );

/* the procedure with the header below finds vertices distant from */
/* the mesh boundary by at most d edges */
void bsm_TagBoundaryZoneVertices ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   char d, char *vtag );

/* the procedures below assign each vertex a distance from a specified vertex */
/* 1. measured with the number of edges */
boolean bsm_FindVertexDistances1 ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   int v, boolean back, int *dist );

/* 2. measured with the number of facets (assuming that facets with a common) */
/*    vertex are neighbours */
boolean bsm_FindVertexDistances2 ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   int v, int *dist );

/* procedures of changing the mesh topology */
int bsm_HalfedgeLoopLength ( int nv, BSMvertex *mv, int *mvhei,
                             int nhe, BSMhalfedge *mhe,
                             int he );

boolean bsm_GlueHalfedgeLoopsd ( int spdimen,
                                 int inv, BSMvertex *imv, int *imvhei, double *ivc,
                                 int inhe, BSMhalfedge *imhe,
                                 int infac, BSMfacet *imfac, int *imfhei,
                                 int he1, int he2,
                                 int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                                 int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei );

/* procedures extracting a part of a mesh */
boolean bsm_ExtractSubmeshVNum ( int inv, const BSMvertex *imv, const int *imvhei,
                                 int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 boolean *vtag,
                                 int *onv, int *onhe, int *onfac );

boolean bsm_ExtractSubmeshVd ( int spdimen, int inv,
                               const BSMvertex *imv, const int *imvhei,
                               double *iptc, int inhe, const BSMhalfedge *imhe,
                               int infac, const BSMfacet *imfac, const int *imfhei,
                               boolean *vtag,
                               int *onv, BSMvertex *omv, int *omvhei,
                               double *optc, int *onhe, BSMhalfedge *omhe,
                               int *onfac, BSMfacet *omfac, int *omfhei );

/* procedures sealing a hole in a mesh with a single facet */
boolean bsm_SealMeshHoleNum ( int inv, const BSMvertex *imv, const int *imvhei,
                              int inhe, const BSMhalfedge *imhe,
                              int infac, const BSMfacet *imfac, const int *imfhei,
                              int nbhe,
                              int *onv, int *onhe, int *onfac );

boolean bsm_SealMeshHoled ( int spdimen,
                            int inv, const BSMvertex *imv, const int *imvhei,
                            const double *iptc,
                            int inhe, const BSMhalfedge *imhe,
                            int infac, const BSMfacet *imfac, const int *imfhei,
                            int nbhe,
                            int *onv, BSMvertex *omv, int *omvhei, double *optc,
                            int *onhe, BSMhalfedge *omhe,
                            int *onfac, BSMfacet *omfac, int *omfhei );

/* procedures glueing two boundary halfedges */
boolean bsm_GlueTwoHalfedgesNum ( int inv, const BSMvertex *imv, const int *imvhei,
                                  int inhe, const BSMhalfedge *imhe,
                                  int infac, const BSMfacet *imfac, const int *imfhei,
                                  int he1, int he2,
                                  int *onv, int *onhe, int *onfac );

boolean bsm_GlueTwoHalfedgesd ( int spdimen,
                                int inv, const BSMvertex *imv, const int *imvhei,
                                const double *ivc,
                                int inhe, const BSMhalfedge *imhe,
                                int infac, const BSMfacet *imfac, const int *imfhei,
                                int he1, int he2,
                                int *onv, BSMvertex *omv, int *omvhei, double *ovc,
                                int *onhe, BSMhalfedge *omhe,
                                int *onfac, BSMfacet *omfac, int *omfhei );

/* splitting a boundary edge to halves */
boolean bsm_SplitBoundaryEdged ( int spdimen, int inv,
                                 const BSMvertex *imv, const int *imvhei,
                                 double *iptc, int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 int splithe,
                                 int *onv, BSMvertex *omv, int *omvhei,
                                 double *optc, int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei );

/* splitting a facet between two vertices */
boolean bsm_DivideFacetd ( int spdimen, int inv,
                           const BSMvertex *imv, const int *imvhei,
                           double *iptc, int inhe, const BSMhalfedge *imhe,
                           int infac, const BSMfacet *imfac, const int *imfhei,
                           int nV0, int nV1,
                           int *onv, BSMvertex *omv, int *omvhei,
                           double *optc, int* onhe, BSMhalfedge *omhe,
                           int *onfac, BSMfacet *omfac, int *omfhei );

/* doubling edges forming a closed loop */
void bsm_EdgeLoopDoublingNum ( int inv, int inhe, int infac,
                               int EdgeLoopLength,
                               int *onv, int *onhe, int *onfac );

boolean bsm_EdgeLoopDoublingd ( int spdimen,
                                int inv, BSMvertex *imv, int *imvhei, double *iptc,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int EdgeLoopLength, int *EdgeLoop,
                                int *onv, BSMvertex *omv, int *omvhei, double *optc,
                                int *onhe, BSMhalfedge *omhe,
                                int *onfac, BSMfacet *omfac, int *omfhei );

/* dividing all facets to triangles */
void bsm_TriangulateFacetsNum ( int inv, BSMvertex *imv, int *imvhei,
                                int inhe, BSMhalfedge *imhe,
                                int infac, BSMfacet *imfac, int *imfhei,
                                int *onv, int *onhe, int *onfac );

boolean bsm_TriangulateFacetsd ( int spdimen, int inv,
                                 const BSMvertex *imv, const int *imvhei,
                                 double *iptc, int inhe, const BSMhalfedge *imhe,
                                 int infac, const BSMfacet *imfac, const int *imfhei,
                                 int *onv, BSMvertex *omv, int *omvhei,
                                 double *optc, int *onhe, BSMhalfedge *omhe,
                                 int *onfac, BSMfacet *omfac, int *omfhei );

/* mesh simplification - experimental */
boolean bsm_SimplifyMeshd ( int spdimen, int inv,
                            const BSMvertex *imv, const int *imvhei,
                            double *iptc, int inhe, const BSMhalfedge *imhe,
                            int infac, const BSMfacet *imfac, const int *imfhei,
                            int *nboxes,
                            int *onv, BSMvertex **omv, int **omvhei,
                            double **optc, int *onhe, BSMhalfedge **omhe,
                            int *onfac, BSMfacet **omfac, int **omfhei );

#ifdef __cplusplus
}
#endif

#endif /*BSMESH_H*/

