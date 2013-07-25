
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef G2MBLENDINGD_H
#define G2MBLENDINGD_H

#ifndef EGHOLED_H
#include "egholed.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define G2MBL_TIME_IT

/* this constant is relevant for the two-level block algorithm */
#define G2MBL_MAX_BLOCKS 32   /* int is assumed at least 32 bits */

/* this constant is relevant for the multilevel block algorithm */
#define G2MBL_MAX_LEVELS 7

/* this variable is initially set to 1, but in multiprocessor systems   */
/* it may be assigned the number of processors; then parallel algoritms */
/* will be used (if they exist - most of them are still to be written)  */
extern int _g2mbl_npthreads;

extern GHoleDomaind *g2mbl_domaind[GH_MAX_K-3];
extern double *g2mbl_patchmatrixd[GH_MAX_K-3];
extern void (*g2mbl_outputnzdistr)( int nbl, int blnum, boolean final,
                                    int nvcp, int n, byte *nzdistr );

void g2mbl_CleanupHoleDomainsd ( void );
boolean g2mbl_SetupHolePatchMatrixd ( int k );


void g2mbl_OptLMTDeallocated ( void **data );

/* procedure counting the variable vertices of the mesh */
int g2mbl_GetNvcp ( int nv, BSMvertex *mv, int *mvhei,
                    int nhe, BSMhalfedge *mhe,
                    int nfac, BSMfacet *mfac, int *mfhei,
                    byte *mkcp );
/* procedures for the simplest version nonblock optimization */
boolean g2mbl_InitBlSurfaceOptLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                     point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                     int nfac, BSMfacet *mfac, int *mfhei,
                                     byte *mkcp,
                                     double C, double dO, double dM,
                                     int nkn1, int nkn2, void **data );
boolean g2mbl_IterBlSurfaceOptLMTd ( void *data, boolean *finished );
boolean g2mbl_FindBlSurfaceLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                  point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  byte *mkcp,
                                  double C, double dO, double dM,
                                  int maxit, int nkn1, int nkn2 );

/* procedures for the simplest version block optimization */
boolean g2mbl_InitBlSurfaceOptAltBLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                         point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                         int nfac, BSMfacet *mfac, int *mfhei,
                                         byte *mkcp,
                                         double C, double dO, double dM,
                                         int nkn1, int nkn2, int nbl,
                                         void **data );

boolean g2mbl_InitBlCMPSurfaceOptd ( /* fine mesh */
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei,
                    byte *fmkcp,
                           /* number of vertices of the coarse mesh */
                    int cnv,
                           /* refinement matrix */
                    int rmnnz, index2 *rmnzi, double *rmnzc,
                           /* optimization parameters */
                    double C, double dO, double dM,
                    int nkn1, int nkn2, int nbl,
                           /* created data structure */
                    void **data );

boolean g2mbl_IterBlSurfaceOptAltBLMTd ( void *data, boolean *finished );

boolean g2mbl_FindBlSurfaceAltBLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                      point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                      int nfac, BSMfacet *mfac, int *mfhei,
                                      byte *mkcp,
                                      double C, double dO, double dM,
                                      int maxit, int nkn1, int nkn2, int nbl );

/* auxiliary procedures, for various tests etc. */
boolean g2mbl_TimeBlSurfaceOptBLMTd ( void *data, int bnum, double *tt,
                                      int *n0, int *n, int **prof,
                                      double ***rows, double ***lrows,
                                      double *_func, double *_grad,
                                      int **iHbl, int **tHbl, double **Hbl );

int g2mbl_GetBLMBlockNumd ( void *data, int *lastblock );
void g2mbl_GetBLMTBlockInfod ( void *data,
                               int bln, int *nv, int *nvcp, int **nncpi,
                               int *c0, int *bnvcp,
                               int **vncpi, int **bvncpi, int **vpermut );

/* multilevel optimization algorithm */
boolean g2mbl_MLOptInitd ( int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                           int nhe, BSMhalfedge *mhe,
                           int nfac, BSMfacet *mfac, int *mfhei,
                           byte *mkcp,
                           double C, double dO, double dM,
                           int nkn1, int nkn2, short nlevels, void **data );
boolean g2mbl_MLCMPOptInitd ( /* fine mesh */
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei,
                    byte *fmkcp,
                           /* number of vertices of the coarse mesh */
                    int cnv,
                           /* refinement matrix */
                    int rmnnz, index2 *rmnzi, double *rmnzc,
                           /* optimization parameters */
                    double C, double dO, double dM,
                    int nkn1, int nkn2, short nlevels,
                           /* created data structure */
                    void **data );
void g2mbl_MLOptDeallocated ( void **data );
boolean g2mbl_MLOptIterd ( void *data, boolean *finished );
boolean g2mbl_MLCOptIterd ( void *data, boolean *finished );

void g2mbl_MLSetLogLeveld ( void *data, short level );
short g2mbl_MLGetInfod ( void *data );

int g2mbl_MLGetLastBlockd ( void *data );
boolean g2mbl_MLGetBlockVCPNumbersd ( void *data, int bl,
                                      int *nvcp, int **vncpi, int *seed );

boolean g2mbl_MLSuggestNLevels ( int nv, BSMvertex *mv, int *mvhei,
                                 int nhe, BSMhalfedge *mhe,
                                 int nfac, BSMfacet *mfac, int *mfhei,
                                 byte *mkcp,
                                 int *minlev, int *maxlev );
boolean g2mbl_MLCPSuggestNLevels ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte *mkcp,
                                   int *minlev, int *maxlev );

void g2mbl_MLSetNextBlock ( void *data, int nbl );

/* shape only multilevel optimization algorithm */
boolean g2mbl_MLSOptInitd ( int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                            int nhe, BSMhalfedge *mhe,
                            int nfac, BSMfacet *mfac, int *mfhei,
                            byte *mkcp,
                            int nkn1, int nkn2, short nlevels, void **data );
boolean g2mbl_MLSCMPOptInitd ( /* fine mesh */
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei,
                    byte *fmkcp,
                           /* coarse mesh */
                    int cnv,
                           /* refinement matrix */
                    int rmnnz, index2 *rmnzi, double *rmnzc,
                           /* optimization parameters */
                    int nkn1, int nkn2, short nlevels,
                           /* created data structure */
                    void **data );
boolean g2mbl_MLSOptIterd ( void *data, boolean *finished );
boolean g2mbl_MLSCOptIterd ( void *data, boolean *finished );

boolean g2mbl_MLSSuggestNLevels ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  byte *mkcp,
                                  int *minlev, int *maxlev );
boolean g2mbl_MLCPSSuggestNLevels ( int nv, BSMvertex *mv, int *mvhei,
                                    int nhe, BSMhalfedge *mhe,
                                    int nfac, BSMfacet *mfac, int *mfhei,
                                    byte *mkcp,
                                    int *minlev, int *maxlev );

#ifdef G2MBL_TIME_IT
void g2mbl_MLGetTimes ( void *data,
                        float *time_prep, float *time_h, float *time_cg );
#endif
#ifdef __CPLUSPLUS
}
#endif

#endif /*G2MBLENDINGD_H*/

