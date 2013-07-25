
/* ///////////////////////////////////////////////////////////////////////// */  
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */   
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* header file for the procedures implementing the multilevel optimization */
/* algorithm */

#ifndef _CONST
#define _CONST const
#endif


/* a divided block must have at least that many variable control points */
#define MIN_NVCP    100
/* a block with at most that many variable control points will use      */
/* band matrix Hessian representation (Cholesky, not CG). Larger blocks */
/* will use band matrix Hessian if they have not been divided.          */
#define MAX_NVCP   3500
/*#define MAX_NVCP   50*/
/* this is an analoguous threshold for the shape only optimization procedure */
#define MAX_NVCPS 10000
/*#define MAX_NVCPS 100*/

/* a subblock of a block with n vertices will have                         */
/* SUBBLOCK_STEP/(2*SUBBLOCK_STEP-1) vertices. At first it was fixed at 3, */
/* experiments might be used to fine-tune it. Perhaps different values are */
/* proper for the methods with and without the coarse mesh preconditioner. */
#define SUBBLOCK_STEP    4
#define SUBBLOCK_STEP_CP 5

#define MYINFINITY    1.0e+308
#define EPS1          1.0e-4
#define EPS2          1.0e-7
#define EPS3          1.0e-8
#define EPS4          1.0e-8
#define DELTA4        1.0e-7
#define BETA          0.125
#define CG_EPS        1.0e-3
#define CG_DELTA      1.0e-20
#define MAXCGN      512
#define MAXNTN       16
#define MAXBTN       20
#define MAXCTN       24
#define MAXGTN       10 /*16*/
#define GRTHR         0.02

#define FLAG_F       0x0001
#define FLAG_G       0x0002
#define FLAG_H       0x0004
#define FLAG_LH      0x0008
#define FLAG_CH      0x0010
#define FLAG_CMH     0x0020
#define FLAG_CMLH    0x0040
#define FLAG_ADVANCE 0x0080

/*#define GRDIV(a,b) ((1.0-TAU)*(a)+TAU*(b)) */
#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

typedef struct {
    int           nvcp, nvars, ndomel, seed;
    int           *vncpi, *domelind;
        /* coarse mesh vertices */
    int           nwcp, *wncpi;
        /* refinement submatrix */
    int           rmnnz;
    index3        *rmnzi;
        /* data for fast matrix multiplications related with */
        /* the coarse mesh preconditioner */
    int           nnz1, nnz2;
        /* Hessian representation for a non-bottom block */
    int           nHbl;  /* total number of nonzero 3x3 or 1x1 blocks */
    nzHbl_rowdesc *iHbl;
    int           *cHbl,
                  *tHbl;
        /* Hessian representation for a bottom block or a sparse grid block */
    int           *hprof, hsize;
    double        **hrows, **lhrows;
        /* function, gradient and Hessian valid flags */
    short         fghflag;
  } mlblock_desc;

typedef struct {
        /* mesh data */
    int          nv, nhe, nfac;
    BSMvertex    *mv;
    int          *mvhei;
    point3d      *mvcp;
    BSMhalfedge  *mhe;
    BSMfacet     *mfac;
    int          *mfhei;
    char         *mvtag;

    int          cnv; /* number of the coarse mesh vertices */

        /* refinement matrix - for the optional preconditioner */
        /* obtained with a coarse mesh */
    int          rmnnz;      /* number of nonzero coefficients */
    index2       *rmnzi;     /* their positions */
    double       *rmnzc;     /* in this array */  

        /* normal vectors associated with the control points, */
        /* used only by the shape-only optimization procedure */
    vector3d     *mvcpn;

        /* domain element description */
    int          tdomelems,   /* total number of elements */
                 ndomelems,   /* number of relevant elements */
                 ndomelcpind; /* number of vertex indices for the elements */
    meshdom_elem *domelem;
    int          *domelcpind;
    boolean      eltypes[GH_MAX_K-2];  /* which special elements are present? */

        /* quadrature knots and values of basis functions and their derivatives */
    int          nkn1, nkn2;
    double       *aqcoeff, *bqcoeff;
    double       *aNitabs[GH_MAX_K-2], *aNijtabs[GH_MAX_K-2],
                 *aMijtabs[GH_MAX_K-2], *aJac[GH_MAX_K-2],
                 *bNitabs[GH_MAX_K-2], *bJac[GH_MAX_K-2];

        /* tables of integrals over the elements */
    int          hti, bls;
    double       *ftab1, *ftab2, *gtab1, *gtab2, *htab1;

        /* 3x3 or 1x1 nonzero blocks of the Hessian matrix */
    int          Hblsize;  /* number of nonzero Hessian blocks */
    double       *Hbl;     /* actual 3x3 or 1x1 blocks */

        /* number of control points to be optimized */
    int          nvcp, nvars;
    int          *nncpi;  /* numbers of variable vertices */
                          /* (-1 for fixed vertices) */

        /* blocks */
    short        nlevels, nblocks;
    mlblock_desc *bd;

        /* iteration specific data */
    short        currentblock, dirtyblock, lastblock, nextlevel;
    double       nu[2];

        /* diagnostic ouput */
    short        log_level;   /* printf: 0 - none, 1 - iterations, 2 - values */
    short        error_code;  /* 0 - none */
#ifdef G2MBL_TIME_IT
    int          time_prep, time_h, time_cg;
#endif
  } mesh_ml_optdata;


typedef struct {
    mesh_ml_optdata *d;
    short           bl, cbl;
    double          nu;
    boolean         failure;
  } mesh_ml_cg_data;

typedef boolean (*cg_mult)(int nvars, void *usrdata, const double *x, double *Ax );


boolean _g2mbl_MLDivideBlock ( int nv, BSMvertex *mv, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               int nvcp, int *vncpi,
                               int *nvbcp1, int **vnbcp1, int *seed1,
                               int *nvbcp2, int **vnbcp2, int *seed2,
                               int step );
boolean _g2mbl_MLAssignMeshd ( mesh_ml_optdata *d,
                               int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               byte *mkcp );
boolean _g2mbl_MLSetupBlocksd ( mesh_ml_optdata *d,
                                short nlevels, short bsize, int step );
boolean _g2mbl_MLOptAllocBFArraysd ( mesh_ml_optdata *d, int nkn1, int nkn2,
                                     boolean reparam );
boolean _g2mbl_MLFindElementsd ( mesh_ml_optdata *d, int bls, boolean dnkn );
boolean _g2mbl_MLFindBlockElementsd ( mesh_ml_optdata *d );
boolean _g2mbl_MLFindVCPNeighboursd ( mesh_ml_optdata *d,       
                                      nzHbl_rowdesc **vn, int **vni );
boolean _g2mbl_MLSetupBlockCGHessiand ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_MLSetupBlockCholHessiand ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_MLSetupBlockHessiansd ( mesh_ml_optdata *d );
boolean _g2mbl_MLSetupElemConstd ( mesh_ml_optdata *d,
                                   double dM, double dO, double C );

double g2mbl_MLFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                       int nv, point3d *mvcp, int ndomel, int *domelind,
                       meshdom_elem *domelem, int *domelcpind,
                       boolean recalc, double *ftab );
boolean g2mbl_MLFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                            int nv, point3d *mvcp, int ndomel, int *domelind,
                            meshdom_elem *domelem, int *domelcpind,
                            int nvcp, int *vncpi,
                            boolean recalc, double *ftab, double *gtab,
                            double *func, double *grad );
boolean g2mbl_MLFuncGradHessianAd ( int nkn, double *qcoeff,
             double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
             int nv, point3d *mvcp, int ndomel, int *domelind,
             meshdom_elem *domelem, int *domelcpind,
             int nvcp, int *vncpi,
             boolean recalc, double *ftab, double *gtab, double *htab,
             double *func, double *grad,
             int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl );
boolean g2mbl_MLFuncGradHessianBd ( int nkn, double *qcoeff,
             double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
             int nv, point3d *mvcp, int ndomel, int *domelind,
             meshdom_elem *domelem, int *domelcpind,
             int nvcp, int *vncpi,
             boolean recalc, double *ftab, double *gtab, double *htab,
             double *func, double *grad,
             int hsize, int *hprof, double **hrows );
boolean g2mbl_MLGetHessianRowsd ( int nv, int nvcp1, int *vncpi1,
                                  int nHbl, nzHbl_rowdesc *iHbl,
                                  int *cHbl, int *tHbl, double *Hbl,
                                  int nvcp, int *vncpi,
                                  int hsize, int *hprof, double **hrows );

boolean _g2mbl_MLmultAxd ( int nvars, void *usrdata,
                           _CONST double *x, double *Ax );
boolean _g2mbl_MLmultQIxd ( int nvars, void *usrdata,
                            _CONST double *x, double *Qix );
boolean _g2mbl_MLDecomposeBlockPrecond ( mesh_ml_optdata *d, short int bl,
                                         double nu, boolean *positive );
boolean g2mbl_MLOptBlockAd ( void *data, int bl );
boolean g2mbl_MLOptBlockBd ( void *data, int bl );

/* ///////////////////////////////////////////////////////////////////////// */
/* procedures used by the multilevel algorithm with the preconditioner using */
/* the coarse mesh */
boolean _g2mbl_CMPAssignMeshd ( mesh_ml_optdata *d,
                    int fnv, BSMvertex *fmv, int *fmvhei, point3d *fmvcp,
                    int fnhe, BSMhalfedge *fmhe,
                    int fnfac, BSMfacet *fmfac, int *fmfhei, byte *fmkcp,
                    int cnv, int rmnnz, index2 *rmnzi, double *rmnzc );
boolean _g2mbl_CMPFindCoarseMeshBlock (
                     int fnv, int cnv, int rmnnz, index2 *rmnzi,
                     int nvbcp, int *vnbcp,
                     int *nwbcp, int **wnbcp );
boolean _g2mbl_CMPFindBlockCGPmatrix (
                             int fnv, int cnv, int rmnnz, index2 *rmnzi,
                             int nvbcp, int *vnbcp, int nwbcp, int *wnbcp,
                             int *smnnz, index3 **smi );
boolean _g2mbl_CMPOrderCoarsePoints ( int nwcp, int nnz, index2 *nzci,
                                      int brmnnz, index3 *brmnzi,
                                      int blsize, int *hsize, int *hprof );
boolean _g2mbl_CMPSetupBlockCGPrecondd ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_CMPSetupBlockHessiansd ( mesh_ml_optdata *d );

boolean _g2mbl_CMPMultRTHR3x3d ( int nrowsa, int nnza, index3 *ai, double *ac,
                                 int ncolsb, int nnzb, index3 *bi, double *bc,
                                 int nnz1,
                                 double nu,
                                 int hsize, int *hprof, double **hrows );
boolean _g2mbl_CMPSetupCoarseHessiand ( mesh_ml_optdata *d, int bl, double nu );

/* ///////////////////////////////////////////////////////////////////////// */
/* procedures used by the shape-only optimization algorithm */
boolean g2mbl_MLSFuncGradd ( int nkn, double *qcoeff,
              double **Nitabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab,
              double *func, double *grad );
boolean g2mbl_MLSFuncGradHessianAd ( int nkn, double *qcoeff,
              double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab, double *htab,
              double *func, double *grad,
              int nHbl, nzHbl_rowdesc *iHbl, int *cHbl, int *tHbl, double *Hbl );
boolean g2mbl_MLSFuncGradHessianBd ( int nkn, double *qcoeff,
              double **Nitabs, double **Nijtabs, double **Mijtabs, double **Jac,
              int nv, point3d *mvcp, vector3d *mvcpn,
              int ndomel, int *domelind, meshdom_elem *domelem, int *domelcpind,
              int nvcp, int *vncpi,
              boolean recalc, double *ftab, double *gtab, double *htab,
              double *func, double *grad,
              int hsize, int *hprof, double **hrows );

boolean g2mbl_MLSGetHessianRowsd ( int nv, int nvcp1, int *vncpi1,
                                   int nHbl, nzHbl_rowdesc *iHbl,
                                   int *cHbl, int *tHbl, double *Hbl,
                                   int nvcp, int *vncpi,
                                   int hsize, int *hprof, double **hrows );

boolean _g2mbl_MLSSetupBlockCGHessiand ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_MLSSetupBlockCholHessiand ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_MLSSetupBlockHessiansd ( mesh_ml_optdata *d );
boolean _g2mbl_MLSFindCPNormalsd ( mesh_ml_optdata *d );
boolean _g2mbl_CMPSSetupCoarseHessiand ( mesh_ml_optdata *d, int bl, double nu );
void _g2mbl_MLInvalSmallBlocks ( mesh_ml_optdata *d, int bl );
boolean _g2mbl_MLSmultAxd ( int nvars, void *usrdata, double *x, double *Ax );
boolean _g2mbl_MLSmultQIxd ( int nvars, void *usrdata, double *x, double *Qix );
boolean _g2mbl_MLSDecomposeBlockPrecond ( mesh_ml_optdata *d, short int bl,
                                          double nu, boolean *positive );
void _g2mbl_MLSAddCPIncrement ( int nvcp, int *vncpi,
                                point3d *mvcp, vector3d *mvcpn,
                                double *incr, point3d *omvcp );
boolean g2mbl_MLSOptBlockAd ( void *data, int bl );
boolean g2mbl_MLSOptBlockBd ( void *data, int bl );

boolean _g2mbl_CMPSSetupBlockCGPrecondd ( mesh_ml_optdata *d, int bn );
boolean _g2mbl_CMPSSetupBlockHessiansd ( mesh_ml_optdata *d );

int _g2mbl_MLNextBlockNumd ( mesh_ml_optdata *d, boolean advance );

