
/* ///////////////////////////////////////////////////////////////////////// */  
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */   
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#define COUNT
#define _DEBUG

#define MYINFINITY    1.0e+308

/*#define GRDIV(a,b) ((1.0-TAU)*(a)+TAU*(b)) */
#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

#define CG_THRESHOLD  10000
/*#define CG_THRESHOLD  1000*/


typedef struct {
    char   type;
    short  ncp;       /* number of control points relevant for this element */
    int    firstcpi;  /* index of the first relevant control point number */
                      /* and (x3) to the gradient integrals table */
    int    hti;       /* index to the Hessian integrals table */
    double C;         /* regularization constant */
    int    cfn;       /* central facet number */
    int    blmask;    /* bit mask for blocks, unused by nonblock procedures */
  } meshdom_elem;

typedef struct {
    int   vnum;      /* number of the mesh vertex */
    int   firstsel;  /* number of the first selected neighbour */
    short nneigh;    /* number of nonzero elements in its Hessian row (neighbours) */
    short nncn;      /* number of neighbours not chosen yet */
    int   fneigh;    /* index of the first neighbour element */
  } vertex_desc;

typedef struct {
    int     n0, n1, m0, m1;     /* block limits and influence interval */
    int     ncp, nvars;
    int     ndomelems, *domelind;
    int     *nncpi;
    int     *vpermut, *vncpi;   /* variable permutation and indices */
    int     hsize;              /* band coefficient array size */
    int     *hprof;             /* Hessian block profile */
    double  **hrows, **lhrows;  /* block rows */
    double  func, gn0, nu;
    boolean accurate;
    char    hblstat;            /* 0 - undefined, */
                                /* 1 - computed,  */
                                /* 2 - decomposed (positive), */
                                /* 3 - to be reused */
  } block_desc;

typedef struct {
    int firsthbl, nhbl;
  } nzHbl_rowdesc;

typedef struct {
        /* mesh data */
    int           nv, nhe, nfac;
    BSMvertex     *mv;
    int           *mvhei;
    point3d       *mvcp;
    BSMhalfedge   *mhe;
    BSMfacet      *mfac;
    int           *mfhei;
    char          *mvtag;
        /* coarse mesh - total number of vertices */
    int           cnv;
        /* the refinement matrix */
    int           rmnnz;
    index2        *rmnzi;
    double        *rmnzc;
        /* relevant vertices of the coarse mesh */
    int           nwcp, *wncpi;
        /* the submatrix used by the preconditioner */
    int           pmnnz;
    index3        *pmnzi;
        /* domain element description */
    int           ndomelems, ndomelcpind;
    meshdom_elem  *domelem;
    int           *domelcpind;
    boolean       eltypes[GH_MAX_K-2];  /* which special elements are present? */

    int           hti;
    double        *ftab, *gtab, *htab;
        /* number of control points to be optimized, their indices */
    int           nvcp, nvars;
    int           *nncpi, *vncpi;
        /* Hessian */
    int           hsize, *hprof;
    double        *Hessian, **hrows, *LHessian, **lhrows;
        /* Hessian coarse mesh preconditioner term */
    int           phsize, *phprof;
    double        **phrows;
    int           nnz1, nnz2;
        /* iterations data */
    int           itn, next_entire;
    double        func, gn0, nu[2];
    boolean       newpoint, accurate;
    int           nkn1, nkn2;
    double        *aqcoeff, *bqcoeff;
    double        *aNitabs[GH_MAX_K-2], *aNijtabs[GH_MAX_K-2],
                  *aMijtabs[GH_MAX_K-2], *aJac[GH_MAX_K-2],
                  *bNitabs[GH_MAX_K-2], *bJac[GH_MAX_K-2];
        /* the data below are used only by the block algorithm */
    boolean       use_cg, last_ntn;
          /* 3x3 nonzero block Hessian representation */
    int           Hblsize;  /* number of nonzero Hessian blocks */
    nzHbl_rowdesc *iHbl;    /* indices for subsequent rows */
    int           *cHbl;    /* indices of columns with nonzero blocks in a row */
    double        *Hbl;     /* actual blocks */
          /* block description */
    int           nbl;                       /* number of blocks */
    block_desc    block[G2MBL_MAX_BLOCKS];   /* block description */
    int           *belind;                   /* indices of domain elements */
                                             /* for each block */
    int           ibl;                       /* number of blocks, which did not */
                                             /* satisfy the criterion for */
                                             /* quadrature switching */
    int           nextblock;
    double        *BHessian;

    int           *bltag;              /* for alternative block setup */
#ifdef COUNT
    int           nLM, nN, nF, nG, nH;
#endif
  } mesh_lmt_optdata;


void _g2mbl_TagBoundaryCondVert ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe, int cvn,
                                  char *vtag );
void _g2mbl_CountRegularFacets ( int d, int *vertnum, int *mtab, void *usrptr );
void _g2mbl_CountSpecialVertices ( int d, int k, int *vertnum, int *mtab,
                                   void *usrptr );
void _g2mbl_GetRegularFacetVertNum ( int d, int *vertnum, int *mtab, void *usrptr );
void _g2mbl_GetSpecialElemVertNum ( int d, int k, int *vertnum, int *mtab,
                                    void *usrptr );
boolean _g2mbl_AssignMeshd ( mesh_lmt_optdata *d,
                             int nv, BSMvertex *mv, int *mvhei, point3d *mvcp,
                             int nhe, BSMhalfedge *mhe,
                             int nfac, BSMfacet *mfac, int *mfhei,
                             byte *mkcp );
boolean _g2mbl_OrderCPoints ( int nv, int nvcp, int fcp,
                              int nzcdsize, byte *nzcdistr,
                              int *nncpi, int *vncpi, vertex_desc *nvcpi,
                              int ndomelems, meshdom_elem *domelem, int *domelcpind,
                              int *domelind,
                              int blsize, int *hsize, int *hprof, int *vpermut );
boolean _g2mbl_SetupHessianProfiled ( mesh_lmt_optdata *d, boolean use_blocks );
boolean _g2mbl_AllocBFArraysd ( mesh_lmt_optdata *d, int nkn1, int nkn2 );
boolean _g2mbl_SetupElemConstd ( mesh_lmt_optdata *d,
                                 double dM, double dO, double C );

boolean _g2mbl_SetupHbl3x3d ( mesh_lmt_optdata *d );
boolean _g2mbl_SetupBlockHessianProfiled ( mesh_lmt_optdata *d, int bnum );
boolean _g2mbl_SetupBlocksd ( mesh_lmt_optdata *d, int nbl );
boolean _g2mbl_SetupBlockRowsd ( mesh_lmt_optdata *d );

int g2mbl_NiSize ( int nkn, int k );
int g2mbl_NijSize ( int nkn, int k );
int g2mbl_MijSize ( int nkn, int k );

boolean g2mbl_TabNid ( int hole_k, int nkn, double *qknots,
                       double *Nitab, double *Jac, boolean reparam );
void g2mbl_TabNijd ( int nf, int nkn, double *Nitab, double *Nijtab );
void g2mbl_TabMijd ( int nf, int nkn, double *Nitab, double *Mijtab );

boolean g2mbl_FindDistances ( int nv, BSMvertex *mv, int *mvhei,
                              int nhe, BSMhalfedge *mhe,
                              int nfac, BSMfacet *mfac, int *mfhei,
                              int *fdist );

void _g2mbl_UCompRDerd ( int nkn2, double *Nitab, int knot,
                         int *cpind, point3d *mvcp, vector3d pder[11] );
void g2mbl_UFuncRSQd ( int nkn, double *qcoeff, double *Nitab,
                       int *cpind, point3d *mvcp,
                       double C, double *ftab );
void _g2mbl_UCompSDerd ( int nkn2, double *Nitab, int knot, int k, int l,
                         int *cpind, point3d *mvcp, vector3d pder[11] );
void g2mbl_UFuncSSQd ( int nkn, double *qcoeff, int k, double *Nitab, double *Jac,
                       int *cpind, point3d *mvcp, double C,
                       double *ftab );
void _g2mbl_UCompDGStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           double *DGstar );
void _g2mbl_UCompDBStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           double *DBstar );
void _g2mbl_UFuncGradSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                    int knot, int k, int l, vector3d pder[11],
                                    int *cpind, int *nncpi,
                                    double Gstar[9], double *DGstar,
                                    double Bstar[11], double *DBstar,
                                    double C,
                                    double *F, double *DF );
void g2mbl_UFuncGradRSQd ( int nkn, double *qcoeff, double *Nitab,
                           int *cpind, int *nncpi, point3d *mvcp,
                           double C,
                           double *ftab, double *gtab );
boolean g2mbl_UFuncGradSSQd ( int nkn, double *qcoeff, int k,
                              double *Nitab, double *Jac,
                              int *cpind, int *nncpi, point3d *mvcp,
                              double C,
                              double *ftab, double *gtab );
void _g2mbl_UCompDDGstard ( int nkn2, int nf, double *Nijtab,
                            int knot, int k, int l,
                            int *cpind, int *nncpi, double *DDGstar );
void _g2mbl_UCompDDBstard ( int nkn2, int nf, double *Mijtab, int knot, int k, int l,
                            vector3d *pder, int *cpind, int *nncpi,
                            double *DDBstar );
boolean _g2mbl_UFuncGradHessSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                int knot, int k, int l, vector3d *pder,
                                int *cpind, int *nncpi,
                                double *Gstar, double *DGstar, double *DDGstar,
                                double *Bstar, double *DBstar, double *DDBstar,
                                double C,
                                double *F, double *DF, double *DDF );
boolean g2mbl_UFuncGradHessRSQd ( int nkn, double *qcoeff,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  int *cpind, int *nncpi, point3d *mvcp,
                                  double C,
                                  double *ftab, double *gtab, double *htab );
boolean g2mbl_UFuncGradHessSSQd ( int nkn, double *qcoeff, int k,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  double *Jac,
                                  int *cpind, int *nncpi, point3d *mvcp,
                                  double C,
                                  double *ftab, double *gtab, double *htab );

void g2mbl_SFuncGradRSQd ( int nkn, double *qcoeff, double *Nitab,
                           int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                           double *ftab, double *gtab );
boolean g2mbl_SFuncGradSSQd ( int nkn, double *qcoeff, int k,
                              double *Nitab, double *Jac,
                              int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                              double *ftab, double *gtab );
boolean g2mbl_SFuncGradHessRSQd ( int nkn, double *qcoeff,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                                  double *ftab, double *gtab, double *htab );
boolean g2mbl_SFuncGradHessSSQd ( int nkn, double *qcoeff, int k,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  double *Jac,
                                  int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                                  double *ftab, double *gtab, double *htab );

/* nonblock algorithm functions */
double g2mbl_UFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                      int nv, point3d *mvcp,
                      int ndomelems, meshdom_elem *domelem, int *domelcpind,
                      boolean force, double *ftab );
boolean g2mbl_UFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                        int nv, point3d *mvcp,
                        int ndomelems, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi,
                        boolean force, double *ftab, double *gtab,
                        double *func, int nvars, double *grad );
boolean g2mbl_UFuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac,
                        int nv, point3d *mvcp,
                        int ndomelems, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int hsize, int *hprof, double **hrows );

void _g2mbl_AddCPIncrement ( int nvcp, int *vncpi, point3d *mvcp,
                             vector3d *incr, point3d *omvcp );
double _g2mbl_AuxNuFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int nvcp, int *vncpi,
                           int ndomelems, meshdom_elem *domelem, int *domelcpind,
                           double *ftab, int nvars, int hsize, int *hprof,
                           double *Hessian, double *LHessian, double **lhrows,
                           double *grad, double *dcoeff, point3d *auxmvcp,
                           double nu );

/* block algorithm functions */
boolean g2mbl_B3x3FuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac,
                        int nv, point3d *mvcp,
                        int ndomelems, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int Hblsize, nzHbl_rowdesc *iHbl, int *cHbl, double *Hbl );
double g2mbl_AFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                      int nv, point3d *mvcp, int bdomelems, int *domelind,
                      meshdom_elem *domelem, int *domelcpind,
                      boolean force, double *ftab );
boolean g2mbl_AFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int bdomelems, int *domelind,
                           meshdom_elem *domelem, int *domelcpind, int *bncpi,
                           boolean force, double *ftab, double *gtab,
                           double *func, int nvars, double *grad );
boolean g2mbl_AFuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac, int nv, point3d *mvcp,
                        int bdomelems, int *domelind, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi, int *bncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int hsize, int *hprof, double **hrows );

double _g2mbl_AuxAltBNuFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int nvcp, int *vncpi,
                           int ndomelems, int *domelind, meshdom_elem *domelem,
                           int *domelcpind, double *ftab,
                           int nvars, int hsize, int *hprof,
                           double *Hessian, double *LHessian, double **lhrows,
                           double *grad, double *dcoeff, point3d *auxmvcp,
                           double nu );

boolean _g2mbl_SetupAltBlocks ( int nv, BSMvertex *mv, int *mvhei,
                                int nhe, BSMhalfedge *mhe,
                                int nfac, BSMfacet *mfac, int *mfhei,
                                int *nncpi, char nbl, int *bltag, int *blseed );
boolean _g2mbl_SetupAltBlockHessianProfiled ( mesh_lmt_optdata *d, int bnum );
boolean _g2mbl_SetupAltBlockDescription ( mesh_lmt_optdata *d, int nbl );

boolean _g2mbl_MultAxd ( int n, void *usrdata, const double *x, double *Ax );
boolean _g2mbl_MultQixAltd ( int n, void *usrdata, const double *x, double *Qix );

boolean _g2mbl_CMPSetupCGPrecondd ( mesh_lmt_optdata *d );

/* auxiliary procedures */
void _g2mbl_OutputNZDistr ( mesh_lmt_optdata *d );

