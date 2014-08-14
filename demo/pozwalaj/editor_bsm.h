
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

  /* maximal degree of a mesh facet or vertex */
#define MAX_BSM_DEGREE     16
#define MAX_BSM_NV     131072
#define MAX_BSM_NHE    524288
#define MAX_BSM_NFAC   131072

/* bit masks for mesh edges to be drawn */
#define MASK_HE_INNER    0x01
#define MASK_HE_BOUNDARY 0x02

/* number of display lists */
#define BSM_NDL       5

/* display list masks and numbers */
#define BSM_DLM_CNET     DLISTMASK_CPOINTS /*0x0001*/
#define BSM_DL_CNET      0
#define BSM_DLM_SURF     0x0002
#define BSM_DL_SURF      1
#define BSM_DLM_ELEM     0x0004
#define BSM_DL_ELEM      2
#define BSM_DLM_SPECIAL  0x0008
#define BSM_DL_SPECIAL   3
#define BSM_DLM_HOLEFILL 0x0010
#define BSM_DL_HOLEFILL  4

#define BSM_NCURRENT_VERTICES 2
#define BSM_NCURRENT_EDGES    2
#define BSM_NCURRENT_FACETS   2

/* blending surface mesh optimization parameter range */
#define BSM_BL_MIN_CPARAM  0.0001
#define BSM_BL_MAX_CPARAM  0.1
#define BSM_BL_MIN_NKN     3
#define BSM_BL_MAX_NKN    10
#define BSM_BL_MAX_ITER   50
#define BSM_BL_MIN_CFACT   0.001
#define BSM_BL_MAX_CFACT   1.0
#define BSM_BL_MAX_SPRAD  20 


typedef struct {
    geom_object me;
    int         maxnv, maxnhe, maxnfac;
    int         degree, nv, nhe, nfac;
    BSMvertex   *meshv;
    int         *meshvhei;
    double      *meshvpc;
    BSMhalfedge *meshhe;
    BSMfacet    *meshfac;
    int         *meshfhei;
    int         savedsize;
    double      *savedcpoints;
    byte        *mkcp;
    boolean     rational, subdivision;
    boolean     view_surf, view_cnet, view_special, view_holefill;
    boolean     integrity_ok;
    boolean     fill_Coons, fill_Bezier, fill_G1, fill_G2, fill_G1Q2;
    bsm_special_elem_list
                spvlist;
    double      sl_g1q2param, g1q2param;
    char        density, subdivl, blending;
    boolean     bl_constr, bl_shape_only, bl_use_coarse, spvlist_ok;
    char        coarse_name[MAX_NAME_LENGTH+1];
    int         current_vertex[BSM_NCURRENT_VERTICES],
                current_edge[BSM_NCURRENT_EDGES],
                current_facet[BSM_NCURRENT_FACETS];
    double      bsm_bl_C;
    int         nkn1, nkn2, nlevels, startfrom, nblocks, maxit;
        /* special patches of a cubic blending surface in R^3, */
        /* which may be optimised separately */
    int         nspecial_patches, special_deg;
    double      *special_patches;
    boolean    special_patches_ok;
  } GO_BSplineMesh;


/* B-spline meshes processing procedures */
void GeomObjectAssignBSplineMesh ( GO_BSplineMesh *obj,
                      int spdimen, boolean rational,
                      int nv, BSMvertex *mv, int *mvhei, double *mvpc,
                      int nhe, BSMhalfedge *mhe,
                      int nfac, BSMfacet *mfac, int *mfhei, byte *mkcp );
boolean GeomObjectExtendBSplineMesh ( GO_BSplineMesh *obj,
                      int spdimen, boolean rational,
                      int nv, BSMvertex *mv, int *mvhei, double *mvpc,
                      int nhe, BSMhalfedge *mhe,
                      int nfac, BSMfacet *mfac, int *mfhei, byte *mkcp );
boolean GeomObjectInitBSplineMesh ( GO_BSplineMesh *obj,
                                    char spdimen, boolean rational );
geom_object *GeomObjectAddBSplineMesh ( const char *name,
                                        char spdimen, boolean rational );
geom_object *GeomObjectCopyBSplineMesh ( GO_BSplineMesh *obj );
void GeomObjectDeleteBSplineMesh ( GO_BSplineMesh *obj );

void GeomObjectDrawBSplineSubdivMesh ( GO_BSplineMesh *obj );
void GeomObjectDrawMeshPatch ( int d, int *vertnum, int *mtab, void *usrptr );
void GeomObjectDrawBSplineMeshPatches ( GO_BSplineMesh *obj );
void GeomObjectDrawBSplineMeshSurf ( GO_BSplineMesh *obj );
void GeomObjectDrawBSplineMNet ( GO_BSplineMesh *obj );
void GeomObjectDrawBSplineMElements ( GO_BSplineMesh *obj );
void GeomObjectBSMDrawVSpecial ( int d, int k, int *vertnum, int *mtab,
                                 void *usrptr );
void GeomObjectBSMDrawFSpecial ( int d, int k, int *vertnum, int *mtab,
                                 void *usrptr );
void GeomObjectDrawBSplineMSpecialElements ( GO_BSplineMesh *obj );
void GeomObjectDisplayBSplineMesh ( GO_BSplineMesh *obj );
void GeomObjectBSplineMeshOutputBezierPatch ( int n, int m, const double *cp,
                                              void *usrptr );

boolean GeomObjectBSplineMeshFindSpecialVertices ( GO_BSplineMesh *obj, int d );

void GeomObjectDrawBSplineMeshCubicHoleFillingPatches ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshSetDensLevel ( GO_BSplineMesh *obj, int dens );
boolean GeomObjectBSplineMeshSetSubdiv ( GO_BSplineMesh *obj, boolean on );

boolean GeomObjectBSplineMeshSetCurrentVertex ( GO_BSplineMesh *obj, int vn, int cvi );
boolean GeomObjectBSplineMeshSetCurrentEdge ( GO_BSplineMesh *obj, int en, int cei );
boolean GeomObjectBSplineMeshSetCurrentFacet ( GO_BSplineMesh *obj, int fn, int cfi );
boolean GeomObjectBSplineMeshSetDegree ( GO_BSplineMesh *obj, int deg );
boolean GeomObjectBSplineMeshSetBlending ( GO_BSplineMesh *obj, boolean bl );
boolean GeomObjectBSplineMeshGetCurrentVertices ( GO_BSplineMesh *obj,
                                                  int n, point3d *p );

void GeomObjectBSplineMeshFindBBox ( GO_BSplineMesh *obj, Box3d *box );
boolean GeomObjectBSplineMeshFindCPoint ( GO_BSplineMesh *obj,
                                          CameraRecd *CPos, short x, short y,
                                          int *dist );
void GeomObjectBSplineMeshSetCPoint ( GO_BSplineMesh *obj,
                                      CameraRecd *CPos, short x, short y );
void GeomObjectBSplineMeshMarkCPoints ( GO_BSplineMesh *obj,
                                        CameraRecd *CPos, Box2s *box,
                                        byte mask, int action );
void GeomObjectBSplineMeshMarkCPoint ( GO_BSplineMesh *obj,
                                       byte mask, int action );
boolean GeomObjectBSplineMeshSaveCPoints ( GO_BSplineMesh *obj );
void GeomObjectBSplineMeshUndoLastTransformation ( GO_BSplineMesh *obj );
void GeomObjectBSplineMeshTransformCPoints ( GO_BSplineMesh *obj,
                                             trans3d *tr, byte mask );

boolean GeomObjectBSplineMeshGetPointCoord ( GO_BSplineMesh *obj, int p,
                                  int *spdimen, int *cpdimen, double **pc );

boolean GeomObjectBSplineMeshRefinement ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshDoubling ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshAveraging ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshExtractSubmesh ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshRemoveCurrentVertex ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshContractCurrentEdge ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshShrinkCurrentEdge ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshGlueEdges ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshGlueEdgeLoops ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshSealHole ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshSplitBoundaryEdge ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshRemoveCurrentFacet ( GO_BSplineMesh *obj );
boolean GeomObjectBSplineMeshDoubleCurrentFacEdges ( GO_BSplineMesh *obj );

boolean GeomObjectBSplineMeshInitTetrahedron ( GO_BSplineMesh *obj, boolean add );
boolean GeomObjectBSplineMeshInitCube ( GO_BSplineMesh *obj, boolean add );
boolean GeomObjectBSplineMeshInitDodecahedron ( GO_BSplineMesh *obj, boolean add );
boolean GeomObjectBSplineMeshInitKGon ( GO_BSplineMesh *obj, int k, boolean add );
boolean GeomObjectBSplineMeshInitKPrism ( GO_BSplineMesh *obj, int k, boolean add );

void GeomObjectBSplineMeshSetHoleFillingOpt ( GO_BSplineMesh *obj,
                        boolean fillCoons, boolean fillBezier,
                        boolean fillG1, boolean fillG2, boolean fillG1Q2 );
void GeomObjectBSplineMeshSetG1Q2param ( GO_BSplineMesh *obj, double param );
GHoleDomaind *GeomObjectBSplineMeshSetupBicubicHoleDomain ( int hole_k );
boolean GeomObjectBSplineMeshFillBicubicHole ( GO_BSplineMesh *obj,
                          int hole_k, int *vertnum,
                          void (*outpatch)( int n, int m, const double *cp,
                                            void *usrptr ) );

boolean GeomObjectWriteBSMAttributes ( GO_BSplineMesh *obj );
boolean GeomObjectWriteBSplineMesh ( GO_BSplineMesh *obj );
boolean GeomObjectBSMResolveDependencies ( GO_BSplineMesh *obj );
void GeomObjectReadBSplineMesh ( void *usrdata,
                    const char *name, int ident, int degree,
                    int nv, const BSMvertex *mv, const int *mvhei,
                    const point4d *vc,
                    int nhe, const BSMhalfedge *mhe,
                    int nfac, const BSMfacet *mfac, const int *mfhei,
                    int spdimen, boolean rational );

void GeomObjectOutputMeshBSPatchToRenderer3D ( int d, int *vertnum, int *mtab,
                                               void *usrptr );
void GeomObjectBSplineMeshOutputBezPatchToRenderer3D ( int n, int m, const double *cp,
                                                       void *usrptr );
void GeomObjectBSplineMeshOutputHoleFillingToRenderer3D ( int d, int k,
               int *vertnum, int *mtab, void *usrptr );
void GeomObjectBSplineMeshOutputToRenderer3D ( GO_BSplineMesh *obj );

boolean GeomObjectBSplineMeshSelectPoint ( GO_BSplineMesh *obj, CameraRecd *CPos,
                                           short x, short y );
boolean GeomObjectBSplineMeshSelectEdge ( GO_BSplineMesh *obj, CameraRecd *CPos,
                                          short x, short y );
void GeomObjectBSplineMeshUnselectPoint ( GO_BSplineMesh *obj );
void GeomObjectBSplineMeshUnselectEdge ( GO_BSplineMesh *obj );

boolean GeomObjectBSplineMeshMarkBetweenVertices ( GO_BSplineMesh *obj,
                                                   boolean mark );
boolean GeomObjectBSplineMeshFilterPolyline ( GO_BSplineMesh *obj );

void GeomObjectBSplineMeshDisplayInfoText ( GO_BSplineMesh *obj );

void GeomObjectBSplineMeshOptSpecialPatches ( GO_BSplineMesh *obj );

boolean GeomObjectBSplineMeshProcessDep ( GO_BSplineMesh *obj, geom_object *go );
void GeomObjectBSplineMeshProcessDeletedDep ( GO_BSplineMesh *obj, geom_object *go );

