
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

/*#define SCRATCH_MEM_SIZE (256*1048576)*/  /* 256 MB */
#define SCRATCH_MEM_SIZE (128*1048576)  /* 128 MB */

#define MASK_CP_SPECIAL  0x04  /* same as in editor.h */

extern jmp_buf escape_place;
extern boolean jump_ready;

extern ipc_data_item ipc_buffer[];
extern int           ipc_buf_count, ipc_data_size;

extern void    *optdata;
extern int     itn;
extern struct  tms start, stop;
extern clock_t time0, time1;   
extern boolean finished;

extern trans3d pretrans_inv;

/* B-spline curve with optimized Menger curvature */
extern ipc_bscmc_size    bscmcsize;
extern ipc_bscmc_options bscmcoptions;
extern double            *bscmcknots, *bscmccp;
extern byte              *bscmcmkcp;
extern int               bscmccpsize, bscmcknsize, bscmcmkcpsize;
extern mengerc_data      bscmcdata;
extern ipc_bscmc_info    bscmcinfo;

/* blending B-spline patch with optimized shape */
extern ipc_blp_size      blpsize;
extern ipc_blp_options   blpoptions;
extern double            *blpcp, *_blpcp;
extern char              *blpmkcp;
extern int               blpcpsize, blpmkcpsize;

/* blending mesh surface with optimized shape */
extern ipc_bsm_size      bsmsize;
extern ipc_bsm_options   bsmoptions;
extern BSMvertex         *meshv;
extern int               *meshvhei;
extern double            *meshvpc, *_meshvpc;
extern char              *meshmkcp;
extern BSMhalfedge       *meshhe;
extern BSMfacet          *meshfac;
extern int               *meshfhei;
extern int               meshvpcsize;
/* coarse mesh to define a preconditioner */
extern ipc_bsm_size      cmeshsize;
extern BSMvertex         *cmeshv;
extern int               *cmeshvhei;
extern BSMhalfedge       *cmeshhe;
extern BSMfacet          *cmeshfac;
extern int               *cmeshfhei;
/* refinement matrix - from coarse to fine mesh */
extern int               rmnnz;
extern index2            *rmnzi;
extern double            *rmnzc;

/* pozwalaj_proc01.c */
void ResetIPCBuffer ( void );
boolean IPCAppendDataItem ( int desc, int size, void *ptr );
void IPCSendData ( void );

/* pozwalaj_proc02.c */
void ReadBLPSize ( int isize );
void ReadBLPCPoints ( int isize );
void ReadBLPMkcp ( int isize );
void ReadBLPOptimizeOptions ( int isize );

/* pozwalaj_proc03.c */
void ReadBSMSize ( int isize );
void ReadBSMVert ( int isize );
void ReadBSMVertHE ( int isize );
void ReadBSMVertPC ( int isize );
void ReadBSMVertMK ( int isize );
void ReadBSMHalfedges ( int isize );
void ReadBSMFacets ( int isize );
void ReadBSMFacetHE ( int isize );
void ReadBSMOptimizeOptions ( int isize );

void ReadBSMCSize ( int isize );
void ReadBSMCVert ( int isize );
void ReadBSMCVertHE ( int isize );
void ReadBSMCHalfedges ( int isize );
void ReadBSMCFacets ( int isize );
void ReadBSMCFacetHE ( int isize );

/* pozwalaj_proc04.c */
void GetData ( int size );
boolean InvertPretransformation ( trans3d *pretrans );
void TransformCPoints ( trans3d *tr, int n, int m, int pitch, point3d *cp );

/* pozwalaj_proc05.c */
void BeginBLPOptimization ( void );
void PrepareBLPOutput ( void );
void ContinueBLPOptimization1 ( void );
void ContinueBLPOptimization2 ( void );

/* pozwalaj_proc06.c */
void SetupRefinementMatrix ( void );
void BeginBSMOptimization ( void );
void PrepareBSMOutput ( void );
void MarkBlockPoints ( boolean mark );
void MarkMLBlockPoints ( void );
void ContinueBSMOptimization ( void );

/* pozwalaj_proc07.c */
void ReadBSCSize ( int isize );
void ReadBSCKnots ( int isize );
void ReadBSCCPoints ( int isize );
void ReadBSCMkcp ( int isize );
void ReadBSCMCOptimizeOptions ( int isize );

/* pozwalaj_proc08.c */
void BeginBSCMCOptimization ( void );
void PrepareBSCMCOutput ( mengerc_data *mcd );
void OutputBSCMCInfo ( mengerc_data *mcd );
void ContinueBSCMCOptimization ( void );

/* pozwalaj_proc.c */
void pozwalaj_proc_callback ( int msg, int size );
void my_signal_handler ( void );
#ifdef _FPU_CONTROL_H
void SetupFPEInterrupt ( void );
#endif
void pozwalaj_proc_signal_handler ( int sig );
int main ( int argc, char **argv );

