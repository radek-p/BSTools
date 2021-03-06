
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

/*#define DEBUG_IPC*/

#define POZWALAJ_IPC_MAGIC      3 /* ought to be unique for each application */

#define IPC_BUFFER_LENGTH      64

/* commands for the child process */
#define ipccmd_GET_DATA         1

/* messages from child to parent process */
#define ipccmd_CHILD_READY      2
#define ipccmd_SUCCESS          3
#define ipccmd_ERROR            4
#define ipccmd_PARTIAL_RESULT   5
#define ipccmd_FINAL_RESULT     6
#define ipccmd_ITER_INFO        7

/* messages from the child to itself */
#define ipccmd_BEGIN_BLP        8
#define ipccmd_CONTINUE_BLP1    9
#define ipccmd_CONTINUE_BLP2   10
#define ipccmd_BEGIN_BSM       11
#define ipccmd_CONTINUE_BSM    12
#define ipccmd_BEGIN_BSCMC     13
#define ipccmd_CONTINUE_BSCMC  14
#define ipccmd_SEND_RESULT     15

/* states of the parent process, related with the communication */
#define ipcstate_NO_CHILD       0
#define ipcstate_CHILD_LAUNCHED 1
#define ipcstate_CHILD_READY    2
#define ipcstate_CHILD_BUSY     3

/* identifiers of data items sent between the parent and child processes */
#define ipcd_PRETRANS          1

  /* blending B-spline patch optimization */
#define ipcd_BLP_SIZE          2
#define ipcd_BLP_CPOINTS       3
#define ipcd_BLP_MKCP          4
#define ipcd_BLP_OPTIMIZE      5

  /* blending mesh surface optimization */
#define ipcd_BSM_SIZE          6
#define ipcd_BSM_VERT          7
#define ipcd_BSM_VHE           8
#define ipcd_BSM_VERTC         9
#define ipcd_BSM_VERTMK       10
#define ipcd_BSM_HALFE        11
#define ipcd_BSM_FAC          12
#define ipcd_BSM_FHE          13
#define ipcd_BSM_COARSE_SIZE  14
#define ipcd_BSM_COARSE_VERT  15
#define ipcd_BSM_COARSE_VHE   16
#define ipcd_BSM_COARSE_HALFE 17
#define ipcd_BSM_COARSE_FAC   18
#define ipcd_BSM_COARSE_FHE   19
#define ipcd_BSM_OPTIMIZE     20

  /* Menger curvature optimization */
#define ipcd_BSC_SIZE         21
#define ipcd_BSC_KNOTS        22
#define ipcd_BSC_CPOINTS      23
#define ipcd_BSC_MKCP         24
#define ipcd_BSC_OPTIMIZE     25
#define ipcd_BSC_MCINFO       26

/* data related to optimization of B-spline blending patches */
typedef struct {
    char gcont, cpdimen, degu, degv, closed_u, clamped;
    int  lastknotu, lastknotv, pitch;
  } ipc_blp_size;

typedef struct {
    int     bl_range[4];
    int     nkn1, nkn2, maxit;
    double  C;
    trans3d pretrans;
  } ipc_blp_options;

/* data related to optimization of blending surfaces represented by a mesh */
typedef struct {
    int nv, nhe, nfac, spdimen, cpdimen, degree;
  } ipc_bsm_size;

typedef struct {
    int     nkn1, nkn2, nlevels, nblocks, maxit, startbl;
    double  C;
    trans3d pretrans;
    boolean use_constraints, shape_only, use_coarse, alt_multilevel;
    byte    constr_mask;
    byte    npthreads;
  } ipc_bsm_options;

/* data related to minimization of integral Menger curvature */
/* of closed B-spline curves */
typedef struct {
    int     spdimen, cpdimen, degree, lkn;
    boolean closed;
  } ipc_bscmc_size;

typedef struct {
    int    nqkn, maxiter, ppopt;
    double exponent, pparam[MENGERC_NPPARAM];
    byte   npthreads;
  } ipc_bscmc_options;

typedef struct {
    double exponent;
    double pparam[MENGERC_NPPARAM];
    double imc, imcp, length;
    double pfv[MENGERC_NPPARAM];
  } ipc_bscmc_info;

/* low-level protocol */
typedef struct {
    int  desc;  /* data item description */
    int  size;  /* size of pointed data in bytes */
    void *ptr;
  } ipc_data_item;


extern ipc_data_item ipc_buffer[IPC_BUFFER_LENGTH];
extern int           ipc_buf_count, ipc_data_size;


#ifdef PARENT_SIDE
extern int ipc_state;

/* parent procedures */
boolean LaunchAChildProcess ( void );
void ProcessChildMessage ( int msg, int size );

void ResetIPCBuffer ( void );
boolean IPCAppendDataItem ( int desc, int size, void *ptr );
void IPCSendData ( void );
void IPCWakeUpChild ( void );
void BindChildToGeomObject ( geom_object *go );
void IPCInterruptTheChild ( void );

void AssignBSCurveData ( GO_BSplineCurve *obj, int size, char *buf );
void AssignBSPatchData ( GO_BSplinePatch *obj, int size, char *buf );
void AssignBSMeshData ( GO_BSplineMesh *obj, int size, char *buf );
#endif

#ifdef CHILD_SIDE
/* child procedures */
void ResetIPCBuffer ( void );
boolean IPCAppendDataItem ( int desc, int size, void *ptr );
void IPCSendData ( void );
#endif

