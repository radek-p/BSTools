
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define POMNIJ_IPC_MAGIC 2 /* ought to be distinct for each application */


#define MAX_BLENDING_CONSTRAINTS 20  /* same as in spl3d.h */

/* commands for the child process */
#define cmdGET_OPTIONS              1
#define cmdGET_BSPATCH              2
#define cmdGET_CONSTRAINTS          3
#define cmdG2BL_OPTIMIZE_LMT        4
#define cmdSEND_BSPATCH             5
#define cmdSEND_CONSTRAINTS         6
#define cmdINTERRUPT                7
#define cmdTERMINATE                8

/* messages from child to parent process */
#define cmdCHILD_READY              9
#define cmdGOT_INTERRUPT           10
#define cmdSUCCESS                 12
#define cmdERROR                   13
#define cmdG2BL_PARTIAL_RESULT     14
#define cmdG2BL_FINAL_RESULT       15
#define cmdG2BL_CONSTRAINTS        16

/* messages from the child to itself */
#define cmdSEND_THE_BSPATCH        17
#define cmdSEND_THE_CONSTRAINTS    18
/* these messages may be sent to the child by itself or by the parent */
#define cmdCONTINUE_OPT            19

/* states of the parent process, related with the communication with the child */
#define ipcSTATE_NOTHING           0
#define ipcSTATE_OPTIONS_SENT      1
#define ipcSTATE_PATCH_SENT        2
#define ipcSTATE_CONSTRAINTS_SENT  3
#define ipcSTATE_G2BLOPT_LAUNCHED  4
#define ipcSTATE_G2BLOPT_FINISHED  5


typedef struct {
    double  NLConst;
    int     maxiter, nkn1, nkn2;
    int     opt_range[4];         /* piece of the patch to optimize */
    int     nconstr;              /* number of constraints */
    char    gcont;                /* 1 or 2 for biquadratic or bicubic patches */
    boolean send_partial;         /* switches */
    boolean closed;
    boolean dumpdata;
    trans3d trans;                /* pre-transformation */
  } ipc_options;

extern int ipc_state;

