
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define SCRATCH_MEM_SIZE 268435456 /* 256MB */

extern jmp_buf escape_place;
extern boolean jump_ready;
extern boolean got_interrupted;

extern ipc_options options;

extern int degu, degv, lknu, lknv, dim;
extern double *knotsu, *knotsv, *cpoints;
extern boolean clu, clv;  /* clamped boundary indicators */       

  /* constraints data: curves of constant "u" parameter */
extern int nucknots;
extern double *ucknots, *ucurvcp;
/* constraint equations matrix and right hand side */
extern double *cmat, *crhs;

extern int     ncp, optlknu, optlknv, fcp;
extern point3d *ccp, *occp;
extern void    *optdata;
extern int     itn;
extern boolean finished;
extern struct tms start;


void DumpData ( void );
void WriteThePatch ( int spdimen, int udeg, int lastknotu, double *knotsu,
                     int vdeg, int lastknotv, double *knotsv,
                     int pitch, double *cpoints );
void WriteTheConstraints ( int spdimen, int nucknots, double *ucknots,           
                           int vdeg, int lastknotv, double *knotsv, 
                           int pitch, double *cp );

boolean InputBSPatch ( int size );
boolean InputOptions ( int size );
boolean InputConstraints ( int size );

int PatchDataSize ( void );
void OutputBSPatch ( void );
int ConstraintsSize ( void );
void OutputConstraints ( void );

void SetupG1BLKnots ( void );
void G1FreeBoundaryToClamped ( point3d *ccp );
void G1ClampedBoundaryToFree ( point3d *ccp );
boolean G1AdjustOptRange ( void );

void SetupG2BLKnots ( void );
void G2FreeBoundaryToClamped ( point3d *ccp );
void G2ClampedBoundaryToFree ( point3d *ccp );
boolean G2AdjustOptRange ( void );

void PreTransformation ( int npcp, point3d *pcp );
boolean PostTransformation ( int npcp, point3d *pcp );

boolean G1OptimizeLMTInit ( void );
boolean G1OptimizeLMTIter ( void );
boolean G2OptimizeLMTInit ( void );
boolean G2OptimizeLMTIter ( void );

void pomnij_proc_callback ( int msg, int size );
void my_signal_handler ( void );
void pomnij_proc_signal_handler ( int sig );
void lib_error_handler ( int module, const char *file, int line,
                         int errcode, const char *errstr );
#ifdef _FPU_CONTROL_H
void SetupFPEInterrupt ( void );
#endif

int main ( int argc, char **argv );

