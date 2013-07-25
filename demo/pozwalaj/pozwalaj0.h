
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#define POZWALAJ0_H

#define SCRATCHMEMSIZE (128*1048576)  /* 128 MB */

#define MAX_PTHREADS       32

#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      63
#define MAX_FILENAME_LGT   64

extern char          **prog_argv;
extern int           ncpu;  /* number of CPUs detected */

