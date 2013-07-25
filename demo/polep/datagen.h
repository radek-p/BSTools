
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */


#define SCRATCH_MEM_SIZE 1048576

/* displaying switches */
extern boolean netk3, netk5, netk6, netk8;

extern float dataparam[3];
extern int hole_k;
extern int hole_np;

#define MAX_K    8
#define MAX_BPTS (12*MAX_K+2)
extern point3f Bpt[MAX_BPTS];

void InitHole ( int kk, float ip1, float ip2, float ip3 );
void GetBspInd ( int i, int j, int *ind );
