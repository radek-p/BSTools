
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpkvaria library of C procedures -              */
/* miscellaneous but useful stuff                                        */

#ifndef PKVARIA_H
#define PKVARIA_H

#ifndef _SYS_TIMES_H
#include <sys/times.h>
#endif

#ifndef MSGPOOL_H
#include "msgpool.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef false
#define false 0
#endif
#ifndef true
#define true  1
#endif

#define MAX_1BYTE 0xFF
#define MAX_2BYTE 0xFFFF
#define MAX_3BYTE 0xFFFFFF
#define MAX_4BYTE 0xFFFFFFFF

#define EXP1   2.7182818284590452353
#define PI     3.1415926535897932384
#define SQRT2  1.4142135623730950488
#define SQRT3  1.7320508075688772935
#define TAU    0.6180339887498948482  /* golden ratio factor */

#ifndef min
#define min(a,b) ((a)<(b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a)>(b) ? (a) : (b))
#endif

typedef unsigned char boolean;
typedef unsigned char byte;

/* point and vector types, for general use */
typedef struct { signed char x, y; }     char2,   point2c, vector2c;
typedef struct { signed short x, y; }    short2,  point2s, vector2s;
typedef struct { signed short x, y, z; } short3,  point3s, vector3s;
typedef struct { signed int x, y; }      int2,    point2i, vector2i;
typedef struct { signed int x, y, z; }   int3,    point3i, vector3i;
typedef struct { float x, y; }           float2,  point2f, vector2f;
typedef struct { double x, y; }          double2, point2d, vector2d;
typedef struct { float x, y, z; }        float3,  point3f, vector3f;
typedef struct { double x, y, z; }       double3, point3d, vector3d;
typedef struct { float x, y, z, w; }     float4,  point4f, vector4f;
typedef struct { double x, y, z, w; }    double4, point4d, vector4d;


/* 2d boxes, for various purposes */
typedef struct Box2i {
    int  x0, x1, y0, y1;
  } Box2i;

typedef struct Box2s {
    short x0, x1, y0, y1;
  } Box2s;


/* error message handlers */
void pkv_SignalError ( int module, const char *file, int line,
                       int errcode, const char *errstr );
#define PKV_SIGNALERROR(module,errcode,errstr) \
  pkv_SignalError ( module, __FILE__, __LINE__, errcode, errstr )

void pkv_SetErrorHandler (
        void (*ehandler)( int module, const char *file, int line,
                          int errcode, const char *errstr ) );


/* scratch memory management procedures */
boolean pkv_InitScratchMem ( size_t size );
void pkv_DestroyScratchMem ( void );
void PrintScratchMemData ( void );

/* pointers to procedures, which by default point to simple  */
/* allocation/deallocation routines in a scratch memory pool */
/* which is a simple stack; in an application using threads  */
/* running in parallel, alternative routines (using          */
/* a separate stack for each thread) have to be assigned     */
extern void *(*pkv_GetScratchMem) ( size_t size );
extern void (*pkv_FreeScratchMem) ( size_t size );
extern void *(*pkv_GetScratchMemTop) ( void );
extern void (*pkv_SetScratchMemTop) ( void *p );
extern size_t (*pkv_ScratchMemAvail) ( void );
extern size_t (*pkv_MaxScratchTaken) ( void );

#define pkv_GetScratchMemi(size) \
  (int*)pkv_GetScratchMem ( (size)*sizeof(int) )
#define pkv_FreeScratchMemi(size) \
  pkv_FreeScratchMem ( (size)*sizeof(int) )
#define pkv_ScratchMemAvaili() \
  (pkv_ScratchMemAvail()/sizeof(int))

#define pkv_GetScratchMemf(size) \
  (float*)pkv_GetScratchMem ( (size)*sizeof(float) )
#define pkv_FreeScratchMemf(size) \
  pkv_FreeScratchMem ( (size)*sizeof(float) )
#define pkv_ScratchMemAvailf() \
  (pkv_ScratchMemAvail()/sizeof(float))

#define pkv_GetScratchMemd(size) \
  (double*)pkv_GetScratchMem ( (size)*sizeof(double) )
#define pkv_FreeScratchMemd(size) \
  pkv_FreeScratchMem ( (size)*sizeof(double) )
#define pkv_ScratchMemAvaild() \
  (pkv_ScratchMemAvail()/sizeof(double))


/* square angle measure */
double pkv_SqAngle ( double x, double y );

/* conversion of radians to degree string */
void pkv_RadToDegreeStr ( double angle, char *txt );
boolean pkv_DegreeStrToRad ( char *txt, double *angle );

/* computing binomial coefficients */
int pkv_Binom ( int n, int k );

/* multiple integer counter for multidimensional tasks */
boolean pkv_IncMultiCounter ( int dim, const int *limits, int *cnt );

/* Sorting arrays of records with numerical keys.               */
/* The method is stable, i.e. it does not change the order      */
/* of records with the same numerical key                       */
/* except that it will reorder +0.0 and -0.0 (float or double). */
/* The procedures need some scratch memory:                     */
/* (256*sizeof(data) + 2*(number of records))*sizeof(int) +     */
/* sizeof(record)                                               */

/* data types for sorting:                */
#define ID_UNSIGNED            0
#define ID_SIGNED_INT          1
#define ID_IEEE754_FLOAT       2
#define ID_IEEE754_DOUBLE      3
#define ID_IEEE754_LONG_DOUBLE 4

/* possible sorting results */

#define SORT_OK        1
#define SORT_NO_MEMORY 0
#define SORT_BAD_DATA  2

char pkv_SortKernel ( size_t key_size, char key_type,
                      size_t item_length, size_t key_offset,
                      unsigned int num_data, void *data,
                      unsigned int *permut );
char pkv_SortPermute ( size_t item_length, unsigned int num_data, void *data,
                       unsigned int *permut );
char pkv_SortPermute2 ( unsigned int num_data,
                        size_t item_length1, void *data1,
                        size_t item_length2, void *data2,
                        unsigned int *permut );
char pkv_SortFast ( size_t key_size, char key_type,
                    size_t item_length, size_t key_offset,
                    unsigned int num_data, void *data );

void pkv_QuickSort ( int n, void *usrptr,
                     boolean (*less) ( int i, int j, void *usrptr ),
                     void (*swap) (int i, int j, void *usrptr ) );


/* queue */
typedef struct {
    int  nmax, itemsize, fr, en;
    char qtab[1];
  } pkv_queue;

pkv_queue *pkv_InitQueue ( int nmax, int itemsize );
/* use PKV_FREE to destroy queues */
void pkv_ResetQueue ( pkv_queue *q );
boolean pkv_QueueEmpty ( pkv_queue *q );
boolean pkv_QueueFull ( pkv_queue *q );
boolean pkv_QueueInsert ( pkv_queue *q, void *item );
boolean pkv_QueueGetFirst ( pkv_queue *q, void *item );
boolean pkv_QueueRemoveFirst ( pkv_queue *q, void *item );


/* heap priority queue */

int pkv_UpHeap ( void *a[], int l, boolean (*cmp)(void*,void*) );
int pkv_DownHeap ( void *a[], int l, int f, boolean (*cmp)(void*,void*) );
int pkv_HeapInsert ( void *a[], int *l, void *newelem,
                     boolean (*cmp)(void*,void*) );
void pkv_HeapRemove ( void *a[], int *l, int el,
                      boolean (*cmp)(void*,void*) );
void pkv_HeapOrder ( void *a[], int n, boolean (*cmp)(void*,void*) );
void pkv_HeapSort ( void *a[], int n, boolean (*cmp)(void*,void*) );


/* char array processing */

void pkv_ZeroMatc ( int nrows, int rowlen, int pitch, char *data );
void pkv_Rearrangec ( int nrows, int rowlen, int inpitch, int outpitch,
                      char *data );
void pkv_Selectc ( int nrows, int rowlen, int inpitch, int outpitch,
                   const char *indata, char *outdata );
void pkv_Movec ( int nrows, int rowlen, int pitch, int shift, char *data );
void pkv_ReverseMatc ( int nrows, int rowlen, int pitch, char *data );

/* float and double array processing - by macros */

#define pkv_ZeroMatf(nrows,rowlen,pitch,data) \
  pkv_ZeroMatc(nrows,(rowlen)*sizeof(float), \
    (pitch)*sizeof(float),(char*)data)
#define pkv_Rearrangef(nrows,rowlen,inpitch,outpitch,data) \
  pkv_Rearrangec(nrows,(rowlen)*sizeof(float),(inpitch)*sizeof(float), \
    (outpitch)*sizeof(float),(char*)data)
#define pkv_Selectf(nrows,rowlen,inpitch,outpitch,indata,outdata) \
  pkv_Selectc(nrows,(rowlen)*sizeof(float),(inpitch)*sizeof(float), \
    (outpitch)*sizeof(float),(char*)indata,(char*)outdata)
#define pkv_Movef(nrows,rowlen,pitch,shift,data) \
  pkv_Movec(nrows,(rowlen)*sizeof(float),(pitch)*sizeof(float), \
    (shift)*sizeof(float),(char*)data)
#define pkv_ReverseMatf(nrows,rowlen,pitch,data) \
  pkv_ReverseMatc ( nrows, (rowlen)*sizeof(float), (pitch)*sizeof(float), \
    (char*)data )

#define pkv_ZeroMatd(nrows,rowlen,pitch,data) \
  pkv_ZeroMatc(nrows,(rowlen)*sizeof(double), \
    (pitch)*sizeof(double),(char*)data)
#define pkv_Rearranged(nrows,rowlen,inpitch,outpitch,data) \
  pkv_Rearrangec(nrows,(rowlen)*sizeof(double),(inpitch)*sizeof(double), \
    (outpitch)*sizeof(double),(char*)data)
#define pkv_Selectd(nrows,rowlen,inpitch,outpitch,indata,outdata) \
  pkv_Selectc(nrows,(rowlen)*sizeof(double),(inpitch)*sizeof(double), \
    (outpitch)*sizeof(double),(char*)indata,(char*)outdata)
#define pkv_Moved(nrows,rowlen,pitch,shift,data) \
  pkv_Movec(nrows,(rowlen)*sizeof(double),(pitch)*sizeof(double), \
    (shift)*sizeof(double),(char*)data)
#define pkv_ReverseMatd(nrows,rowlen,pitch,data) \
  pkv_ReverseMatc ( nrows, (rowlen)*sizeof(double), (pitch)*sizeof(double), \
    (char*)data )


/* matrix transposition */

void pkv_TransposeMatrixc ( int nrows, int ncols, int elemsize,
                            int inpitch, const char *indata, 
                            int outpitch, char *outdata );

#define pkv_TransposeMatrixf(nrows,ncols,inpitch,indata,outpitch,outdata) \
  pkv_TransposeMatrixc ( nrows, ncols, sizeof(float), \
    (inpitch)*sizeof(float), (char*)indata, \
    (outpitch)*sizeof(float), (char*)outdata )
#define pkv_TransposeMatrixd(nrows,ncols,inpitch,indata,outpitch,outdata) \
  pkv_TransposeMatrixc ( nrows, ncols, sizeof(double), \
    (inpitch)*sizeof(double), (char*)indata, \
    (outpitch)*sizeof(double), (char*)outdata )


/* float <-> double array conversion */

void pkv_Selectfd ( int nrows, int rowlen, int inpitch, int outpitch,
                    const float *indata, double *outdata );
void pkv_Selectdf ( int nrows, int rowlen, int inpitch, int outpitch,
                    const double *indata, float *outdata );


/* computing integer powers of real numbers */

double pkv_rpower ( double x, int e );


/* exchanging variables */

void pkv_Exchange ( void *x, void *y, int size );
void pkv_Sort2f ( float *a, float *b );
void pkv_Sort2d ( double *a, double *b );


/* sign inspection */

char pkv_signf ( float x );
char pkv_signd ( double x );


/* miscellanea */
void pkv_HexByte ( byte b, char *s );


/* malloc and free wrappers for concurrent processes */
/* If the application is supposed to intercept signals, it */
/* must provide a signal procedure, which, if the variable */
/* pkv_critical is true, must set the variable pkv_signal  */
/* and return; then it will be called after completing the */
/* task of memory allocation/deallocation. Only if         */
/* pkv_critical is false, long jump may be done. Similar   */
/* macros may be used to wrap any noninterruptible action. */
extern boolean pkv_critical, pkv_signal;
extern void (*pkv_signal_handler)( void );
extern void (*pkv_register_memblock)( void *ptr, boolean alloc );

#define PKV_MALLOC(ptr,size) \
  { \
    if ( pkv_signal_handler ) { \
      pkv_signal = false; \
      pkv_critical = true; \
      (ptr) = malloc ( size ); \
      if ( pkv_register_memblock ) \
        pkv_register_memblock ( (void*)(ptr), true ); \
      pkv_critical = false; \
      if ( pkv_signal ) \
        pkv_signal_handler (); \
    } \
    else \
      (ptr) = malloc ( size ); \
/*printf ( "malloc: size = %d, ptr = 0x%x\n", size, (size_t)(ptr) ); */ \
  }

#define PKV_FREE(ptr) \
  { \
/*printf ( "free: ptr = 0x%x\n", (size_t)(ptr) ); */ \
    if ( pkv_signal_handler ) { \
      pkv_signal = false; \
      pkv_critical = true; \
      free ( (void*)(ptr) ); \
      if ( pkv_register_memblock ) \
        pkv_register_memblock ( (void*)(ptr), false ); \
      (ptr) = NULL; \
      pkv_critical = false; \
      if ( pkv_signal ) \
        pkv_signal_handler (); \
    } \
    else { \
      free ( (void*)(ptr) ); \
      (ptr) = NULL; \
    } \
  }


/* straight line rasterization*/

typedef struct {
  short x, y;
} xpoint;         /* intentionally identical to XPoint of Xlib.h */

#define PKV_BUFSIZE 256

extern void   (*_pkv_OutputPixels)(const xpoint *buf, int n);
extern xpoint *_pkv_pixbuf;
extern int    _pkv_npix;

#define PKV_FLUSH \
  { _pkv_OutputPixels ( _pkv_pixbuf, _pkv_npix );  _pkv_npix = 0; }
#define PKV_PIXEL(p,px) \
  { (px).x = (short)((p).x/(p).z+0.5);  (px).y = (short)((p).y/(p).z+0.5); }
#define PKV_SETPIXEL(xx,yy) \
  { if ( _pkv_npix == PKV_BUFSIZE ) { PKV_FLUSH; _pkv_npix = 0; } \
    _pkv_pixbuf[_pkv_npix].x = (short)xx;  _pkv_pixbuf[_pkv_npix].y = (short)yy; \
    _pkv_npix++; }

void _pkv_InitPixelBuffer ( void );
void _pkv_DestroyPixelBuffer ( void );
void _pkv_DrawLine ( int x1, int y1, int x2, int y2 );
void pkv_DrawLine ( int x1, int y1, int x2, int y2,
                    void (*output)(const xpoint *buf, int n) );

/* time measurement */
void pkv_Tic ( clock_t *_tic );
int pkv_Toc ( clock_t *_toc );
float pkv_Seconds ( clock_t ticks );

/* some stuff useful for debugging */

void WriteArrayf ( const char *name, int lgt, const float *tab );
void WriteArrayd ( const char *name, int lgt, const double *tab );

void *DMalloc ( size_t size );
void DFree ( void *ptr );

#ifdef __cplusplus
}
#endif

#endif /*PKVARIA_H*/

