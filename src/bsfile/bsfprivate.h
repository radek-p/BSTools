
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef BSFPRIVATE_H
#define BSFPRIVATE_H

/* ////////////////////////////////////////////////////////////////////////// */
/* tokens */
#define BSF_SYMB_EOF               0
#define BSF_SYMB_ERROR             1
#define BSF_SYMB_INTEGER           2
#define BSF_SYMB_FLOAT             3
#define BSF_SYMB_LBRACE            4
#define BSF_SYMB_RBRACE            5
#define BSF_SYMB_PLUS              6
#define BSF_SYMB_MINUS             7
#define BSF_SYMB_STRING            8
#define BSF_SYMB_COMMA             9

/* the subsequent keyword identifiers below must */
/* match the strings in the table in bsfile00r.c */
#define BSF_FIRST_KEYWORD          10  /* number of the first keyword */
#define BSF_SYMB_BCURVE            10  /* consecutive numbers of keywords */
#define BSF_SYMB_BPATCH            11  /* sorted alphabetically, with uppercase */
#define BSF_SYMB_BSCURVE           12  /* preceding all lowercase letters (ASCII) */
#define BSF_SYMB_BSHOLE            13
#define BSF_SYMB_BSMESH            14
#define BSF_SYMB_BSPATCH           15
#define BSF_SYMB_EULERANGLES       16
#define BSF_SYMB_CAMERA            17
#define BSF_SYMB_CLOSED            18
#define BSF_SYMB_COLOR             19  /* it can read both, AE and BE, */
#define BSF_SYMB_COLOUR            20  /* but it will write in BE */
#define BSF_SYMB_CPOINTS           21
#define BSF_SYMB_CPOINTSMK         22
#define BSF_SYMB_DEGREE            23
#define BSF_SYMB_DEPTH             24
#define BSF_SYMB_DIM               25
#define BSF_SYMB_DOMAIN            26
#define BSF_SYMB_FACETMK           27
#define BSF_SYMB_FACETS            28
#define BSF_SYMB_FRAME             29
#define BSF_SYMB_HALFEDGEMK        30
#define BSF_SYMB_HALFEDGES         31
#define BSF_SYMB_IDENT             32
#define BSF_SYMB_KNOTS             33
#define BSF_SYMB_KNOTS_U           34
#define BSF_SYMB_KNOTS_V           35
#define BSF_SYMB_NAME              36
#define BSF_SYMB_PARALLEL          37
#define BSF_SYMB_PERSPECTIVE       38
#define BSF_SYMB_POINTS            39
#define BSF_SYMB_POLYLINE          40
#define BSF_SYMB_POSITION          41
#define BSF_SYMB_RATIONAL          42
#define BSF_SYMB_SIDES             43
#define BSF_SYMB_SPHERICAL_PRODUCT 44
#define BSF_SYMB_TRIMMED           45
#define BSF_SYMB_UNIFORM           46
#define BSF_SYMB_VERTICES          47
#define BSF_LAST_KEYWORD           48  /* number of the last keyword */

/* namely, this table */
#define BSF_NKEYWORDS (BSF_LAST_KEYWORD-BSF_FIRST_KEYWORD+1)

extern const char *bsf_keyword[BSF_NKEYWORDS];

/* ////////////////////////////////////////////////////////////////////////// */
/* static variables */
extern FILE    *bsf_input, *bsf_output;

extern int     bsf_nextsymbol;
extern int     bsf_nextint;
extern double  bsf_nextfloat;
extern char    *bsf_namebuffer;

/* variables for output formatting */
extern boolean bsf_newline;
extern int     bsf_current_length;
extern int     bsf_current_indentation;

/* ////////////////////////////////////////////////////////////////////////// */
/* auxiliary, low-level reading procedures */
void bsf_GetNextSymbol ( void );

boolean bsf_ReadIntNumber ( int *number );
boolean bsf_ReadDoubleNumber ( double *number );
boolean bsf_ReadPointd ( int maxcpdimen, double *point, int *cpdimen );
boolean bsf_ReadIdent ( int *ident );
int bsf_ReadPointsd ( int maxcpdimen, int maxnpoints,
                      double *points, int *cpdimen );
int bsf_ReadPointsMK ( int maxnpoints, unsigned int *mk );
int bsf_ReadHEdgeMK ( int maxnhe, unsigned int *hemk );
int bsf_ReadFacetMK ( int maxnfac, unsigned int *fmk );

boolean bsf_ReadSpaceDim ( int maxdim, int *spdimen );
boolean bsf_ReadCurveDegree ( int maxdeg, int *degree );
boolean bsf_ReadPatchDegree ( int maxdeg, int *udeg, int *vdeg );

boolean bsf_ReadKnotSequenced ( int maxlastknot, int *lastknot, double *knots,
                                boolean *closed );

boolean _bsf_ReadCPMark ( bsf_UserReaders *readers, int maxnpoints );
boolean _bsf_ReadHEdgeMark ( bsf_UserReaders *readers, int maxnhe );
boolean _bsf_ReadFacetMark ( bsf_UserReaders *readers, int maxnfac );
boolean _bsf_ReadColour ( bsf_UserReaders *readers );

/* ////////////////////////////////////////////////////////////////////////// */
/* output file formatting procedures and macros */

#define BSF_OUTPUT_LINE_LENGTH 64

void bsf_EndOutputLine ( void );
void bsf_DivideOutputLine ( void );
void bsf_WriteCurrentIndentation ( void );

#define BSFeol bsf_EndOutputLine ();
#define BSFdol bsf_DivideOutputLine ();
#define BSFwci bsf_WriteCurrentIndentation ();

#endif

