
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef PKVSCANNER_H
#include "pkvscanner.h"
#endif

#ifndef BSFPRIVATE_H
#define BSFPRIVATE_H

/* ////////////////////////////////////////////////////////////////////////// */
/* tokens */
#define BSF_SYMB_EOF               PKV_SYMB_EOF
#define BSF_SYMB_ERROR             PKV_SYMB_ERROR
#define BSF_SYMB_INTEGER           PKV_SYMB_INTEGER
#define BSF_SYMB_FLOAT             PKV_SYMB_FLOAT
#define BSF_SYMB_LBRACE            PKV_SYMB_LBRACE
#define BSF_SYMB_RBRACE            PKV_SYMB_RBRACE
#define BSF_SYMB_PLUS              PKV_SYMB_PLUS
#define BSF_SYMB_MINUS             PKV_SYMB_MINUS
#define BSF_SYMB_STRING            PKV_SYMB_STRING
#define BSF_SYMB_COMMA             PKV_SYMB_COMMA

/* the subsequent keyword identifiers below must */
/* match the strings in the table in bsfile00r.c */
#define BSF_FIRST_KEYWORD          PKV_SYMB_FIRSTOTHER  /* number of the first keyword */
#define BSF_SYMB_BCURVE            PKV_SYMB_FIRSTOTHER  /* consecutive numbers of keywords */
#define BSF_SYMB_BPATCH            (PKV_SYMB_FIRSTOTHER+1)   /* sorted alphabetically, with uppercase */
#define BSF_SYMB_BSCURVE           (PKV_SYMB_FIRSTOTHER+2)   /* preceding all lowercase letters (ASCII) */
#define BSF_SYMB_BSHOLE            (PKV_SYMB_FIRSTOTHER+3)
#define BSF_SYMB_BSMESH            (PKV_SYMB_FIRSTOTHER+4)
#define BSF_SYMB_BSPATCH           (PKV_SYMB_FIRSTOTHER+5)
#define BSF_SYMB_EULERANGLES       (PKV_SYMB_FIRSTOTHER+6)
#define BSF_SYMB_CAMERA            (PKV_SYMB_FIRSTOTHER+7)
#define BSF_SYMB_CLOSED            (PKV_SYMB_FIRSTOTHER+8)
#define BSF_SYMB_COLOR             (PKV_SYMB_FIRSTOTHER+9)   /* it can read both, AE and BE, */
#define BSF_SYMB_COLOUR            (PKV_SYMB_FIRSTOTHER+10)  /* but it will write in BE */
#define BSF_SYMB_CPOINTS           (PKV_SYMB_FIRSTOTHER+11)
#define BSF_SYMB_CPOINTSMK         (PKV_SYMB_FIRSTOTHER+12)
#define BSF_SYMB_DEGREE            (PKV_SYMB_FIRSTOTHER+13)
#define BSF_SYMB_DEPTH             (PKV_SYMB_FIRSTOTHER+14)
#define BSF_SYMB_DIM               (PKV_SYMB_FIRSTOTHER+15)
#define BSF_SYMB_DOMAIN            (PKV_SYMB_FIRSTOTHER+16)
#define BSF_SYMB_FACETMK           (PKV_SYMB_FIRSTOTHER+17)
#define BSF_SYMB_FACETS            (PKV_SYMB_FIRSTOTHER+18)
#define BSF_SYMB_FRAME             (PKV_SYMB_FIRSTOTHER+19)
#define BSF_SYMB_HALFEDGEMK        (PKV_SYMB_FIRSTOTHER+20)
#define BSF_SYMB_HALFEDGES         (PKV_SYMB_FIRSTOTHER+21)
#define BSF_SYMB_IDENT             (PKV_SYMB_FIRSTOTHER+22)
#define BSF_SYMB_KNOTS             (PKV_SYMB_FIRSTOTHER+23)
#define BSF_SYMB_KNOTS_U           (PKV_SYMB_FIRSTOTHER+24)
#define BSF_SYMB_KNOTS_V           (PKV_SYMB_FIRSTOTHER+25)
#define BSF_SYMB_NAME              (PKV_SYMB_FIRSTOTHER+26)
#define BSF_SYMB_PARALLEL          (PKV_SYMB_FIRSTOTHER+27)
#define BSF_SYMB_PERSPECTIVE       (PKV_SYMB_FIRSTOTHER+28)
#define BSF_SYMB_POINTS            (PKV_SYMB_FIRSTOTHER+29)
#define BSF_SYMB_POLYLINE          (PKV_SYMB_FIRSTOTHER+30)
#define BSF_SYMB_POSITION          (PKV_SYMB_FIRSTOTHER+31)
#define BSF_SYMB_RATIONAL          (PKV_SYMB_FIRSTOTHER+32)
#define BSF_SYMB_SIDES             (PKV_SYMB_FIRSTOTHER+33)
#define BSF_SYMB_SPHERICAL_PRODUCT (PKV_SYMB_FIRSTOTHER+34)
#define BSF_SYMB_TRIMMED           (PKV_SYMB_FIRSTOTHER+35)
#define BSF_SYMB_UNIFORM           (PKV_SYMB_FIRSTOTHER+36)
#define BSF_SYMB_VERTICES          (PKV_SYMB_FIRSTOTHER+37)
#define BSF_LAST_KEYWORD           (PKV_SYMB_FIRSTOTHER+38)  /* number of the last keyword */

/* namely, this table */
#define BSF_NKEYWORDS (BSF_LAST_KEYWORD-BSF_FIRST_KEYWORD+1)

extern const char *bsf_keyword[BSF_NKEYWORDS];

/* ////////////////////////////////////////////////////////////////////////// */
/* static variables */
extern FILE    *bsf_input, *bsf_output;

extern int     bsf_nextsymbol;
extern int     bsf_nextint;
extern double  bsf_nextfloat;
extern char    *bsf_name;

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

