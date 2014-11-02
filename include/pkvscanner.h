
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef PKVSCANNER_H
#define PKVSCANNER_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* symbol identifiers - somewhat random, to be rethought */
#define PKV_SYMB_ERROR      -2
#define PKV_SYMB_EOF        -1
#define PKV_SYMB_NONE        0
#define PKV_SYMB_PLUS        1
#define PKV_SYMB_MINUS       2
#define PKV_SYMB_STAR        3
#define PKV_SYMB_SLASH       4
#define PKV_SYMB_LPAREN      5
#define PKV_SYMB_RPAREN      6
#define PKV_SYMB_LBRACKET    7
#define PKV_SYMB_RBRACKET    8
#define PKV_SYMB_LBRACE      9
#define PKV_SYMB_RBRACE     10
#define PKV_SYMB_LANGLE     11
#define PKV_SYMB_RANGLE     12
#define PKV_SYMB_EQUAL      13
#define PKV_SYMB_INTEGER    14
#define PKV_SYMB_FLOAT      15
#define PKV_SYMB_CHAR       16
#define PKV_SYMB_STRING     17
#define PKV_SYMB_COMMA      18
#define PKV_SYMB_DOT        19
#define PKV_SYMB_SEMICOLON  20
#define PKV_SYMB_COLON      21
#define PKV_SYMB_PERCENT    22
#define PKV_SYMB_HASH       23
#define PKV_SYMB_EXCLAM     24
#define PKV_SYMB_TILDE      25
#define PKV_SYMB_AT         26
#define PKV_SYMB_DOLLAR     27
#define PKV_SYMB_ET         28
#define PKV_SYMB_BACKSLASH  29
#define PKV_SYMB_QUESTION   30
#define PKV_SYMB_VERTBAR    31
#define PKV_SYMB_OTHER      32
#define PKV_SYMB_IDENT      33
#define PKV_SYMB_FIRSTOTHER 34

#define PKV_SCANNER_BUFLENGTH 1024

typedef struct pkv_scanner {
    int     inbuflength, namebuflength;
    char    *inbuffer;
    char    *nextname;
    boolean alloc_inb, alloc_nb;
    int     bcomment, ecomment;
    int     inbufpos, inbufcount, namebufcount;
    int     linenum, colnum;
    int     nextchar;
    int     nextsymbol;
    int     nextinteger;
    double  nextfloat;
    int     (*InputData)( void *usrdata, int buflength, char *buffer );
    void    *userdata;
  } pkv_scanner;


boolean pkv_InitScanner ( pkv_scanner *sc,
                          int inbuflength, char *inbuffer,
                          int maxnamelength, char *namebuffer,
                          int bcomment, int ecomment,
                          int (*InputData)(void *userdata,int buflength,char *buffer),
                          void *userdata );
void pkv_GetNextChar ( pkv_scanner *sc );
int pkv_GetNextSymbol ( pkv_scanner *sc );
void pkv_ShutDownScanner ( pkv_scanner *sc );

int pkv_BinSearchNameTable ( int nnames, int firstname, const char **names,
                             const char *name );

#ifdef __cplusplus
}
#endif

#endif

