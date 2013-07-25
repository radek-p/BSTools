
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */


typedef enum {
  _drawline, _dcircle, _fcircle, _mcircle, _drawrect, _fillrect,
  _hatchrect, _moveto, _lineto, _shcone, _newpath
} ps_proc;

extern FILE *_ps_f;
extern unsigned short _ps_dpi;   /* dots per inch */
extern long _ps_written;
extern float _ps_cgray;   /* current grey level */   
extern float _ps_cred, _ps_cgreen, _ps_cblue;  /* current colour components */
extern float _ps_cwidth;   /* current line width */   

extern short _ps_bmp_y;   /* lines to the end of bitmap */
extern unsigned short _ps_bmp_w;  /* size of the bitmap */
extern byte _ps_bmp_p;   /* bitmap pixel size */


void _ps_ExtendRect ( float x, float y );
void _ps_OutProc ( ps_proc p );

extern boolean _psl_trmk, _psl_btrmk, _psl_htrmk, _psl_bhtrmk;

void _psl_InitPSLib ( void );
