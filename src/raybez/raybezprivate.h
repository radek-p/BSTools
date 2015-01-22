
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

typedef void (*divide_vertex_proc)(void *tree,void *vertex);

void raybez_DivideTreeVertex ( void *tree, void *vertex, unsigned char *tag,
                               divide_vertex_proc DivideVertex );

