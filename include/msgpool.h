
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef MSGPOOL_H
#define MSGPOOL_H

/* BSTools library identifiers (for error handling) */
#define LIB_PKVARIA     0
#define LIB_PKNUM       1
#define LIB_GEOM        2
#define LIB_CAMERA      3
#define LIB_PSOUT       4
#define LIB_MULTIBS     5
#define LIB_RAYBEZ      6
#define LIB_EGHOLE      7
#define LIB_G1BLENDING  8
#define LIB_G2BLENDING  9
#define LIB_BSFILE     10
#define LIB_BSMESH     11
#define LIB_XGEDIT     12

/* error codes and descriptions */
#define ERRCODE_0 0
#define ERRMSG_0  "Scratch memory pool not initialized"
#define ERRCODE_1 1
#define ERRMSG_1  "Cannot allocate that large scratch memory pool"
#define ERRCODE_2 2
#define ERRMSG_2  "Not enough scratch memory"
#define ERRCODE_3 3
#define ERRMSG_3  "Invalid parameter"
#define ERRCODE_4 4
#define ERRMSG_4  "Invalid attribute dim_case"
#define ERRCODE_5 5
#define ERRMSG_5  "Invalid data"
#define ERRCODE_6 6
#define ERRMSG_6  "Internal error"
#define ERRCODE_7 7
#define ERRMSG_7  "Point at infinity"
#define ERRCODE_8 8
#define ERRMSG_8  "Degree too high"
#define ERRCODE_9 9
#define ERRMSG_9  "Cannot malloc"
/*
#define ERRCODE_10 10
#define ERRMSG_10 "Invalid data"
*/
#define ERRCODE_11 11
#define ERRMSG_11 "Cannot evaluate basis functions"
#define ERRCODE_12 12
#define ERRMSG_12 "Iterations limit exceeded"
#define ERRCODE_13 13
#define ERRMSG_13 "Bracketing failed"
#define ERRCODE_14 14
#define ERRMSG_14 "Failed to improve a point"
#define ERRCODE_15 15
#define ERRMSG_15 "Cannot decompose constraint equations matrix"
#define ERRCODE_16 16
#define ERRMSG_16 "Cannot transform the Hessian"
#define ERRCODE_17 17
#define ERRMSG_17 "Cannot find rows"
#define ERRCODE_18 18
#define ERRMSG_18 "Could not find neighbours"
#define ERRCODE_19 19
#define ERRMSG_19 "Could not order points"
#define ERRCODE_20 20
#define ERRMSG_20 "Failed to set up a small block"
#define ERRCODE_21 21
#define ERRMSG_21 "Failed to set up a big block"
#define ERRCODE_22 22
#define ERRMSG_22 "Something wrong with the mesh"
#define ERRCODE_23 23
#define ERRMSG_23 "Failed to set up the elements"
#define ERRCODE_24 24
#define ERRMSG_24 "Failed to set up normal vectors"
#define ERRCODE_25 25
#define ERRMSG_25 "Failed to extract block Hessian rows"
#define ERRCODE_26 26
#define ERRMSG_26 "Block Hessian not positive-definite"
#define ERRCODE_27 27
#define ERRMSG_27 "Could not compute the Hessian matrix"
#define ERRCODE_28 28
#define ERRMSG_28 "Could not compute the gradient"
#define ERRCODE_29 29
#define ERRMSG_29 "CG method failed"
#define ERRCODE_30 30
#define ERRMSG_30 "No progress"
#define ERRCODE_31 31
#define ERRMSG_31 "Cannot begin bracketing"
#define ERRCODE_32 32
#define ERRMSG_32 "Failed to set up the coarse mesh preconditioner"

#endif
