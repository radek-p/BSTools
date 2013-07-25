
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is private for g2hole library functions; */
/* it is NOT intended to be #included in application source files */

#define G2H_ERROR_NO_SCRATCH_MEMORY       1
#define G2H_ERROR_CANNOT_MALLOC           2
#define G2H_ERROR_INVALID_OPTION          3
#define G2H_ERROR_INVALID_PARTITION       4
#define G2H_ERROR_INVALID_JUNC_FUNC       5
#define G2H_ERROR_NONPOSITIVE_MATRIX      6
#define G2H_ERROR_NONPOSITIVE_EXT_MATRIX  7
#define G2H_ERROR_UNDEFINED_CONSTR        8
#define G2H_ERROR_INCONSISTENT_CONSTR     9
#define G2H_ERROR_NL_CANNOT_PROJECT      10
#define G2H_ERROR_NL_JACOBIAN            11
#define G2H_ERROR_NL_MINIMIZATION        12
#define G2H_ERROR_SPLINE_BASIS_NOT_READY 13

char *_g2h_GetErrorString ( int errorcode );

