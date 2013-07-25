
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void g1h_SetOptionProcd ( GHoleDomaind *domain,  
    int (*OptionProc)( GHoleDomaind *domain, int query, int qn,
                       int *ndata, int **idata, double **fdata ) )
{
  G1HolePrivateRecd *privateG1;

  privateG1 = domain->privateG1;
  if ( !OptionProc )
    privateG1->GetOption = gh_GetDefaultOptiond;
  else
    privateG1->GetOption = OptionProc;
} /*g1h_SetOptionProcd*/

