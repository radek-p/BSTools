
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void g2h_SetOptionProcd ( GHoleDomaind *domain,  
    int (*OptionProc)( GHoleDomaind *domain, int query, int qn,
                       int *ndata, int **idata, double **fdata ) )
{
  G2HolePrivateRecd *privateG2;

  privateG2 = domain->privateG2;
  if ( !OptionProc )
    privateG2->GetOption = gh_GetDefaultOptiond;
  else
    privateG2->GetOption = OptionProc;
} /*g2h_SetOptionProcd*/

