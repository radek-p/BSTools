
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void g2h_SetOptionProcf ( GHoleDomainf *domain,  
    int (*OptionProc)( GHoleDomainf *domain, int query, int qn,
                       int *ndata, int **idata, float **fdata ) )
{
  G2HolePrivateRecf *privateG2;

  privateG2 = domain->privateG2;
  if ( !OptionProc )
    privateG2->GetOption = gh_GetDefaultOptionf;
  else
    privateG2->GetOption = OptionProc;
} /*g2h_SetOptionProcf*/

