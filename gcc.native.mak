
# //////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# //////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2006, 2013
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# //////////////////////////////////////////////////////////////////////

# the definitions below cause using the gcc compiler with the
# optimization options on, and debugging off

CC = gcc

CFLAGS = -ansi -pedantic -pthread -Wall \
  -O2 -funroll-loops -finline-functions -funswitch-loops -Wall -march=native

AR = ar
ARFLAGS = rv

