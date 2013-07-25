
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/bookg1holef.o $(SRCDIR)/bookg1holed.o

LIBS = -lmultibs -lcamera -lgeom -lpknum -lpkvaria -lX11 -lm

default: $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/bookg1holef.o: $(SRCDIR)/bookg1holef.c $(INCDIR)/bookg1holef.h

$(SRCDIR)/bookg1holed.o: $(SRCDIR)/bookg1holed.c $(INCDIR)/bookg1holed.h

clean:

mrproper:
	rm -f *.o

