
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2008, 2009
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

LDFLAGS =

LIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lGLU -lGL -lX11 -lpthread -lm

default: $(BINDIR)/glsimplest $(BINDIR)/gl2Dwin $(BINDIR)/gl3Dwin

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o

$(BINDIR)/glsimplest: $(SRCDIR)/glsimplest.c $(SRCDIR)/glsimplest.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/glsimplest.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/glsimplest

$(BINDIR)/gl2Dwin: $(SRCDIR)/gl2Dwin.c $(SRCDIR)/gl2Dwin.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/gl2Dwin.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/gl2Dwin

$(BINDIR)/gl3Dwin: $(SRCDIR)/gl3Dwin.c $(SRCDIR)/gl3Dwin.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/gl3Dwin.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/gl3Dwin

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/glsimplest $(BINDIR)/gl2Dwin $(BINDIR)/gl3Dwin

