
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2012, 2014
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
  -lX11 -lXext -lpthread -lm

GLIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lGL -lX11 -lXext -lpthread -lm

default: $(BINDIR)/ifs

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o

$(BINDIR)/ifs: $(SRCDIR)/ifs.c $(SRCDIR)/ifs.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ifs.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/ifs
	rm -f $(SRCDIR)/ifs
	ln -s $(BINDIR)/ifs $(SRCDIR)/ifs

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/ifs

