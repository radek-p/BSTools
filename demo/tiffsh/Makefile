
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2011, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

#LDFLAGS = -static
LDFLAGS =

LIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -ltiff -lX11 -lXext -lpthread -lm

default: $(BINDIR)/tiffsh

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(BINDIR)/tiffsh: $(SRCDIR)/tiffsh.c $(SRCDIR)/tiffsh.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(SRCDIR)/tiffsh.c \
  -I$(INCDIR) $(LIBDIR) $(LIBS) -o $(BINDIR)/tiffsh
	rm -f $(SRCDIR)/tiffsh
	ln -s $(BINDIR)/tiffsh $(SRCDIR)/tiffsh

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/tiffsh $(SRCDIR)/tiffsh

