
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2009
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
G1DIR  = ../../src/projects
BINDIR = ../bin

include ../../options.mak

LDFLAGS =

LIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria -lX11 -lpthread -lm

INCLUDES = $(SRCDIR)/datagen.h $(SRCDIR)/drawbez.h $(SRCDIR)/drawg1hole.h

OBJECTS = $(SRCDIR)/datagen.o $(SRCDIR)/drawbez.o $(SRCDIR)/drawg1hole.o \
  $(G1DIR)/bookg1holef.o

default: $(BINDIR)/polep

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(BINDIR)/polep: $(SRCDIR)/edg1hole.c $(SRCDIR)/edg1hole.h \
  $(INCLUDES) $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/edg1hole.c $(OBJECTS) \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/polep
	rm -f $(SRCDIR)/polep
	ln -s $(BINDIR)/polep $(SRCDIR)/polep

$(SRCDIR)/datagen.o: $(SRCDIR)/datagen.c $(INCLUDES)

$(SRCDIR)/drawbez.o: $(SRCDIR)/drawbez.c $(INCLUDES)

$(SRCDIR)/drawg1hole.o: $(SRCDIR)/drawg1hole.c $(INCLUDES)

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/polep $(SRCDIR)/polep

