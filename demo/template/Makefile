
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

#LDFLAGS = -static
LDFLAGS =

OBJECTS = edtempwidgets.o edtemplate01.o edtemplate02.o edtemplate03.o \
  edtemplate04.o edtemplate05.o edtemplate06.o spltemplate.o

LIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria -lX11 -lpthread -lm

EDINCL = $(SRCDIR)/spltemplate.h $(SRCDIR)/edtemplate.h $(SRCDIR)/edtempwidgets.h

default: $(BINDIR)/template

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/edtempwidgets.o: $(SRCDIR)/edtempwidgets.c $(EDINCL)

$(SRCDIR)/edtemplate01.o: $(SRCDIR)/edtemplate01.c $(EDINCL)

$(SRCDIR)/edtemplate02.o: $(SRCDIR)/edtemplate02.c $(EDINCL)

$(SRCDIR)/edtemplate03.o: $(SRCDIR)/edtemplate03.c $(EDINCL)

$(SRCDIR)/edtemplate04.o: $(SRCDIR)/edtemplate04.c $(EDINCL)

$(SRCDIR)/edtemplate05.o: $(SRCDIR)/edtemplate05.c $(EDINCL)

$(SRCDIR)/edtemplate06.o: $(SRCDIR)/edtemplate06.c $(EDINCL)

$(SRCDIR)/spltemplate.o: $(SRCDIR)/spltemplate.c $(EDINCL)

$(BINDIR)/template: $(SRCDIR)/edtemplate.c $(EDINCL) $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/edtemplate.c \
  $(OBJECTS) $(LIBDIR) $(LIBS) -o $(BINDIR)/template

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/template

