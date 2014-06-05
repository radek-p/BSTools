
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/mengerc00.o $(SRCDIR)/mengerc01.o \
  $(SRCDIR)/mengerc02.o $(SRCDIR)/mengerc03.o $(SRCDIR)/mengerc04.o \
  $(SRCDIR)/mengerc05.o

TARGET = $(LIBDIR)/libmengerc.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) -c $(CFLAGS) -I$(INCDIR) $< -o $*.o

$(SRCDIR)/mengerc00.o: $(SRCDIR)/mengerc00.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

$(SRCDIR)/mengerc01.o: $(SRCDIR)/mengerc01.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

$(SRCDIR)/mengerc02.o: $(SRCDIR)/mengerc02.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

$(SRCDIR)/mengerc03.o: $(SRCDIR)/mengerc03.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

$(SRCDIR)/mengerc04.o: $(SRCDIR)/mengerc04.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

$(SRCDIR)/mengerc05.o: $(SRCDIR)/mengerc05.c $(INCDIR)/mengerc.h \
  $(SRCDIR)/mengercprivate.h

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

