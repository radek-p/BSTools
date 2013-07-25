
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2010, 2013
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/g1blsize.o $(SRCDIR)/g1blhprofile.o

OBJECTSF = $(SRCDIR)/g1blbiharmf.o $(SRCDIR)/g1blcbiharmf.o

OBJECTSD = $(SRCDIR)/g1blbiharmd.o $(SRCDIR)/g1blcbiharmd.o \
  $(SRCDIR)/g1blutabfd.o $(SRCDIR)/g1blufgh01d.o $(SRCDIR)/g1blufgh02d.o \
  $(SRCDIR)/g1blopt01d.o $(SRCDIR)/g1blopt02d.o $(SRCDIR)/g1blopt03d.o

TARGET = $(LIBDIR)/libg1blending.a

default: $(TARGET)

$(TARGET): $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/g1blbiharmf.o: $(SRCDIR)/g1blbiharmf.c $(INCDIR)/g1blendingf.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blcbiharmf.o: $(SRCDIR)/g1blcbiharmf.c $(INCDIR)/g1blendingf.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blbiharmd.o: $(SRCDIR)/g1blbiharmd.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blcbiharmd.o: $(SRCDIR)/g1blcbiharmd.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blsize.o: $(SRCDIR)/g1blsize.c \
  $(INCDIR)/g1blendingf.h $(INCDIR)/g1blendingd.h $(SRCDIR)/g1blprivated.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blhprofile.o: $(SRCDIR)/g1blhprofile.c \
  $(INCDIR)/g1blendingf.h $(INCDIR)/g1blendingd.h $(SRCDIR)/g1blprivated.h \
  $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blutabfd.o: $(SRCDIR)/g1blutabfd.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blufgh01d.o: $(SRCDIR)/g1blufgh01d.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blufgh02d.o: $(SRCDIR)/g1blufgh02d.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blopt01d.o: $(SRCDIR)/g1blopt01d.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blopt02d.o: $(SRCDIR)/g1blopt02d.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h

$(SRCDIR)/g1blopt03d.o: $(SRCDIR)/g1blopt03d.c $(INCDIR)/g1blendingd.h \
  $(SRCDIR)/g1blprivated.h $(SRCDIR)/msgpool.h


clean:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

mrproper:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	rm -f $(TARGET)

