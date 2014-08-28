
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2009, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/bsfile00r.o \
  $(SRCDIR)/bsfile01r.o $(SRCDIR)/bsfile02r.o $(SRCDIR)/bsfile03r.o \
  $(SRCDIR)/bsfile04r.o $(SRCDIR)/bsfile05r.o $(SRCDIR)/bsfile06r.o \
  $(SRCDIR)/bsfile07r.o $(SRCDIR)/bsfile08r.o $(SRCDIR)/bsfile09r.o \
  $(SRCDIR)/bsfile10r.o $(SRCDIR)/bsfile11r.o $(SRCDIR)/bsfile12r.o \
  $(SRCDIR)/bsfile13r.o $(SRCDIR)/bsfile14r.o $(SRCDIR)/bsfile15r.o \
  $(SRCDIR)/bsfile16r.o $(SRCDIR)/bsfile17r.o \
  $(SRCDIR)/bsfile00w.o \
  $(SRCDIR)/bsfile01w.o $(SRCDIR)/bsfile02w.o $(SRCDIR)/bsfile03w.o \
  $(SRCDIR)/bsfile04w.o $(SRCDIR)/bsfile05w.o $(SRCDIR)/bsfile06w.o \
  $(SRCDIR)/bsfile07w.o $(SRCDIR)/bsfile08w.o $(SRCDIR)/bsfile09w.o \
  $(SRCDIR)/bsfile10w.o $(SRCDIR)/bsfile11w.o $(SRCDIR)/bsfile12w.o \
  $(SRCDIR)/bsfile13w.o $(SRCDIR)/bsfile14w.o $(SRCDIR)/bsfile15w.o \
  $(SRCDIR)/bsfile16w.o $(SRCDIR)/bsfile17w.o

TARGET = $(LIBDIR)/libbsfile.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/bsfile00r.o: $(SRCDIR)/bsfile00r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile01r.o: $(SRCDIR)/bsfile01r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile02r.o: $(SRCDIR)/bsfile02r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile03r.o: $(SRCDIR)/bsfile03r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile04r.o: $(SRCDIR)/bsfile04r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile05r.o: $(SRCDIR)/bsfile05r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile06r.o: $(SRCDIR)/bsfile06r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile07r.o: $(SRCDIR)/bsfile07r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile08r.o: $(SRCDIR)/bsfile08r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile09r.o: $(SRCDIR)/bsfile09r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile10r.o: $(SRCDIR)/bsfile10r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile11r.o: $(SRCDIR)/bsfile11r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile12r.o: $(SRCDIR)/bsfile12r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile13r.o: $(SRCDIR)/bsfile13r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile14r.o: $(SRCDIR)/bsfile14r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile15r.o: $(SRCDIR)/bsfile15r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile16r.o: $(SRCDIR)/bsfile16r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile17r.o: $(SRCDIR)/bsfile17r.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile00w.o: $(SRCDIR)/bsfile00w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile01w.o: $(SRCDIR)/bsfile01w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile02w.o: $(SRCDIR)/bsfile02w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile03w.o: $(SRCDIR)/bsfile03w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile04w.o: $(SRCDIR)/bsfile04w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile05w.o: $(SRCDIR)/bsfile05w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile06w.o: $(SRCDIR)/bsfile06w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile07w.o: $(SRCDIR)/bsfile07w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile08w.o: $(SRCDIR)/bsfile08w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile09w.o: $(SRCDIR)/bsfile09w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile10w.o: $(SRCDIR)/bsfile10w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile11w.o: $(SRCDIR)/bsfile11w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile12w.o: $(SRCDIR)/bsfile12w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile13w.o: $(SRCDIR)/bsfile13w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile14w.o: $(SRCDIR)/bsfile14w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile16w.o: $(SRCDIR)/bsfile16w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

$(SRCDIR)/bsfile17w.o: $(SRCDIR)/bsfile17w.c $(INCDIR)/bsfile.h $(SRCDIR)/bsfprivate.h

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

