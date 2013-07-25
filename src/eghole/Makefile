
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2010
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/eghbspind.o $(SRCDIR)/eg1herrstr.o $(SRCDIR)/eg2herrstr.o \
  $(SRCDIR)/eg1hsymmatrix.o $(SRCDIR)/eg2hsymmatrix.o

OBJECTSF = $(SRCDIR)/eghcoref.o $(SRCDIR)/eghgraphf.o $(SRCDIR)/eghddiamf.o \
  $(SRCDIR)/egholeareaf.o \
  $(SRCDIR)/eg1hcoref.o $(SRCDIR)/eg1hextf.o \
  $(SRCDIR)/eg1hconstrf.o $(SRCDIR)/eg1hextconstrf.o $(SRCDIR)/eg1hfunclf.o \
  $(SRCDIR)/eg1hgraphf.o $(SRCDIR)/eg1herrorf.o $(SRCDIR)/eg1hdiagf.o \
  $(SRCDIR)/eg1hnonlin1f.o $(SRCDIR)/eg1hnonlin2f.o $(SRCDIR)/eg1hnonlin3f.o \
  $(SRCDIR)/eg1hnonlin4f.o $(SRCDIR)/eg1hpcurvesf.o $(SRCDIR)/eg1hsplpcurvesf.o \
  $(SRCDIR)/eg1hsplgraphf.o $(SRCDIR)/eg1hsplcore1f.o \
  $(SRCDIR)/eg1hsplcore2f.o $(SRCDIR)/eg1hsplcore3f.o \
  $(SRCDIR)/eg1hsplconstrf.o $(SRCDIR)/eg1hsymmatrixf.o \
  $(SRCDIR)/eg1hq2coref.o $(SRCDIR)/eg1hq2extf.o $(SRCDIR)/eg1hq2splf.o \
  $(SRCDIR)/eg1hq2constrf.o $(SRCDIR)/eg1hq2extconstrf.o \
  $(SRCDIR)/eg1hq2graphf.o $(SRCDIR)/eg1hq2nonlin1f.o $(SRCDIR)/eg1hq2nonlin2f.o \
  $(SRCDIR)/eg1hq2nonlin3f.o $(SRCDIR)/eg1hq2nonlin4f.o \
  $(SRCDIR)/eg2hcoref.o $(SRCDIR)/eg2hextf.o \
  $(SRCDIR)/eg2hconstrf.o $(SRCDIR)/eg2hextconstrf.o $(SRCDIR)/eg2hfunclf.o \
  $(SRCDIR)/eg2hgraphf.o $(SRCDIR)/eg2herrorf.o $(SRCDIR)/eg2hdiagf.o \
  $(SRCDIR)/eg2hnonlin1f.o $(SRCDIR)/eg2hnonlin2f.o $(SRCDIR)/eg2hnonlin3f.o \
  $(SRCDIR)/eg2hnonlin4f.o $(SRCDIR)/eg2hpcurvesf.o $(SRCDIR)/eg2hsplpcurvesf.o \
  $(SRCDIR)/eg2hsplgraphf.o $(SRCDIR)/eg2hsplcore1f.o \
  $(SRCDIR)/eg2hsplcore2f.o $(SRCDIR)/eg2hsplcore3f.o \
  $(SRCDIR)/eg2hsplconstrf.o $(SRCDIR)/eg2hsymmatrixf.o \
  $(SRCDIR)/egheigenmeshf.o $(SRCDIR)/egheigenvalf.o

OBJECTSD = $(SRCDIR)/eghcored.o $(SRCDIR)/eghgraphd.o $(SRCDIR)/eghddiamd.o \
  $(SRCDIR)/egholearead.o \
  $(SRCDIR)/eg1hcored.o $(SRCDIR)/eg1hextd.o \
  $(SRCDIR)/eg1hconstrd.o $(SRCDIR)/eg1hextconstrd.o $(SRCDIR)/eg1hfuncld.o \
  $(SRCDIR)/eg1hgraphd.o $(SRCDIR)/eg1herrord.o $(SRCDIR)/eg1hdiagd.o \
  $(SRCDIR)/eg1hnonlin1d.o $(SRCDIR)/eg1hnonlin2d.o $(SRCDIR)/eg1hnonlin3d.o \
  $(SRCDIR)/eg1hnonlin4d.o $(SRCDIR)/eg1hpcurvesd.o $(SRCDIR)/eg1hsplpcurvesd.o \
  $(SRCDIR)/eg1hsplgraphd.o $(SRCDIR)/eg1hsplcore1d.o \
  $(SRCDIR)/eg1hsplcore2d.o $(SRCDIR)/eg1hsplcore3d.o \
  $(SRCDIR)/eg1hsplconstrd.o $(SRCDIR)/eg1hsymmatrixd.o \
  $(SRCDIR)/eg1hq2cored.o $(SRCDIR)/eg1hq2extd.o $(SRCDIR)/eg1hq2spld.o \
  $(SRCDIR)/eg1hq2constrd.o $(SRCDIR)/eg1hq2extconstrd.o \
  $(SRCDIR)/eg1hq2graphd.o $(SRCDIR)/eg1hq2nonlin1d.o $(SRCDIR)/eg1hq2nonlin2d.o \
  $(SRCDIR)/eg1hq2nonlin3d.o $(SRCDIR)/eg1hq2nonlin4d.o \
  $(SRCDIR)/eg2hcored.o $(SRCDIR)/eg2hextd.o \
  $(SRCDIR)/eg2hconstrd.o $(SRCDIR)/eg2hextconstrd.o $(SRCDIR)/eg2hfuncld.o \
  $(SRCDIR)/eg2hgraphd.o $(SRCDIR)/eg2herrord.o $(SRCDIR)/eg2hdiagd.o \
  $(SRCDIR)/eg2hnonlin1d.o $(SRCDIR)/eg2hnonlin2d.o $(SRCDIR)/eg2hnonlin3d.o \
  $(SRCDIR)/eg2hnonlin4d.o $(SRCDIR)/eg2hpcurvesd.o $(SRCDIR)/eg2hsplpcurvesd.o \
  $(SRCDIR)/eg2hsplgraphd.o $(SRCDIR)/eg2hsplcore1d.o \
  $(SRCDIR)/eg2hsplcore2d.o $(SRCDIR)/eg2hsplcore3d.o \
  $(SRCDIR)/eg2hsplconstrd.o $(SRCDIR)/eg2hsymmatrixd.o \
  $(SRCDIR)/egheigenmeshd.o $(SRCDIR)/egheigenvald.o


TARGET = $(LIBDIR)/libeghole.a

default: $(TARGET)

$(TARGET): $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	$(AR) $(ARFLAGS) $(LIBDIR)/libeghole.a $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/eghbspind.o: $(SRCDIR)/eghbspind.c \
  $(INCDIR)/eg1holef.h $(INCDIR)/eg1holed.h \
  $(INCDIR)/eg2holef.h $(INCDIR)/eg2holed.h

$(SRCDIR)/eg1herrstr.o: $(SRCDIR)/eg1herrstr.c \
  $(INCDIR)/eg1holef.h $(INCDIR)/eg1holed.h

$(SRCDIR)/eg2herrstr.o: $(SRCDIR)/eg2herrstr.c \
  $(INCDIR)/eg2holef.h $(INCDIR)/eg2holed.h

$(SRCDIR)/eg1hsymmatrix.o: $(SRCDIR)/eg1hsymmatrix.c \
  $(INCDIR)/eg1holef.h $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1hprivated.h

$(SRCDIR)/eg2hsymmatrix.o: $(SRCDIR)/eg2hsymmatrix.c \
  $(INCDIR)/eg2holef.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2hprivated.h


$(SRCDIR)/eghcoref.o: $(SRCDIR)/eghcoref.c $(INCDIR)/egholef.h \
  $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eghcored.o: $(SRCDIR)/eghcored.c $(INCDIR)/egholed.h \
  $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eghgraphf.o: $(SRCDIR)/eghgraphf.c $(INCDIR)/egholef.h \
  $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eghgraphd.o: $(SRCDIR)/eghgraphd.c $(INCDIR)/egholed.h \
  $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eghddiamf.o: $(SRCDIR)/eghddiamf.c $(INCDIR)/egholef.h \
  $(INCDIR)/eg1holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eghddiamd.o: $(SRCDIR)/eghddiamd.c $(INCDIR)/egholed.h \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/egholeareaf.o: $(SRCDIR)/egholeareaf.c $(INCDIR)/egholef.h \
  $(SRCDIR)/eghprivatef.h

$(SRCDIR)/egholearead.o: $(SRCDIR)/egholearead.c $(INCDIR)/egholed.h \
  $(SRCDIR)/eghprivated.h

$(SRCDIR)/eg1hcoref.o: $(SRCDIR)/eg1hcoref.c $(SRCDIR)/eg1hcore1f.c \
  $(SRCDIR)/eg1hcore2f.c $(SRCDIR)/eg1hcore3f.c $(SRCDIR)/eg1hcore4f.c \
  $(SRCDIR)/eg1h9core4f.c $(SRCDIR)/eg1h9core5f.c \
  $(SRCDIR)/eg1hcore6f.c $(SRCDIR)/eg1hcore7f.c \
  $(SRCDIR)/eg1hcore8f.c $(SRCDIR)/eg1hcore9f.c $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hcored.o: $(SRCDIR)/eg1hcored.c $(SRCDIR)/eg1hcore1d.c \
  $(SRCDIR)/eg1hcore2d.c $(SRCDIR)/eg1hcore3d.c $(SRCDIR)/eg1hcore4d.c \
  $(SRCDIR)/eg1h9core4d.c $(SRCDIR)/eg1h9core5d.c \
  $(SRCDIR)/eg1hcore6d.c $(SRCDIR)/eg1hcore7d.c \
  $(SRCDIR)/eg1hcore8d.c $(SRCDIR)/eg1hcore9d.c $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hdiagf.o: $(SRCDIR)/eg1hdiagf.c $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hdiagd.o: $(SRCDIR)/eg1hdiagd.c $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hextf.o: $(SRCDIR)/eg1hextf.c $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hextd.o: $(SRCDIR)/eg1hextd.c $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hconstrf.o: $(SRCDIR)/eg1hconstrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hconstrd.o: $(SRCDIR)/eg1hconstrd.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hextconstrf.o: $(SRCDIR)/eg1hextconstrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hextconstrd.o: $(SRCDIR)/eg1hextconstrd.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hfunclf.o: $(SRCDIR)/eg1hfunclf.c $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hfuncld.o: $(SRCDIR)/eg1hfuncld.c $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hgraphf.o: $(SRCDIR)/eg1hgraphf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hgraphd.o: $(SRCDIR)/eg1hgraphd.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1herrorf.o: $(SRCDIR)/eg1herrorf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1herrord.o: $(SRCDIR)/eg1herrord.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin1f.o: $(SRCDIR)/eg1hnonlin1f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin1d.o: $(SRCDIR)/eg1hnonlin1d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin2f.o: $(SRCDIR)/eg1hnonlin2f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin2d.o: $(SRCDIR)/eg1hnonlin2d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin3f.o: $(SRCDIR)/eg1hnonlin3f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin3d.o: $(SRCDIR)/eg1hnonlin3d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin4f.o: $(SRCDIR)/eg1hnonlin4f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hnonlin4d.o: $(SRCDIR)/eg1hnonlin4d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore1f.o: $(SRCDIR)/eg1hsplcore1f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore2f.o: $(SRCDIR)/eg1hsplcore2f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore3f.o: $(SRCDIR)/eg1hsplcore3f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplconstrf.o: $(SRCDIR)/eg1hsplconstrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore1d.o: $(SRCDIR)/eg1hsplcore1d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore2d.o: $(SRCDIR)/eg1hsplcore2d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplcore3d.o: $(SRCDIR)/eg1hsplcore3d.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplconstrd.o: $(SRCDIR)/eg1hsplconstrd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplgraphf.o: $(SRCDIR)/eg1hsplgraphf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplgraphd.o: $(SRCDIR)/eg1hsplgraphd.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2coref.o: $(SRCDIR)/eg1hq2coref.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2cored.o: $(SRCDIR)/eg1hq2cored.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2extf.o: $(SRCDIR)/eg1hq2extf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2extd.o: $(SRCDIR)/eg1hq2extd.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2splf.o: $(SRCDIR)/eg1hq2splf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2spld.o: $(SRCDIR)/eg1hq2spld.c \
  $(INCDIR)/eg1holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2constrf.o: $(SRCDIR)/eg1hq2constrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2constrd.o: $(SRCDIR)/eg1hq2constrd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2extconstrf.o: $(SRCDIR)/eg1hq2extconstrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2extconstrd.o: $(SRCDIR)/eg1hq2extconstrd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2graphf.o: $(SRCDIR)/eg1hq2graphf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(SRCDIR)/eg1hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2graphd.o: $(SRCDIR)/eg1hq2graphd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(SRCDIR)/eg1hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin1f.o: $(SRCDIR)/eg1hq2nonlin1f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin1d.o: $(SRCDIR)/eg1hq2nonlin1d.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin2f.o: $(SRCDIR)/eg1hq2nonlin2f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin2d.o: $(SRCDIR)/eg1hq2nonlin2d.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin3f.o: $(SRCDIR)/eg1hq2nonlin3f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin3d.o: $(SRCDIR)/eg1hq2nonlin3d.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin4f.o: $(SRCDIR)/eg1hq2nonlin4f.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hq2nonlin4d.o: $(SRCDIR)/eg1hq2nonlin4d.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hpcurvesf.o: $(SRCDIR)/eg1hpcurvesf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hpcurvesd.o: $(SRCDIR)/eg1hpcurvesd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplpcurvesf.o: $(SRCDIR)/eg1hsplpcurvesf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsplpcurvesd.o: $(SRCDIR)/eg1hsplpcurvesd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg1holed.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg1hsymmatrixf.o: $(SRCDIR)/eg1hsymmatrixf.c $(INCDIR)/eg1holef.h \
  $(SRCDIR)/eg1hprivatef.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg2hsymmatrixf.o: $(SRCDIR)/eg2hsymmatrixf.c $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h


$(SRCDIR)/eg2hcoref.o: $(SRCDIR)/eg2hcoref.c $(SRCDIR)/eg2hcore1f.c \
  $(SRCDIR)/eg2hcore2f.c $(SRCDIR)/eg2hcore3f.c $(SRCDIR)/eg2hcore4f.c \
  $(SRCDIR)/eg2hcore5f.c $(SRCDIR)/eg2hcore6f.c $(SRCDIR)/eg2hcore7f.c \
  $(SRCDIR)/eg2hcore8f.c $(SRCDIR)/eg2hcore9f.c $(SRCDIR)/eg2hcore10f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hcored.o: $(SRCDIR)/eg2hcored.c $(SRCDIR)/eg2hcore1d.c \
  $(SRCDIR)/eg2hcore2d.c $(SRCDIR)/eg2hcore3d.c $(SRCDIR)/eg2hcore4d.c \
  $(SRCDIR)/eg2hcore5d.c $(SRCDIR)/eg2hcore6d.c $(SRCDIR)/eg2hcore7d.c \
  $(SRCDIR)/eg2hcore8d.c $(SRCDIR)/eg2hcore9d.c $(SRCDIR)/eg2hcore10d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hdiagf.o: $(SRCDIR)/eg2hdiagf.c $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hdiagd.o: $(SRCDIR)/eg2hdiagd.c $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hextf.o: $(SRCDIR)/eg2hextf.c $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hextd.o: $(SRCDIR)/eg2hextd.c $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hconstrf.o: $(SRCDIR)/eg2hconstrf.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hconstrd.o: $(SRCDIR)/eg2hconstrd.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hextconstrf.o: $(SRCDIR)/eg2hextconstrf.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hextconstrd.o: $(SRCDIR)/eg2hextconstrd.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hfunclf.o: $(SRCDIR)/eg2hfunclf.c $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hfuncld.o: $(SRCDIR)/eg2hfuncld.c $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hgraphf.o: $(SRCDIR)/eg2hgraphf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg2holef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hgraphd.o: $(SRCDIR)/eg2hgraphd.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2herrorf.o: $(SRCDIR)/eg2herrorf.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2herrord.o: $(SRCDIR)/eg2herrord.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin1f.o: $(SRCDIR)/eg2hnonlin1f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin1d.o: $(SRCDIR)/eg2hnonlin1d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin2f.o: $(SRCDIR)/eg2hnonlin2f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin2d.o: $(SRCDIR)/eg2hnonlin2d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin3f.o: $(SRCDIR)/eg2hnonlin3f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin3d.o: $(SRCDIR)/eg2hnonlin3d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin4f.o: $(SRCDIR)/eg2hnonlin4f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hnonlin4d.o: $(SRCDIR)/eg2hnonlin4d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore1f.o: $(SRCDIR)/eg2hsplcore1f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore2f.o: $(SRCDIR)/eg2hsplcore2f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore3f.o: $(SRCDIR)/eg2hsplcore2f.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplconstrf.o: $(SRCDIR)/eg2hsplconstrf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg2holef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore1d.o: $(SRCDIR)/eg2hsplcore1d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore2d.o: $(SRCDIR)/eg2hsplcore2d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplcore3d.o: $(SRCDIR)/eg2hsplcore3d.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplconstrd.o: $(SRCDIR)/eg2hsplconstrd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg2holed.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplgraphf.o: $(SRCDIR)/eg2hsplgraphf.c \
  $(INCDIR)/eg2holef.h $(SRCDIR)/eghprivatef.h $(SRCDIR)/eg2hprivatef.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplgraphd.o: $(SRCDIR)/eg2hsplgraphd.c \
  $(INCDIR)/eg2holed.h $(SRCDIR)/eghprivated.h $(SRCDIR)/eg2hprivated.h \
  $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hpcurvesf.o: $(SRCDIR)/eg2hpcurvesf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hpcurvesd.o: $(SRCDIR)/eg2hpcurvesd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg2hprivated.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplpcurvesf.o: $(SRCDIR)/eg2hsplpcurvesf.c \
  $(SRCDIR)/eghprivatef.h $(INCDIR)/eg2holef.h \
  $(SRCDIR)/eg2hprivatef.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg2hsplpcurvesd.o: $(SRCDIR)/eg2hsplpcurvesd.c \
  $(SRCDIR)/eghprivated.h $(INCDIR)/eg2holed.h \
  $(SRCDIR)/eg2hprivated.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/eg1hsymmatrix2.o: $(SRCDIR)/eg1hsymmatrix2.c $(INCDIR)/eg1hole2.h \
  $(SRCDIR)/eg1hprivated.h $(SRCDIR)/eg1herror.h

$(SRCDIR)/eg2hsymmatrix2.o: $(SRCDIR)/eg2hsymmatrix2.c $(INCDIR)/eg2hole2.h \
  $(SRCDIR)/eg2hprivate2.h $(SRCDIR)/eg2herror.h

$(SRCDIR)/egheigenmeshf.o: $(SRCDIR)/egheigenmeshf.c $(INCDIR)/egholef.h

$(SRCDIR)/egheigenvalf.o: $(SRCDIR)/egheigenvalf.c $(INCDIR)/egholef.h

$(SRCDIR)/egheigenmeshd.o: $(SRCDIR)/egheigenmeshd.c $(INCDIR)/egholef.h

$(SRCDIR)/egheigenvald.o: $(SRCDIR)/egheigenvald.c $(INCDIR)/egholed.h


clean:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

mrproper:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	rm -f $(TARGET)
