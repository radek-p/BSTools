
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
CPICDIR = ../cpict
PICDIR = ../cpict
INCDIR = ../include
LIBDIR = ../lib

TEXSRC = bsfile.tex bsmesh.tex bstools.tex camera.tex demo.tex eghole.tex \
  experim.tex g1blending.tex g2blending.tex header.tex intro.tex \
  multibs1.tex multibs2.tex multibs3.tex multibs4.tex multibs5.tex \
  multibs6.tex multibs7.tex multibs.tex pkgeom.tex pknum.tex pkvaria.tex \
  projects.tex psout.tex raybez.tex xgedit.tex

default: bstoolsA4.ps bstoolsB5.ps

bstoolsA4.ps: bstools.dvi
	dvips -t a4 -o bstools.ps -O 1.5cm,2.2cm bstools
	psnup -pa4 -2 -b-90 bstools.ps bstoolsA4.ps
	rm bstools.ps
	ps2pdf -sPAPERSIZE=a4 bstoolsA4.ps bstoolsA4.pdf

bstoolsB5.ps: bstools.dvi
	dvips -t a4 -o bstoolsB5.ps -O 1.84cm,-0.6cm bstools 
	ps2pdf -sPAPERSIZE=a4 bstoolsB5.ps bstoolsB5.pdf

bstools.dvi: $(TEXSRC) bstools.aux memory.ps
	latex bstools

bstools.aux: $(TEXSRC) memory.ps
	latex bstools

memory.ps:
	(cd ../cpict;$(MAKE))
	rm -f $(SRCDIR)/*.ps
	ln -s $(PICDIR)/memory.ps $(SRCDIR)/memory.ps
	ln -s $(PICDIR)/memshift.ps $(SRCDIR)/memshift.ps
	ln -s $(PICDIR)/camera.ps $(SRCDIR)/camera.ps
	ln -s $(PICDIR)/stereo.ps $(SRCDIR)/stereo.ps
	ln -s $(PICDIR)/psout.ps $(SRCDIR)/psout.ps
	ln -s $(PICDIR)/bsder123.ps $(SRCDIR)/bsder123.ps
	ln -s $(PICDIR)/patchder.ps $(SRCDIR)/patchder.ps
	ln -s $(PICDIR)/patchpdir.ps $(SRCDIR)/patchpdir.ps
	ln -s $(PICDIR)/knotrem.ps $(SRCDIR)/knotrem.ps
	ln -s $(PICDIR)/oslo.ps $(SRCDIR)/oslo.ps
	ln -s $(PICDIR)/bspbez.ps $(SRCDIR)/bspbez.ps
	ln -s $(PICDIR)/degel.ps $(SRCDIR)/degel.ps
	ln -s $(PICDIR)/bspdegel.ps $(SRCDIR)/bspdegel.ps
	ln -s $(PICDIR)/splmult.ps $(SRCDIR)/splmult.ps
	ln -s $(PICDIR)/nvtest1.ps $(SRCDIR)/nvtest1.ps
	ln -s $(PICDIR)/nvtest2.ps $(SRCDIR)/nvtest2.ps
	ln -s $(PICDIR)/chknots.ps $(SRCDIR)/chknots.ps
	ln -s $(PICDIR)/intbsc.ps $(SRCDIR)/intbsc.ps
	ln -s $(PICDIR)/bsapprox.ps $(SRCDIR)/bsapprox.ps
	ln -s $(PICDIR)/rasterbc.ps $(SRCDIR)/rasterbc.ps
	ln -s $(PICDIR)/rasterbs.ps $(SRCDIR)/rasterbs.ps
	ln -s $(PICDIR)/trimpatch.ps $(SRCDIR)/trimpatch.ps
	ln -s $(PICDIR)/raybezp.ps $(SRCDIR)/raybezp.ps
	ln -s $(PICDIR)/g1patches1.ps $(SRCDIR)/g1patches1.ps
	ln -s $(PICDIR)/g1patches2.ps $(SRCDIR)/g1patches2.ps
	ln -s $(PICDIR)/g2hdomain.ps $(SRCDIR)/g2hdomain.ps
	ln -s $(PICDIR)/degred1.ps $(SRCDIR)/degred1.ps
	ln -s $(PICDIR)/degred2.ps $(SRCDIR)/degred2.ps
	ln -s $(PICDIR)/addspl.ps $(SRCDIR)/addspl.ps
	ln -s $(PICDIR)/bsplane.ps $(SRCDIR)/bsplane.ps
	ln -s $(PICDIR)/bsmesh.ps $(SRCDIR)/bsmesh.ps

clean:
	(cd ../cpict;$(MAKE) clean)
	mv bstoolsB5.ps bstoolsB5
	mv bstoolsA4.ps bstoolsA4
	rm -f *.ps
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *.toc
	mv bstoolsB5 bstoolsB5.ps
	mv bstoolsA4 bstoolsA4.ps

mrproper:
	(cd ../cpict;$(MAKE) mrproper)
	rm -f *.ps
	rm -f *.pdf
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *.toc

