
HOME = /home/przemek
INCDIR = -I../../include
LIBDIR = -L../../lib
LIBS =-lcamera -lbsfile -lg2blending -lbsmesh -leghole -lmultibs \
  -lpkgeom -lpknum -lpsout  -lpkvaria -lm

include ../../options.mak

TARGET = optblmesh

default: $(TARGET)

.c.o:
	$(CC) -c $(CFLAGS) $(INCDIR) $< -o $*.o

optblmesh: optblmesh.c
	$(CC) $(CFLAGS) optblmesh.c $(INCDIR) $(LIBDIR) $(LIBS) -o optblmesh
	./optblmesh -i trojnik3.bs -o trojnik3f.bs
	./optblmesh -s -i trojnik3f.bs -o trojnik3s.bs

# To test all optimization procedures built in this program, type make all.
# Be warned that these tests may take some time.

all: optblmesh
	./optblmesh -e -b 4 -i trojnik5.bs -o trojnik5f.bs
	./optblmesh -e -b 12 -i trojnik6.bs -c trojnik3.bs -o trojnik6f.bs
	./optblmesh -e -i trojnik6.bs -c trojnik3.bs -o trojnik6f.bs
	./optblmesh -e -s -i trojnik6f.bs -c trojnik3.bs -o trojnik6s.bs

clean:
	rm -f *.o *f.bs *s.bs

mrproper:
	rm -f *.o $(TARGET) *.ps *.txt *f.bs *s.bs

