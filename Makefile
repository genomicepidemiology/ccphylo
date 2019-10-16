CFLAGS = -Wall -O3 -std=c99
LIBS = dbparse.o dist.o fbseek.o filebuff.o hashmapstr.o ltdmatrix.o matcmp.o matparse.o matrix.o nj.o nwck.o pherror.o phy.o qseqs.o resparse.o stdstat.o tree.o union.o vector.o
PROGS = ccphylo

.c .o:
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROGS)

ccphylo: main.c libccphylo.a
	$(CC) $(CFLAGS) -o $@ main.c libccphylo.a -lm -lpthread -lz

libccphylo.a: $(LIBS)
	$(AR) -csru $@ $(LIBS)

clean:
	$(RM) $(LIBS) $(PROGS) libccphylo.a

dbparse.o: dbparse.h pherror.h qseqs.h
dist.o: dist.h ltdmatrix.h matcmp.h
fbseek.o: fbseek.h filebuff.h pherror.h
filebuff.o: filebuff.h pherror.h
hashmapstr.o: hashmapstr.h pherror.h
ltdmatrix.o: ltdmatrix.h filebuff.h ltdmatrix.h matcmp.h matparse.h pherror.h
matcmp.o: matcmp.h filebuff.h matparse.h stdstat.h
matparse.o: matparse.h filebuff.h pherror.h qseqs.h
matrix.o: matrix.h pherror.h
nj.o: nj.h matrix.h pherror.h vector.h
nwck.o: nwck.h qseqs.h pherror.h
pherror.o: pherror.h
phy.o: phy.h matrix.h
qseqs.o: qseqs.h pherror.h
resparse.o: resparse.h filebuff.h qseqs.h pherror.h
stdstat.o: stdstat.h
union.o: union.h filebuff.h hashmapstr.h pherror.h resparse.h
tree.o: tree.h filebuff.h matrix.h nj.h pherror.h phy.h qseqs.h vector.h
vector.o: vector.h pherror.h
