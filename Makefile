CFLAGS ?= -Wall -O3
CFLAGS += -std=c99
LIBS = dbparse.o dist.o fbseek.o filebuff.o hashmapstr.o hashmapstrindex.o ltdmatrix.o matcmp.o matparse.o matrix.o merge.o nj.o nwck.o nwck2phy.o pherror.o phy.o qseqs.o rarify.o resparse.o stdstat.o str.o tree.o ulist.o union.o unionparse.o vector.o
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
dist.o: dist.h ltdmatrix.h matcmp.h matrix.h pherror.h phy.h unionparse.h
fbseek.o: fbseek.h filebuff.h pherror.h
filebuff.o: filebuff.h pherror.h
hashmapstr.o: hashmapstr.h pherror.h
hashmapstrindex.o: hashmapstrindex.h hashmapstr.h pherror.h
ltdmatrix.o: ltdmatrix.h filebuff.h ltdmatrix.h matcmp.h matparse.h pherror.h
matcmp.o: matcmp.h filebuff.h matparse.h stdstat.h
matparse.o: matparse.h filebuff.h pherror.h qseqs.h
matrix.o: matrix.h pherror.h
merge.o: merge.h filebuff.h hashmapstr.h hashmapstrindex.h matrix.h phy.h qseqs.h ulist.h
nj.o: nj.h matrix.h pherror.h vector.h
nwck.o: nwck.h qseqs.h pherror.h
nwck2phy.o: nwck2phy.h
pherror.o: pherror.h
phy.o: phy.h matrix.h
qseqs.o: qseqs.h pherror.h
rarify.o: rarify.h
resparse.o: resparse.h filebuff.h qseqs.h pherror.h
stdstat.o: stdstat.h
str.o: str.h
tree.o: tree.h filebuff.h matrix.h nj.h pherror.h phy.h qseqs.h vector.h
ulist.o: ulist.h pherror.h
union.o: union.h filebuff.h hashmapstr.h pherror.h resparse.h
unionparse.o: unionparse.h filebuff.h pherror.h
vector.o: vector.h pherror.h

