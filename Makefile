CFLAGS ?= -Wall -O3
CFLAGS += -std=c99
LIBS = cdist.o dbscan.o dbparse.o dist.o fbseek.o filebuff.o fsacmp.o fsacmpthrd.o hashmapstr.o hashmapstrindex.o ltdmatrix.o ltdmatrixthrd.o matcmp.o matparse.o matrix.o merge.o meth.o methparse.o nj.o nwck.o nwck2phy.o pherror.o phy.o qseqs.o rarify.o resparse.o seqparse.o seq2fasta.o stdnuc.o stdstat.o str.o tmp.o tree.o trim.o ulist.o union.o unionparse.o vector.o
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

cdist.o: cdist.h filebuff.h matrix.h meth.h pherror.h phy.h seqparse.h
dbscan.o: dbscan.h filebuff.h matrix.h pherror.h phy.h qseqs.h tmp.h
dbparse.o: dbparse.h pherror.h qseqs.h
dist.o: dist.h ltdmatrix.h matcmp.h matrix.h meth.h methparse.h pherror.h phy.h unionparse.h
fbseek.o: fbseek.h filebuff.h pherror.h
filebuff.o: filebuff.h pherror.h
fsacmp.o: fsacmp.h matrix.h pherror.h qseqs.h threader.h
fsacmpthrd.o: fsacmpthrd.h fsacmp.h matrix.h pherror.h threader.h
hashmapstr.o: hashmapstr.h pherror.h
hashmapstrindex.o: hashmapstrindex.h hashmapstr.h pherror.h
ltdmatrix.o: ltdmatrix.h filebuff.h ltdmatrix.h matcmp.h matparse.h pherror.h
ltdmatrixthrd.o: ltdmatrixthrd.h filebuff.h matcmp.h matparse.h matrix.h pherror.h threader.h
matcmp.o: matcmp.h filebuff.h matparse.h stdstat.h
matparse.o: matparse.h filebuff.h pherror.h qseqs.h
matrix.o: matrix.h pherror.h tmp.h
merge.o: merge.h filebuff.h hashmapstr.h hashmapstrindex.h matrix.h phy.h qseqs.h ulist.h
meth.o: meth.h pherror.h
methparse.o: methparse.h filebuff.h meth.h pherror.h qseqs.h
nj.o: nj.h matrix.h pherror.h threader.h vector.h
nwck.o: nwck.h filebuff.h qseqs.h pherror.h
nwck2phy.o: nwck2phy.h filebuff.h matrix.h nwck.h pherror.h phy.h qseqs.h
pherror.o: pherror.h
phy.o: phy.h matrix.h
qseqs.o: qseqs.h pherror.h
rarify.o: rarify.h filebuff.h matparse.h pherror.h
resparse.o: resparse.h filebuff.h pherror.h qseqs.h
seqparse.o: seqparse.h filebuff.h pherror.h qseqs.h
seq2fasta.o: seq2fasta.h dbparse.h pherror.h qseqs.h stdnuc.h
stdnuc.o: stdnuc.h
stdstat.o: stdstat.h
str.o: str.h
tmp.o: tmp.h pherror.h threader.h
tree.o: tree.h filebuff.h matrix.h nj.h pherror.h phy.h qseqs.h tmp.h vector.h
trim.o: trim.h filebuff.h fsacmp.h matrix.h meth.h methparse.h pherror.h phy.h seqparse.h
ulist.o: ulist.h pherror.h
union.o: union.h filebuff.h hashmapstr.h pherror.h resparse.h
unionparse.o: unionparse.h filebuff.h pherror.h
vector.o: vector.h pherror.h
