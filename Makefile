CFLAGS ?= -Wall -O3
CFLAGS += -std=c99
LIBS = bytescale.o cdist.o cmdline.o dat.o datclust.o dbscan.o dbparse.o distcmp.o dist.o dnj.o fbseek.o filebuff.o fsacmp.o fsacmpthrd.o fullphy.o hashmapstr.o hashmapstrindex.o hclust.o jobs.o ltdmatrix.o ltdmatrixthrd.o machines.o makespan.o matcmp.o matparse.o matrix.o merge.o meth.o methparse.o nj.o nwck.o nwck2phy.o pherror.o phy.o qseqs.o rarify.o resparse.o seqparse.o seq2fasta.o stdnuc.o stdstat.o str.o tabusearch.o tmp.o tree.o trim.o tsv.o tsv2nwck.o tsv2phy.o ulist.o union.o unionparse.o vector.o
PROGS = ccphylo

.c .o:
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROGS)

ccphylo: main.c libccphylo.a
	$(CC) $(CFLAGS) -o $@ main.c libccphylo.a -lm -lpthread -lz

libccphylo.a: $(LIBS)
	$(AR) -csr $@ $(LIBS)

clean:
	$(RM) $(LIBS) $(PROGS) libccphylo.a

bytescale.o: bytescale.h
cdist.o: cdist.h filebuff.h matrix.h meth.h pherror.h phy.h seqparse.h
cmdline.o: cmdline.h
dat.c: dat.h pherror.h tmp.h
datclust.o: datclust.h dat.h distcmp.h hclust.h nwck.h qseqs.h vector.h
dbscan.o: dbscan.h bytescale.h cmdline.h filebuff.h matrix.h pherror.h phy.h qseqs.h tmp.h
dbparse.o: dbparse.h pherror.h qseqs.h
dist.o: dist.h cmdline.h ltdmatrix.h matcmp.h matrix.h meth.h methparse.h pherror.h phy.h unionparse.h
distcmp.o: distcmp.h bytescale.h
dnj.o: dnj.h hclust.h matrix.h nj.h nwck.h qseqs.h pherror.h str.h threader.h vector.h
fbseek.o: fbseek.h filebuff.h pherror.h
filebuff.o: filebuff.h pherror.h
fsacmp.o: fsacmp.h matrix.h pherror.h qseqs.h threader.h
fsacmpthrd.o: fsacmpthrd.h fsacmp.h matrix.h pherror.h threader.h
fullphy.o: fullphy.h bytescale.h cmdline.h filebuff.h matrix.h pherror.h phy.h qseqs.h tmp.h
hashmapstr.o: hashmapstr.h pherror.h
hashmapstrindex.o: hashmapstrindex.h hashmapstr.h pherror.h
hclust.o: hclust.h matrix.h nj.h nwck.h qseqs.h pherror.h str.h vector.h
jobs.o: jobs.h filebuff.h pherror.h
ltdmatrix.o: ltdmatrix.h filebuff.h ltdmatrix.h matcmp.h matparse.h pherror.h
ltdmatrixthrd.o: ltdmatrixthrd.h filebuff.h matcmp.h matparse.h matrix.h pherror.h threader.h
machines.o: machines.h jobs.h pherror.h
makespan.o: makespan.h filebuff.h jobs.h machines.h pherror.h tabusearch.h tsv.h
matcmp.o: matcmp.h filebuff.h matparse.h stdstat.h
matparse.o: matparse.h filebuff.h pherror.h qseqs.h
matrix.o: matrix.h pherror.h tmp.h
merge.o: merge.h cmdline.h filebuff.h hashmapstr.h hashmapstrindex.h matrix.h phy.h qseqs.h ulist.h
meth.o: meth.h pherror.h
methparse.o: methparse.h filebuff.h meth.h pherror.h qseqs.h
nj.o: nj.h matrix.h pherror.h threader.h vector.h
nwck.o: nwck.h filebuff.h qseqs.h pherror.h
nwck2phy.o: nwck2phy.h cmdline.h filebuff.h matrix.h nwck.h pherror.h phy.h qseqs.h
pherror.o: pherror.h
phy.o: phy.h matrix.h
qseqs.o: qseqs.h pherror.h
rarify.o: rarify.h cmdline.h filebuff.h matparse.h pherror.h
resparse.o: resparse.h filebuff.h pherror.h qseqs.h
seqparse.o: seqparse.h filebuff.h pherror.h qseqs.h
seq2fasta.o: seq2fasta.h dbparse.h pherror.h qseqs.h stdnuc.h
stdnuc.o: stdnuc.h
stdstat.o: stdstat.h
str.o: str.h
tabusearch.o: tabusearch.h jobs.h machines.h
tmp.o: tmp.h pherror.h threader.h
tree.o: tree.h bytescale.h cmdline.h dnj.h filebuff.h hclust.h matrix.h nj.h pherror.h phy.h qseqs.h tmp.h vector.h
trim.o: trim.h cmdline.h filebuff.h fsacmp.h matrix.h meth.h methparse.h pherror.h phy.h seqparse.h
tsv.o: tsv.h bytescale.h dat.h filebuff.h jobs.h pherror.h stdstat.h
tsv2nwck.o: tsv2nwck.h bytescale.h dat.h datclust.h distcmp.h filebuff.h hclust.h nwck.h pherror.h qseqs.h tmp.h tsv.h vector.h
tsv2phy.o: tsv2phy.h bytescale.h cmdline.h dat.h distcmp.h filebuff.h pherror.h tmp.h tsv.h
ulist.o: ulist.h pherror.h
union.o: union.h cmdline.h filebuff.h hashmapstr.h pherror.h resparse.h
unionparse.o: unionparse.h filebuff.h pherror.h
vector.o: vector.h pherror.h
