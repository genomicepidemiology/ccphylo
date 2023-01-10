/* Philip T.L.C. Clausen Apr 2021 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <stdio.h>
#include "bytescale.h"
#include "cmdline.h"
#include "dbscan.h"
#include "filebuff.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"

int dbscan(Matrix *D, int *N, int *C, double maxDist, int minN) {
	
	int i, j, c, N_i, nClust, *Nptr, *Cptr;
	double d, **Dmat, *Dptr;
	float **Dfmat, *Dfptr;
	short unsigned **Dsmat, *Dsptr;
	unsigned char **Dbmat, *Dbptr;
	
	/*
	D	Distance matrix, given
	N	Number of neighbors, output
	C	Cluster number of each node, output
	*/
	
	/* get number of neighbors pr. node */
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	i = -1;
	while(++i < D->n) {
		N_i = 0;
		Nptr = N - 1;
		j = -1;
		while(++j < i) {
			d = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(d <= maxDist) {
				++N_i;
				++*++Nptr;
			} else {
				++Nptr;
			}
		}
		N[i] = N_i;
		C[i] = i;
	}
	
	/* assign cluster numbers to nodes */
	nClust = 0;
	Dmat = 0;
	Dfmat = 0;
	Dsmat = 0;
	Dbmat = 0;
	if(Dptr) {
		Dmat = D->mat;
	} else if(Dfptr) {
		Dfmat = D->fmat;
	} else if(Dsptr) {
		Dsmat = D->smat;
	} else {
		Dbmat = D->bmat;
	}
	Nptr = N - 1;
	Cptr = C - 1;
	c = 0;
	i = -1;
	while(++i < D->n) {
		if(minN <= *++Nptr) {
			if(Dmat) {
				Dptr = Dmat[i] - 1;
			} else if(Dfmat) {
				Dfptr = Dfmat[i] - 1;
			} else if(Dsmat) {
				Dsptr = Dsmat[i] - 1;
			} else {
				Dbptr = Dbmat[i] - 1;
			}
			
			c = i;
			j = -1;
			while(++j < c) {
				d = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
				if(d <= maxDist) {
					/* assign node i to the same cluster as node j */
					c = C[j];
				}
			}
			if(i != c) {
				*++Cptr = c;
			} else {
				++nClust;
				++Cptr;
			}
		} else if(*Nptr) {
			if(Dmat) {
				Dptr = Dmat[i] - 1;
			} else if(Dfmat) {
				Dfptr = Dfmat[i] - 1;
			} else if(Dsmat) {
				Dsptr = Dsmat[i] - 1;
			} else {
				Dbptr = Dbmat[i] - 1;
			}
			
			N_i = *Nptr;
			c = i;
			j = -1;
			while(++j < c) {
				d = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
				if(d <= maxDist) {
					if(minN <= N[j]) {
						/* assign node i to the same cluster as node j */
						c = C[j];
					} else if(!--N_i) {
						/* no more neighbors */
						j = c;
					}
				}
			}
			if(i != c) {
				*++Cptr = c;
			} else {
				++nClust;
				++Cptr;
			}
		} else {
			++nClust;
			++Cptr;
		}
	}
	
	/* return number of clusters */
	return nClust;
}

void print_dbscan(Qseqs **names, int *N, int *C, int Dn, int nClust, double maxDist, int minN, FILE *out) {
	
	/* print header */
	fprintf(out, "## %d\t%d\t%lf\t%d\n", Dn, nClust, maxDist, minN);
	fprintf(out, "#%s\t%s\t%s\n", "Sample", "Neighbors", "Cluster");
	
	/* name, #Neighbors, Cluster */
	--names;
	--N;
	--C;
	++Dn;
	while(--Dn) {
		fprintf(out, "%s\t%d\t%d\n", (*++names)->seq, *++N, *++C);
	}
}

void make_dbscan(char *inputfilename, char *outputfilename, double maxDist, int minN, char sep, char quotes) {
	
	int i, size, nClust, *N, *C;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	size = 32;
	D = ltdMatrix_init(size);
	N = smalloc(2 * size * sizeof(unsigned));
	C = N + size;
	names = smalloc(size * sizeof(Qseqs *));
	names += size;
	i = size + 1;
	while(--i) {
		*--names = setQseqs(64);
	}
	
	/* set */
	openAndDetermine(infile, inputfilename);
	
	/* generate dbscans */
	while((names = loadPhy(D, names, header, infile, sep, quotes)) && D->n) {
		/* realloc */
		if(size < D->n) {
			size = D->n;
			N = smalloc(2 * size * sizeof(unsigned));
			C = N + size;
		}
		
		/* dbscan */
		nClust = dbscan(D, N, C, maxDist, minN);
		
		/* output */
		if(header->len) {
			fprintf(outfile, "#%s\n", header->seq);
		}
		print_dbscan(names, N, C, D->n, nClust, maxDist, minN, outfile);
	}
	
	/* clean */
	fclose(outfile);
	closeFileBuff(infile);
	free(N);
	Matrix_destroy(D);
	free(names);
	destroyFileBuff(infile);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo make a DBSCAN given a set of phylip distance matrices.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'q', "quotes", "Quote taxa", "\\0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'N', "min_neighbors", "Minimum neighbors", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'e', "max_distance", "Maximum distance", "10.0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'p', "float_precision", "Float precision on distance matrix", "double");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 's', "short_precision", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "byte_precision", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'H', "mmap", "Allocate matrix on the disk", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tmp", "Set directory for temporary files", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	n	N	min_neighbors
	d	e	max_distance
	fp	p	float_precision
	sp	s	short_precision
	bp	b	byte_precision
	mm	H	mmap
	tmp	T	tmp
	h	h	help
	*/
}

int main_dbscan(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, minNum, size, len, offset;
	double maxDist;
	char **Arg, *arg, *inputfilename, *outputfilename, *tmp, opt, sep, quotes;
	
	/* set defaults */
	size = sizeof(double);
	minNum = 1;
	maxDist = 10.0;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	tmp = 0;
	sep = '\t';
	quotes = '\0';
	
	/* parse cmd-line */
	args = argc - 1;
	Arg = argv;
	if(args && **++Arg == '-') {
		len = 1;
		--Arg;
	} else {
		len = 0;
	}
	while(args && len) {
		arg = *++Arg;
		if(*arg++ == '-') {
			if(*arg == '-') {
				/* check if argument is included */
				len = getOptArg(++arg);
				offset = 2 + (arg[len] ? 1 : 0);
				
				/* long option */
				if(*arg == 0) {
					/* terminate cmd-line */
					++Arg;
				} else if(cmdcmp(arg, "input") == 0) {
					inputfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "quotes") == 0) {
					quotes = getcArgDie(&Arg, &args, len + offset, "quotes");
				} else if(cmdcmp(arg, "min_neighbors") == 0) {
					minNum = getNumArg(&Arg, &args, len + offset, "min_neighbors");
				} else if(cmdcmp(arg, "max_distance") == 0) {
					maxDist = getdArg(&Arg, &args, len + offset, "max_distance");
				} else if(cmdcmp(arg, "float_precision") == 0) {
					size = sizeof(float);
				} else if(cmdcmp(arg, "short_precision") == 0) {
					size = sizeof(short unsigned);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "short_precision");
				} else if(cmdcmp(arg, "byte_precision") == 0) {
					size = sizeof(unsigned char);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "byte_precision");
				} else if(cmdcmp(arg, "mmap") == 0) {
					ltdMatrix_init = &ltdMatrixMinit;
				} else if(cmdcmp(arg, "tmp") == 0) {
					tmp = getArgDie(&Arg, &args, len + offset, "tmp");
				} else if(cmdcmp(arg, "help") == 0) {
					return helpMessage(stdout);
				} else {
					unknArg(arg - 2);
				}
			} else {
				/* multiple option */
				len = 1;
				opt = *arg;
				while(opt && (opt = *arg++)) {
					++len;
					if(opt == 'i') {
						inputfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'q') {
						quotes = getcArgDie(&Arg, &args, len, "q");
						opt = 0;
					} else if(opt == 'N') {
						minNum = getNumArg(&Arg, &args, len, "N");
						opt = 0;
					} else if(opt == 'e') {
						maxDist = getdArg(&Arg, &args, len, "e");
						opt = 0;
					} else if(opt == 'p') {
						size = sizeof(float);
					} else if(opt == 's') {
						size = sizeof(short unsigned);
						ByteScale = getdDefArg(&Arg, &args, len, ByteScale, "p");
						opt = 0;
					} else if(opt == 'b') {
						size = sizeof(unsigned char);
						ByteScale = getdDefArg(&Arg, &args, len, ByteScale, "b");
						opt = 0;
					} else if(opt == 'H') {
						ltdMatrix_init = &ltdMatrixMinit;
						opt = 0;
					} else if(opt == 'T') {
						tmp = getArgDie(&Arg, &args, len, "T");
						opt = 0;
					} else if(opt == 'h') {
						return helpMessage(stdout);
					} else {
						*arg = 0;
						unknArg(arg - 1);
					}
				}
			}
		} else {
			/* terminate cmd-line */
			--arg;
			++args;
			len = 0;
		}
		--args;
	}
	
	
	/* non-options */
	if(args) {
		inputfilename = *Arg;
		if(--args) {
			nonOptError();
		}
	}
	
	/* tmp dir */
	if(tmp) {
		tmpF(tmp);
	}
	
	/* set precision */
	ltdMatrixInit(-size);
	ltdMatrixMinit(-size);
	
	/* make tree */
	make_dbscan(inputfilename, outputfilename, maxDist, minNum, sep, quotes);
	
	return 0;
}
