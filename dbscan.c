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
#include "dbscan.h"
#include "filebuff.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

int dbscan(Matrix *D, int *N, int *C, double maxDist, int minN) {
	
	int i, j, c, N_i, nClust, *Nptr, *Cptr;
	double **Dmat, *Dptr;
	
	/*
	D	Distance matrix, given
	N	Number of neighbours, output
	C	Cluster number of each node, output
	*/
	
	/* get number of neighbours pr. node */
	Dptr = *(D->mat) - 1;
	i = -1;
	while(++i < D->n) {
		N_i = 0;
		Nptr = N - 1;
		j = -1;
		while(++j < i) {
			if(*++Dptr <= maxDist) {
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
	nClust = 1;
	Dmat = D->mat;
	Dptr = *(D->mat) - 1;
	Nptr = N - 1;
	Cptr = C - 1;
	i = -1;
	while(++i < D->n) {
		c = i;
		if(minN <= *++Nptr) {
			j = -1;
			while(++j < c) {
				if(Dmat[i][j] <= maxDist) {
					/* assign node i to the same cluster as node j */
					c = C[j];
				}
			}
		} else {
			++nClust;
		}
		*++Cptr = c;
	}
	
	/* return number of clusters */
	return nClust;
}

void print_dbscan(Qseqs **names, int *N, int *C, int Dn, int nClust, double maxDist, int minN, FILE *out) {
	
	/* print header */
	fprintf(out, "# %d\t%d\t%lf\t%d\n", Dn, nClust, maxDist, minN);
	
	/* name, #Neighbours, Cluster */
	--names;
	--N;
	--C;
	++Dn;
	while(--Dn) {
		fprintf(out, "%s\t%d\t%d\n", (*++names)->seq, *++N, *++C);
	}
}

void make_dbscan(char *inputfilename, char *outputfilename, double maxDist, int minN) {
	
	int i, size, nClust, *N, *C;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) ? stdout : sfopen(outputfilename, "wb");
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
	while((names = loadPhy(D, names, header, infile)) && D->n) {
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file.", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-n", "Minimum neighbours", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Maximum distance", "10.0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_dbscan(int argc, char *argv[]) {
	
	int args, minNum;
	double maxDist;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	
	/* set defaults */
	minNum = 1;
	maxDist = 10.0;
	inputfilename = "--";
	outputfilename = "--";
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					inputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "n") == 0) {
				if(++args < argc) {
					minNum = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-n\"");
					}
				} else {
					missArg("\"-n\"");
				}
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					maxDist = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || maxDist < 0) {
						invaArg("\"-d\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "tmp") == 0) {
				if(++args < argc) {
					if(argv[args][0] != '-') {
						tmpF(argv[args]);
					} else {
						invaArg("\"-tmp\"");
					}
				}
			} else if(strcmp(arg, "h") == 0) {
				return helpMessage(stdout);
			} else {
				fprintf(stderr, "Unknown option:%s\n", arg - 1);
				return helpMessage(stderr);
			}
		} else {
			fprintf(stderr, "Unknown argument:%s\n", arg - 1);
			return helpMessage(stderr);
		}
		++args;
	}
	
	/* make tree */
	make_dbscan(inputfilename, outputfilename, maxDist, minNum);
	
	return 0;
}
