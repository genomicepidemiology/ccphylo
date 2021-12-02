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
	double d, **Dmat, *Dptr;
	float **Dfmat, *Dfptr;
	short unsigned **Dsmat, *Dsptr;
	unsigned char **Dbmat, *Dbptr;
	
	/*
	D	Distance matrix, given
	N	Number of neighbours, output
	C	Cluster number of each node, output
	*/
	
	/* get number of neighbours pr. node */
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-sp", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_dbscan(int argc, char *argv[]) {
	
	int args, minNum, size;
	double maxDist;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	
	/* set defaults */
	size = sizeof(double);
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
			} else if(strcmp(arg, "fp") == 0) {
				size = sizeof(float);
			} else if(strcmp(arg, "sp") == 0) {
				if(++args < argc && argv[args][0] != '-') {
					ByteScale = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || ByteScale == 0) {
						invaArg("\"-sp\"");
					}
				} else {
					--args;
				}
				size = sizeof(short unsigned);
			} else if(strcmp(arg, "bp") == 0) {
				if(++args < argc && argv[args][0] != '-') {
					ByteScale = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || ByteScale == 0) {
						invaArg("\"-bp\"");
					}
				} else {
					--args;
				}
				size = sizeof(unsigned char);
			} else if(strcmp(arg, "mm") == 0) {
				ltdMatrix_init = &ltdMatrixMinit;
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
	
	/* set precision */
	ltdMatrixInit(-size);
	ltdMatrixMinit(-size);
	
	/* make tree */
	make_dbscan(inputfilename, outputfilename, maxDist, minNum);
	
	return 0;
}
