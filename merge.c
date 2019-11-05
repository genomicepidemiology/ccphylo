/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

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
#include "filebuff.h"
#include "hashmapstr.h"
#include "hashmapstrindex.h"
#include "matrix.h"
#include "merge.h"
#include "phy.h"
#include "qseqs.h"
#include "ulist.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static void syncMatrices(HashMapStr *names_index, uList *dist_index, Qseqs **names, unsigned n, Matrix *dist, Matrix *num) {
	
	int i, index;
	
	dist_index->n = 0; /* clear dist_index */
	for(i = 0; i < n; ++i) {
		if((index = HashMapStrindex_add(names_index, names[i]->seq)) < 0) {
			index = ltdMatrix_add(dist);
			ltdMatrix_add(num);
		}
		uList_push(dist_index, index);
	}
}

static void normalize_ltdMatrix(Matrix *D, Matrix *N) {
	
	int i;
	double *Dptr, *Nptr;
	
	/* normalize new distance matrix */
	Dptr = *(D->mat) - 1;
	Nptr = *(N->mat) - 1;
	i = (D->n * (D->n - 1)) / 2 + 1;
	while(--i) {
		if(*++Nptr != 0) {
			*++Dptr /= *Nptr;
		} else {
			*++Dptr = -1.0;
		}
	}
}

static char ** getOrderedNames(HashMapStr *names_index) {
	
	int i;
	unsigned char **distNames;
	BucketStr *node, *node_next;
	
	distNames = smalloc(names_index->n * sizeof(char *));
	for(i = 0; i < names_index->size && names_index->n; ++i) {
		for(node = names_index->table[i]; node; node = node_next) {
			node_next = node->next;
			distNames[node->num] = node->str;
			free(node);
			--names_index->n;
		}
	}
	HashMapStr_destroy(names_index);
	
	return (char **) distNames;
}

char ** merge(Matrix *dist, Matrix *num, FileBuff *phyfile, FileBuff *numfile) {
	
	unsigned i, j, m, n;
	unsigned *distIndex, *distIndexM, *distIndexN;
	double *Dptr, *Nptr, **distMat, **numMat;
	HashMapStr *names_index;
	Matrix *D, *N;
	Qseqs **names;
	uList *dist_index;
	
	/* init */
	if(num->size < dist->size) {
		ltdMatrix_realloc(num, dist->size);
	}
	D = ltdMatrix_init(dist->size);
	N = ltdMatrix_init(dist->size);
	names_index = HashMapStr_init(128);
	dist_index = uList_init(32);
	names = smalloc(dist->size * sizeof(char *));
	names = loadPhy(dist, 0, 0, phyfile);
	loadPhy(num, names, 0, numfile);
	if(dist->n != num->n) {
		fprintf(stderr, "Distance and included nucleotides does not concur!\n");
		exit(1);
	}
	/* weigh the distances */
	i = (dist->n * (dist->n - 1)) / 2 + 1;
	Dptr = *(dist->mat) - 1;
	Nptr = *(num->mat) - 1;
	while(--i) {
		*++Dptr *= *++Nptr;
	}
	
	/* keep track of names */
	i = -1;
	while(++i < dist->n) {
		HashMapStrindex_add(names_index, names[i]->seq);
	}
	
	/* iterate input */
	while((names = loadPhy(D, names, 0, phyfile)) && D->n) {
		if(!loadPhy(N, names, 0, numfile) || N->n != D->n) {
			fprintf(stderr, "Distance and included nucleotides does not concur!\n");
			exit(1);
		}
		
		/* get index of row w.r.t. merged distance matrix */
		syncMatrices(names_index, dist_index, names, D->n, dist, num);
		
		/* add new matrix to the merged matrix */
		distMat = dist->mat;
		numMat = num->mat;
		Dptr = *(D->mat) - 1;
		Nptr = *(N->mat) - 1;
		distIndex = dist_index->list - 1;
		distIndexM = distIndex;
		i = 0;
		while(i != D->n) {
			m = *++distIndexM;
			distIndexN = distIndex;
			j = ++i;
			while(--j) {
				n = *++distIndexN;
				/* update cell */
				distMat[m][n] += *++Dptr * *++Nptr;
				numMat[m][n] += *Nptr;
			}
		}
	}
	
	/* normalize new distance matrix */
	normalize_ltdMatrix(dist, num);
	fprintf(stderr, "YO\n");
	
	/* clean up */
	i = D->size;
	while(i--) {
		destroyQseqs(names[i]);
	}
	free(names);
	fprintf(stderr, "NO\n");
	Matrix_destroy(D);
	Matrix_destroy(N);
	uList_destroy(dist_index);
	
	/* return names w.r.t. merged matrix */
	return getOrderedNames(names_index);
}

char ** jl_merge(Matrix *dist, Matrix *num, FileBuff *phyfile) {
	
	unsigned i, j, m, n;
	unsigned *distIndex, *distIndexM, *distIndexN;
	double *Dptr, **distMat, **numMat;
	HashMapStr *names_index;
	Matrix *D;
	Qseqs **names;
	uList *dist_index;
	
	/* init */
	D = ltdMatrix_init(dist->size);
	names_index = HashMapStr_init(128);
	dist_index = uList_init(32);
	names = loadPhy(dist, 0, 0, phyfile);
	if(num->size < dist->size) {
		ltdMatrix_realloc(num, dist->size);
	}
	num->n = dist->n;
	i = (num->n * (num->n - 1)) / 2 + 1;
	Dptr = *(num->mat) - 1;
	while(--i) {
		*++Dptr = 1;
	}
	
	/* keep track of names */
	i = -1;
	while(++i < dist->n) {
		HashMapStrindex_add(names_index, names[i]->seq);
	}
	
	/* iterate input */
	while((names = loadPhy(D, names, 0, phyfile)) && D->n) {
		/* get index of row w.r.t. merged distance matrix */
		syncMatrices(names_index, dist_index, names, D->n, dist, num);
		
		/* add new matrix to the merged matrix */
		distMat = dist->mat;
		numMat = num->mat;
		Dptr = *(D->mat) - 1;
		distIndex = dist_index->list - 1;
		distIndexM = distIndex;
		i = 0;
		while(i != D->n) {
			m = *++distIndexM;
			distIndexN = distIndex;
			j = ++i;
			while(--j) {
				n = *++distIndexN;
				/* update cell */
				distMat[m][n] += *++Dptr;
				++numMat[m][n];
			}
		}
	}
	
	/* normalize new distance matrix */
	normalize_ltdMatrix(dist, num);
	
	/* clean up */
	i = D->size;
	while(i--) {
		destroyQseqs(names[i]);
	}
	free(names);
	Matrix_destroy(D);
	uList_destroy(dist_index);
	
	/* return names w.r.t. merged matrix */
	return getOrderedNames(names_index);
}

int merger(char *phyfilename, char *numfilename, char *outphyfilename, char *outnumfilename, unsigned format) {
	
	char **names;
	FILE *outphy, *outnum;
	FileBuff *phyfile, *numfile;
	Matrix *dist, *num;
	
	/* init */
	phyfile = setFileBuff(1048576);
	dist = ltdMatrix_init(64);
	num = ltdMatrix_init(64);
	
	/* merge matrices */
	openAndDetermine(phyfile, phyfilename);
	if(numfilename) {
		numfile = setFileBuff(1048576);
		openAndDetermine(numfile, numfilename);
		names = merge(dist, num, phyfile, numfile);
		closeFileBuff(numfile);
		destroyFileBuff(numfile);
	} else {
		names = jl_merge(dist, num, phyfile);
	}
	closeFileBuff(phyfile);
	destroyFileBuff(phyfile);
	
	/* output new matrices */
	if(outphyfilename == 0 || (*outphyfilename == '-' && outphyfilename[1] == '-' && outphyfilename[2] == 0)) {
		outphy = stdout;
	} else {
		outphy = sfopen(outphyfilename, "wb");
	}
	printphy(outphy, dist, names, 0, "Merged", format);
	if(numfilename) {
		if(*outnumfilename == '-' && outnumfilename[1] == '-' && outnumfilename[2] == 0) {
			outnum = stdout;
		} else {
			outnum = sfopen(outnumfilename, "wb");
		}
		printphy(outnum, num, names, 0, "Merged", format);
		if(outnum != outphy) {
			fclose(outnum);
		}
	}
	fclose(outphy);
	
	/* clean */
	Matrix_destroy(dist);
	Matrix_destroy(num);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo merges matrices from a multi Phylip file into one matrix\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-n", "Weigh distance with this Phylip file", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-no", "Output number of nucleotides included", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output format", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_merge(int argc, char *argv[]) {
	
	unsigned args, format;
	char *arg, *phyfilename, *numfilename, *outphyfilename, *outnumfilename;
	
	format = 0;
	phyfilename = "--";
	numfilename = 0;
	outphyfilename = "--";
	outnumfilename = "--";
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					phyfilename = argv[args];
				} else {
					missArg("\"-i\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outphyfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "n") == 0) {
				if(++args < argc) {
					numfilename = argv[args];
				} else {
					numfilename = "--";
				}
			} else if(strcmp(arg, "no") == 0) {
				if(++args < argc) {
					outnumfilename = argv[args];
				} else {
					outnumfilename = "--";
				}
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					format = strtoul(argv[args], &arg, 10);
					if(*arg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-f\"");
				}
			} else if(strcmp(arg, "fh") == 0) {
				fprintf(stdout, "Format flags output format, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "# 1:\tRelaxed Phylip\n");
				fprintf(stdout, "#\n");
				return 0;
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
	
	/* merge matrices */
	return merger(phyfilename, numfilename, outphyfilename, outnumfilename, format);
}
