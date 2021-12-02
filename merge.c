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
#include "bytescale.h"
#include "filebuff.h"
#include "hashmapstr.h"
#include "hashmapstrindex.h"
#include "matrix.h"
#include "merge.h"
#include "phy.h"
#include "qseqs.h"
#include "ulist.h"
#include "tmp.h"
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
	float *Dfptr, *Nfptr;
	short unsigned *Dsptr, *Nsptr;
	unsigned char *Dbptr, *Nbptr;
	
	/* normalize new distance matrix */
	i = (D->n * (D->n - 1)) / 2 + 1;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
		Nptr = *(N->mat) - 1;
		while(--i) {
			if(*++Nptr != 0) {
				*++Dptr /= *Nptr;
			} else {
				*++Dptr = -1.0;
			}
		}
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
		Nfptr = *(N->fmat) - 1;
		while(--i) {
			if(*++Nfptr != 0) {
				*++Dfptr /= *Nfptr;
			} else {
				*++Dfptr = -1.0;
			}
		}
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
		Nsptr = *(N->smat) - 1;
		while(--i) {
			if(*++Nsptr != 0) {
				++Dsptr;
				*Dsptr = dtouc(((uctod(*Dsptr)) / (uctod(*Nsptr))));
			} else {
				*++Dsptr = uctod(-1.0);
			}
		}
	} else {
		Dbptr = *(D->bmat) - 1;
		Nbptr = *(N->bmat) - 1;
		while(--i) {
			if(*++Nbptr != 0) {
				++Dbptr;
				*Dbptr = dtouc(((uctod(*Dbptr)) / (uctod(*Nbptr))));
			} else {
				*++Dbptr = uctod(-1.0);
			}
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
	
	unsigned i, j, m, n, *distIndex, *distIndexM, *distIndexN;
	double *Dptr, *Nptr, **distMat, **numMat;
	float *Dfptr, *Nfptr, **distfMat, **numfMat;
	short unsigned *Dsptr, *Nsptr, **distsMat, **numsMat;
	unsigned char *Dbptr, *Nbptr, **distbMat, **numbMat;
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
	if(dist->mat) {
		Dptr = *(dist->mat) - 1;
		Nptr = *(num->mat) - 1;
		while(--i) {
			*++Dptr *= *++Nptr;
		}
	} else if(dist->fmat) {
		Dfptr = *(dist->fmat) - 1;
		Nfptr = *(num->fmat) - 1;
		while(--i) {
			*++Dfptr *= *++Nfptr;
		}
	} else if(dist->smat) {
		Dsptr = *(dist->smat) - 1;
		Nsptr = *(num->smat) - 1;
		while(--i) {
			*++Dsptr *= *++Nsptr;
		}
	} else {
		Dbptr = *(dist->bmat) - 1;
		Nbptr = *(num->bmat) - 1;
		while(--i) {
			*++Dbptr *= *++Nbptr;
		}
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
		if(dist->mat) {
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
					if(m < n) {
						distMat[n][m] += *++Dptr * *++Nptr;
						numMat[n][m] += *Nptr;
					} else {
						distMat[m][n] += *++Dptr * *++Nptr;
						numMat[m][n] += *Nptr;
					}
				}
			}
		} else if(dist->fmat) {
			distfMat = dist->fmat;
			numfMat = num->fmat;
			Dfptr = *(D->fmat) - 1;
			Nfptr = *(N->fmat) - 1;
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
					if(m < n) {
						distfMat[n][m] += *++Dfptr * *++Nfptr;
						numfMat[n][m] += *Nfptr;
					} else {
						distfMat[m][n] += *++Dfptr * *++Nfptr;
						numfMat[m][n] += *Nfptr;
					}
				}
			}
		} else if(dist->smat) {
			distsMat = dist->smat;
			numsMat = num->smat;
			Dsptr = *(D->smat) - 1;
			Nsptr = *(N->smat) - 1;
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
					if(m < n) {
						distsMat[n][m] += dtouc((uctod(*++Dsptr) * uctod(*++Nsptr)));
						numsMat[n][m] += *Nsptr;
					} else {
						distsMat[m][n] += dtouc((uctod(*++Dsptr) * uctod(*++Nsptr)));
						numsMat[m][n] += *Nsptr;
					}
				}
			}
		} else {
			distbMat = dist->bmat;
			numbMat = num->bmat;
			Dbptr = *(D->bmat) - 1;
			Nbptr = *(N->bmat) - 1;
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
					if(m < n) {
						distbMat[n][m] += dtouc((uctod(*++Dbptr) * uctod(*++Nbptr)));
						numbMat[n][m] += *Nbptr;
					} else {
						distbMat[m][n] += dtouc((uctod(*++Dbptr) * uctod(*++Nbptr)));
						numbMat[m][n] += *Nbptr;
					}
				}
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
	Matrix_destroy(N);
	uList_destroy(dist_index);
	
	/* return names w.r.t. merged matrix */
	return getOrderedNames(names_index);
}

char ** jl_merge(Matrix *dist, Matrix *num, FileBuff *phyfile) {
	
	unsigned i, j, m, n, *distIndex, *distIndexM, *distIndexN;
	double *Dptr, **distMat, **numMat;
	float *Dfptr, **distfMat, **numfMat;
	short unsigned *Dsptr, **distsMat, **numsMat;
	unsigned char *Dbptr, **distbMat, **numbMat;
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
	if(num->mat) {
		Dptr = *(num->mat) - 1;
		while(--i) {
			*++Dptr = 1;
		}
	} else if(num->fmat) {
		Dfptr = *(num->fmat) - 1;
		while(--i) {
			*++Dfptr = 1;
		}
	} else if(num->smat) {
		Dsptr = *(num->smat) - 1;
		while(--i) {
			*++Dsptr = 1;
		}
	} else {
		Dbptr = *(num->bmat) - 1;
		while(--i) {
			*++Dbptr = 1;
		}
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
		if(dist->mat) {
			distMat = dist->mat;
			numMat = num->mat;
			Dptr = *(D->mat) - 1;
			distIndexM = dist_index->list;
			distIndex = distIndexM - 1;
			i = 1;
			while(i != D->n) {
				m = *++distIndexM;
				distIndexN = distIndex;
				j = ++i;
				while(--j) {
					n = *++distIndexN;
					/* update cell */
					if(m < n) {
						distMat[n][m] += *++Dptr;
						++numMat[n][m];
					} else {
						distMat[m][n] += *++Dptr;
						++numMat[m][n];
					}
				}
			}
		} else if(dist->fmat) {
			distfMat = dist->fmat;
			numfMat = num->fmat;
			Dfptr = *(D->fmat) - 1;
			distIndexM = dist_index->list;
			distIndex = distIndexM - 1;
			i = 1;
			while(i != D->n) {
				m = *++distIndexM;
				distIndexN = distIndex;
				j = ++i;
				while(--j) {
					n = *++distIndexN;
					/* update cell */
					if(m < n) {
						distfMat[n][m] += *++Dfptr;
						++numfMat[n][m];
					} else {
						distfMat[m][n] += *++Dfptr;
						++numfMat[m][n];
					}
				}
			}
		} else if(dist->smat) {
			distsMat = dist->smat;
			numsMat = num->smat;
			Dsptr = *(D->smat) - 1;
			distIndexM = dist_index->list;
			distIndex = distIndexM - 1;
			i = 1;
			while(i != D->n) {
				m = *++distIndexM;
				distIndexN = distIndex;
				j = ++i;
				while(--j) {
					n = *++distIndexN;
					/* update cell */
					if(m < n) {
						distsMat[n][m] += *++Dsptr;
						++numsMat[n][m];
					} else {
						distsMat[m][n] += *++Dsptr;
						++numsMat[m][n];
					}
				}
			}
		} else {
			distbMat = dist->bmat;
			numbMat = num->bmat;
			Dbptr = *(D->bmat) - 1;
			distIndexM = dist_index->list;
			distIndex = distIndexM - 1;
			i = 1;
			while(i != D->n) {
				m = *++distIndexM;
				distIndexN = distIndex;
				j = ++i;
				while(--j) {
					n = *++distIndexN;
					/* update cell */
					if(m < n) {
						distbMat[n][m] += *++Dbptr;
						++numbMat[n][m];
					} else {
						distbMat[m][n] += *++Dbptr;
						++numbMat[m][n];
					}
				}
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input multi phylip distance file", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-n", "Weigh distance with this Phylip file", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-no", "Output number of nucleotides included", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output format", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	//fprintf(out, "# %16s\t%-32s\t%s\n", "-sp", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_merge(int argc, char *argv[]) {
	
	int size;
	unsigned args, format;
	char *arg, *phyfilename, *numfilename, *outphyfilename, *outnumfilename;
	char *errorMsg;
	
	/* set defaults */
	size = sizeof(double);
	format = 1;
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
				fprintf(stdout, "# Format flags output, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#   1:\tRelaxed Phylip\n");
				fprintf(stdout, "#   2:\tDistances are pairwise, always true on *.mat files\n");
				fprintf(stdout, "#   4:\tInclude template name in phylip file\n");
				fprintf(stdout, "#   8:\tInclude insignificant bases in distance calculation, only affects fasta input\n");
				fprintf(stdout, "#  16:\tInclude reference, only affects fasta input\n");
				fprintf(stdout, "#  32:\tDistances based on fasta input\n");
				fprintf(stdout, "#\n");
				return 0;
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
	
	/* merge matrices */
	return merger(phyfilename, numfilename, outphyfilename, outnumfilename, format);
}
