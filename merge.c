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
#include "cmdline.h"
#include "filebuff.h"
#include "hashmapstr.h"
#include "hashmapstrindex.h"
#include "matrix.h"
#include "merge.h"
#include "phy.h"
#include "qseqs.h"
#include "ulist.h"
#include "tmp.h"

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
				*Dsptr = dtouc(((uctod(*Dsptr)) / (uctod(*Nsptr))), 0.5);
			} else {
				*++Dsptr = dtouc(-1.0, 0);
			}
		}
	} else {
		Dbptr = *(D->bmat) - 1;
		Nbptr = *(N->bmat) - 1;
		while(--i) {
			if(*++Nbptr != 0) {
				++Dbptr;
				*Dbptr = dtouc(((uctod(*Dbptr)) / (uctod(*Nbptr))), 0.5);
			} else {
				*++Dbptr = dtouc(-1.0, 0);
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

char ** merge(Matrix *dist, Matrix *num, FileBuff *phyfile, FileBuff *numfile, char sep, char quotes) {
	
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
	names = loadPhy(dist, 0, 0, phyfile, sep, quotes);
	loadPhy(num, names, 0, numfile, sep, quotes);
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
	while((names = loadPhy(D, names, 0, phyfile, sep, quotes)) && D->n) {
		if(!loadPhy(N, names, 0, numfile, sep, quotes) || N->n != D->n) {
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
						distsMat[n][m] += dtouc((uctod(*++Dsptr) * uctod(*++Nsptr)), 0.5);
						numsMat[n][m] += *Nsptr;
					} else {
						distsMat[m][n] += dtouc((uctod(*++Dsptr) * uctod(*++Nsptr)), 0.5);
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
						distbMat[n][m] += dtouc((uctod(*++Dbptr) * uctod(*++Nbptr)), 0.5);
						numbMat[n][m] += *Nbptr;
					} else {
						distbMat[m][n] += dtouc((uctod(*++Dbptr) * uctod(*++Nbptr)), 0.5);
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

char ** jl_merge(Matrix *dist, Matrix *num, FileBuff *phyfile, char sep, char quotes) {
	
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
	names = loadPhy(dist, 0, 0, phyfile, sep, quotes);
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
	while((names = loadPhy(D, names, 0, phyfile, sep, quotes)) && D->n) {
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

int merger(char *phyfilename, char *numfilename, char *outphyfilename, char *outnumfilename, unsigned format, char sep, char quotes) {
	
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
		names = merge(dist, num, phyfile, numfile, sep, quotes);
		closeFileBuff(numfile);
		destroyFileBuff(numfile);
	} else {
		names = jl_merge(dist, num, phyfile, sep, quotes);
	}
	closeFileBuff(phyfile);
	destroyFileBuff(phyfile);
	
	/* output new matrices */
	if(outphyfilename == 0 || (*outphyfilename == '-' && outphyfilename[1] == 0)) {
		outphy = stdout;
	} else {
		outphy = sfopen(outphyfilename, "wb");
	}
	printphy(outphy, dist, names, 0, "Merged", format);
	if(numfilename) {
		if(*outnumfilename == '-' && outnumfilename[1] == 0) {
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
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input multi phylip distance file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'w', "nucleotides_weights", "Weigh distance with this Phylip file", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'n', "nucleotide_numbers", "Output number of nucleotides included", "False/None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'x', "print_precision", "Floating point print precision", "9");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
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
	n	w	nucleotides_weights
	no	n	nucleotide_numbers
	f	f	flag
	fh	F	flag_help
	fp	p	float_precision
	sp	s	short_precision
	bp	b	byte_precision
	mm	H	mmap
	tmp	T	tmp
	h	h	help
	*/
}

int main_merge(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, flag, size, len, offset, precision;
	char **Arg, *arg, *phyfilename, *numfilename, *outphyfilename;
	char *outnumfilename, *tmp, opt, sep, quotes;
	
	/* set defaults */
	size = sizeof(double);
	precision = 9;
	flag = 1;
	phyfilename = (char *)(stdstream);
	numfilename = 0;
	outphyfilename = (char *)(stdstream);
	outnumfilename = (char *)(stdstream);
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
					phyfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(cmdcmp(arg, "output") == 0) {
					outphyfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "print_precision") == 0) {
					precision = getNumArg(&Arg, &args, len + offset, "print_precision");
				} else if(cmdcmp(arg, "nucleotides_weights") == 0) {
					numfilename = getArgDie(&Arg, &args, len + offset, "nucleotides_weights");
				} else if(cmdcmp(arg, "nucleotide_numbers") == 0) {
					outnumfilename = getArgDie(&Arg, &args, len + offset, "nucleotide_numbers");
				} else if(cmdcmp(arg, "flag") == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(cmdcmp(arg, "flag_help") == 0) {
					flag = -1;
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
						phyfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outphyfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'x') {
						precision = getNumArg(&Arg, &args, len, "x");
						opt = 0;
					} else if(opt == 'w') {
						numfilename = getArgDie(&Arg, &args, len, "w");
						opt = 0;
					} else if(opt == 'n') {
						outnumfilename = getArgDie(&Arg, &args, len, "n");
						opt = 0;
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
					} else if(opt == 'p') {
						size = sizeof(float);
					} else if(opt == 's') {
						size = sizeof(short unsigned);
						ByteScale = getdDefArg(&Arg, &args, len, ByteScale, "s");
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
		phyfilename = *Arg;
		if(--args) {
			nonOptError();
		}
	}
	
	/* set print precision */
	setPrecisionPhy(precision);
	
	/* tmp dir */
	if(tmp) {
		tmpF(tmp);
	}
	
	if(flag == -1) {
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
	}
	
	/* set precision */
	ltdMatrixInit(-size);
	ltdMatrixMinit(-size);
	
	/* merge matrices */
	return merger(phyfilename, numfilename, outphyfilename, outnumfilename, flag, sep, quotes);
}
