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
#include <stdlib.h>
#include <string.h>
#include "bytescale.h"
#include "cdist.h"
#include "cmdline.h"
#include "dist.h"
#include "fbseek.h"
#include "fsacmp.h"
#include "fsacmpthrd.h"
#include "ltdmatrix.h"
#include "ltdmatrixthrd.h"
#include "matcmp.h"
#include "matrix.h"
#include "meth.h"
#include "methparse.h"
#include "pherror.h"
#include "phy.h"
#include "str.h"
#include "tmp.h"
#include "unionparse.h"

static void makeMatrix(unsigned numFile, char **filenames, char *outputfilename, char *noutputfilename, char *diffilename, char *targetTemplate, double minCov, double alpha, unsigned norm, unsigned minDepth, unsigned minLength, unsigned proxi, unsigned flag, double (*veccmp)(short unsigned*, short unsigned*, int, int), char *methfilename, int tnum) {
	
	int i, informat, cSize, unionin;
	unsigned n, pos, len, **includes;
	long unsigned **seqs;
	unsigned char *include, *trans;
	FILE *outfile, *noutfile, *tmp, *diffile;
	FileBuff *infile, *unionfile;
	Matrix *distMat, *nDest;
	MatrixCounts *mat1;
	MethMotif *motif;
	NucCount *mat2;
	Qseqs *ref, *seq, *header;
	TimeStamp **targetStamps;
	UnionEntry *entry;
	
	/* init */
	distMat = 0;
	nDest = 0;
	targetStamps = 0;
	include = 0;
	motif = 0;
	unionin = 0;
	seqs = 0;
	includes = 0;
	
	/* open output */
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	if(noutputfilename) {
		if(strcmp(noutputfilename, outputfilename) == 0) {
			noutfile = outfile;
		} else if(*noutputfilename == '-' && noutputfilename[1] == 0) {
			noutfile = stdout;
		} else {
			noutfile = sfopen(noutputfilename, "wb");
		}
	} else {
		noutfile = 0;
	}
	if(diffilename) {
		if(strcmp(diffilename, outputfilename) == 0) {
			diffile = outfile;
		} else if(*diffilename == '-' && diffilename[1] == 0) {
			diffile = stdout;
		} else {
			diffile = sfopen(diffilename, "wb");
		}
	} else {
		diffile = 0;
	}
	
	/* determine format of input */
	infile = setFileBuff(1048576);
	if(flag & 16) {
		informat = '>';
	} else if(numFile) {
		informat = fileExist(infile, *filenames);
	} else {
		informat = getc(stdout);
		ungetc(informat, stdout);
	}
	if(informat != '>') {
		informat = '#';
		unionin = 1;
	}
	
	if(informat == '#') {
		mat1 = initMat(1048576, 128);
		mat2 = initNucCount(128);
		trans = 0;
		ref = 0;
		seq = 0;
		header = 0;
	} else if(informat == '>') {
		mat1 = 0;
		mat2 = 0;
		trans = get2BitTable(flag);
		ref = setQseqs(1048576);
		seq = setQseqs(1048576);
		header = setQseqs(64);
		
		/* get methylationsites */
		if(methfilename) {
			openAndDetermine(infile, methfilename);
			motif = getMethMotifs(infile, seq);
			closeFileBuff(infile);
		}
	} else {
		fprintf(stderr, "Format of inputfiles are unsupported.\n");
		exit(1);
	}
	
	if(targetTemplate && 1 < numFile) {
		/* init */
		distMat = ltdMatrix_init(numFile);
		nDest = ltdMatrix_init(numFile);
		if(informat == '>') {
			seqs = smalloc(numFile * sizeof(long unsigned *));
			includes = (flag & 2) ? smalloc(numFile * sizeof(unsigned *)) : smalloc(sizeof(unsigned *));
			cSize = 32768;
			*seqs = smalloc(cSize * sizeof(long unsigned));
			*includes = smalloc(cSize * sizeof(unsigned));
			i = numFile;
			while(--i) {
				seqs[i] = smalloc(cSize * sizeof(long unsigned));
				if(flag & 2) {
					includes[i] = smalloc(cSize * sizeof(unsigned));
				}
			}
		} else {
			seqs = 0;
			includes = 0;
		}
		
		if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
			ERROR();
		}
		include = smalloc(numFile);
		memset(include, 1, numFile);
		
		/* make ltd matrix */
		if(informat == '#') {
			ltdMatrixThrd(distMat, nDest, mat1, targetStamps, include, targetTemplate, filenames, numFile, norm, minDepth, minLength, minCov, veccmp, tnum);
			//ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, targetTemplate, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
		} else if(informat == '>') {
			cSize = ltdFsaMatrix_get(distMat, nDest, numFile, seqs, cSize, infile, targetStamps, include, includes, targetTemplate, filenames, trans, ref, seq, header, norm, minLength, minCov, flag, proxi, motif, diffile, tnum);
		}
		
		/* print ltd matrix */
		if(1 < distMat->n) {
			printphy(outfile, distMat, filenames, include, targetTemplate, flag);
			if(noutputfilename && 1 < nDest->n) {
				printphy(noutfile, nDest, filenames, include, targetTemplate, flag);
			}
		}
	} else if(numFile < 2 && unionin) {
		/* create many matrices */
		/* requires *.union */
		unionfile = setFileBuff(524288);
		if(filenames) {
			openAndDetermine(unionfile, *filenames);
		} else {
			openAndDetermine(unionfile, "-");
		}
		
		/* get samplenames */
		if(!(filenames = UnionEntry_getHeader(unionfile, &numFile))) {
			fprintf(stderr, "Malformed union input.\n");
			exit(1);
		}
		/* init */
		entry = UnionEntry_init(32, numFile);
		distMat = ltdMatrix_init(numFile);
		nDest = ltdMatrix_init(numFile);
		if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
			ERROR();
		}
		include = smalloc(numFile);
		if(informat == '>') {
			seqs = smalloc(numFile * sizeof(long unsigned *));
			includes = (flag & 2) ? smalloc(numFile * sizeof(unsigned *)) : smalloc(sizeof(unsigned *));
			cSize = 32768;
			*seqs = smalloc(cSize * sizeof(long unsigned));
			*includes = smalloc(cSize * sizeof(unsigned));
			i = numFile;
			while(--i) {
				seqs[i] = smalloc(cSize * sizeof(long unsigned));
				if(flag & 2) {
					includes[i] = smalloc(cSize * sizeof(unsigned));
				}
			}
		} else {
			seqs = 0;
			includes = 0;
		}
		
		/* get filenames */
		n = numFile;
		while(n--) {
			/* truncate */
			len = strlen(filenames[n]);
			if((pos = rstrpos(filenames[n], '.', len)) != -1) {
				filenames[n][pos] = 0;
				len = pos;
			}
			
			/* add file suffix */
			if(flag & 16) {
				/* .fsa.gz */
				strcpy(filenames[n] + len, ".fsa.gz");
				len += 7;
			} else {
				/* .mat.gz */
				strcpy(filenames[n] + len, ".mat.gz");
				len += 7;
			}
			
			/* check if file exists */
			if(!(tmp = fopen(filenames[n], "rb"))) {
				len -= 3;
				filenames[n][len] = 0;
			} else {
				fclose(tmp);
			}
		}
		
		/* parse candidates */
		while(UnionEntry_get(unionfile, entry)) {
			/* get included templates */
			memset(include, 0, numFile);
			n = entry->num;
			while(n--) {
				include[entry->filenames[n]] = 1;
			}
			
			/* get matrix */
			if(informat == '#') {
				//ltdMatrixThrd(distMat, nDest, mat1, targetStamps, include, entry->target, filenames, numFile, norm, minDepth, minLength, minCov, veccmp, tnum);
				ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, entry->target, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
			} else if(informat == '>') {
				cSize = ltdFsaMatrix_get(distMat, nDest, numFile, seqs, cSize, infile, targetStamps, include, includes, entry->target, filenames, trans, ref, seq, header, norm, minLength, minCov, flag, proxi, motif, diffile, tnum);
			}
			
			/* print ltd matrix */
			if(1 < distMat->n) {
				printphy(outfile, distMat, filenames, include, entry->target, flag);
				if(noutputfilename) {
					printphy(noutfile, nDest, filenames, include, entry->target, flag);
				}
			}
		}
		closeFileBuff(unionfile);
		destroyFileBuff(unionfile);
		UnionEntry_destroy(entry);
	} else if(numFile < 2) {
		/* msa input */
		if(filenames) {
			openAndDetermine(infile, *filenames);
		} else {
			openAndDetermine(infile, "-");
		}
		cSize = ltdMsaMatrix_get(infile, outfile, noutfile, ref,seq, header, trans, norm, minLength, minCov, flag, proxi, motif, diffile, tnum);
	} else {
		fprintf(stderr, "Invalid argument combination.\n");
		exit(1);
	}
	
	if(diffile && diffile != stdout) {
		fclose(diffile);
	}
	if(informat == '#') {
		destroyMat(mat1);
		destroyNucCount(mat2);
	} else if(informat == '>') {
		if(unionin || 1 < numFile) {
			i = numFile;
			while(--i) {
				free(seqs[i]);
				if(flag & 2) {
					free(includes[i]);
				}
			}
			free(*seqs);
			free(*includes);
			free(seqs);
			free(includes);
		}
		destroyQseqs(ref);
		destroyQseqs(seq);
		destroyQseqs(header);
	}
	destroyMethMotifs(motif);
	Matrix_destroy(distMat);
	Matrix_destroy(nDest);
	destroyFileBuff(infile);
	free(targetStamps);
	free(include);
	
	/* close */
	fclose(outfile);
	if(noutputfilename && strcmp(noutputfilename, outputfilename) != 0) {
		fclose(noutfile);
	}
}

static int add2Matrix(char *path, char *addfilename, char *outputfilename, char *noutputfilename, char *diffilename, char *targetTemplate, double minCov, unsigned norm, unsigned minDepth, unsigned minLength, unsigned proxi, unsigned flag, char sep, double (*veccmp)(short unsigned*, short unsigned*, int, int), int tnum) {
	
	int n, pos, informat;
	double *D, *N;
	char *ptr;
	FILE *outfile, *noutfile;
	FileBuff *infile;
	Qseqs **names;
	
	/* determine input format */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, outputfilename);
	
	/* convert path to dir */
	pos = -1;
	n = 0;
	ptr = path - 1;
	while(*++ptr) {
		++n;
		if(*ptr == '/') {
			pos = n;
		}
	}
	if(0 <= pos) {
		path[pos] = 0;
	}
	
	/* get n */
	n = getSizePhy(infile);
	D = smalloc(n * sizeof(double));
	N = smalloc(n * sizeof(double));
	
	/* get names */
	if(!(names = getFilenamesPhy(path, n, infile, sep))) {
		ERROR();
	} else if(infile->bytes) {
		fprintf(stderr, "Cannot update a multi distance phylip file.\n");
		return 1;
	}
	
	/* close phyfile */
	closeFileBuff(infile);
	
	/* add final row */
	informat = fileExist(infile, addfilename);
	destroyFileBuff(infile);
	
	if(informat == '>') {
		if(ltdFsaRowThrd(D, N, targetTemplate, addfilename, diffilename, names, n, norm, minLength, minCov, flag, proxi, tnum)) {
			fprintf(stderr, "Distance measures failed and thus the matrix was not updated.\n");
			return 1;
		}
	} else {
		/* calculate new row */
		if(ltdRowThrd(D, N, targetTemplate, addfilename, names, n, norm, minDepth, minLength, minCov, veccmp, tnum)) {
			fprintf(stderr, "Distance measures failed and thus the matrix was not updated.\n");
			return 1;
		}
	}
	
	/* open output and init new row(s) */
	outfile = sfopen(outputfilename, "rb+");
	printphyUpdate(outfile, ++n, addfilename, D, flag);
	fclose(outfile);
	if(noutputfilename) {
		noutfile = sfopen(noutputfilename, "rb+");
		printphyUpdate(noutfile, n, addfilename, N, flag);
		fclose(noutfile);
	}
	
	/* clean */
	pos = n - 1;
	while(pos--) {
		destroyQseqs(names[pos]);
	}
	free(names);
	free(D);
	free(N);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo dist calculates distances between samples based on overlaps between nucleotide count matrices created by e.g. KMA.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file(s)", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'n', "nucleotide_numbers", "Output number of nucleotides included", "False/None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'x', "print_precision", "Floating point print precision", "9");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'y', "methylation_motifs", "Mask methylation motifs from <file>", "False/None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'V', "nucleotide_variations", "Output nucleotide variations", "False/None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'r', "reference", "Target reference", "None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'a', "add", "Add file to existing matrix", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'E', "min_depth", "Minimum depth", "15");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'C', "min_cov", "Minimum coverage", "50.0%");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'L', "min_len", "Minimum overlapping length", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'W', "normalization_weight", "Normalization weight", "0 / None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'P', "proximity", "Minimum proximity between SNPs", "0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'd', "distance", "Distance method", "cos");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'D', "distance_help", "Help on option \"-d\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'l', "significance_lvl", "Minimum lvl. of signifiacnce", "0.05");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'p', "float_precision", "Float precision on distance matrix", "double");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 's', "short_precision", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "byte_precision", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'H', "mmap", "Allocate matrix on the disk", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tmp", "Set directory for temporary files", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 't', "threads", "Number of threads", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	n	n	nucleotide_numbers
	m	y	methylation_motifs
	nv	V	nucleotide_variations
	r	r	reference
	a	a	add
	md	E	min_depth
	mc	C	min_cov
	ml	L	min_len
	nm	W	normalization_weight
	pr	P	proximity
	f	f	flag
	fh	F	flag_help
	d	d	distance
	dh	D	distance_help
	p	l	significance_lvl
	fp	p	float_precision
	sp	s	short_precision
	bp	b	byte_precision
	mm	H	mmap
	tmp	T	tmp
	t	t	threads
	h	h	help
	*/
}

int main_dist(int argc, char **argv) {
	
	const char *stdstream = "-";
	int size, len, offset, args, precision;
	unsigned numFile, flag, norm, minDepth, minLength, proxi, n, t;
	char **Arg, *arg, *targetTemplate, **filenames, *addfilename, *errorMsg;
	char *outputfilename, *noutputfilename, *methfilename, *diffilename;
	char *method, *tmp, opt, sep;
	double minCov, alpha;
	double (*veccmp)(short unsigned*, short unsigned*, int, int);
	
	/* set defaults */
	precision = 9;
	size = sizeof(double);
	numFile = 0;
	flag = 1;
	norm = 0;
	minDepth = 15;
	minLength = 1;
	proxi = 0;
	t = 1;
	targetTemplate = 0;
	filenames = 0;
	addfilename = 0;
	outputfilename = (char *)(stdstream);
	noutputfilename = 0;
	methfilename = 0;
	diffilename = 0;
	minCov = 0.5;
	alpha = 0.05;
	method = "cos";
	veccmp = &coscmp;
	tmp = 0;
	sep = '\t';
	
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
					filenames = getArgListDie(&Arg, &args, len + offset, "input");
					numFile = getArgListLen(&Arg, &args);
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "nucleotide_numbers") == 0) {
					noutputfilename = getArgDie(&Arg, &args, len + offset, "nucleotide_numbers");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "print_precision") == 0) {
					precision = getNumArg(&Arg, &args, len + offset, "print_precision");
				} else if(cmdcmp(arg, "methylation_motifs") == 0) {
					methfilename = getArgDie(&Arg, &args, len + offset, "methylation_motifs");
				} else if(cmdcmp(arg, "nucleotide_variations") == 0) {
					diffilename = getArgDie(&Arg, &args, len + offset, "nucleotide_variations");
				} else if(cmdcmp(arg, "reference") == 0) {
					targetTemplate = getArgDie(&Arg, &args, len + offset, "reference");
				} else if(cmdcmp(arg, "add") == 0) {
					addfilename = getArgDie(&Arg, &args, len + offset, "add");
				} else if(cmdcmp(arg, "min_depth") == 0) {
					minDepth = getdArg(&Arg, &args, len + offset, "min_depth");
				} else if(cmdcmp(arg, "min_cov") == 0) {
					minCov = getdArg(&Arg, &args, len + offset, "min_cov") / 100;
				} else if(cmdcmp(arg, "min_len") == 0) {
					minLength = getNumArg(&Arg, &args, len + offset, "min_len");
				} else if(cmdcmp(arg, "normalization_weight") == 0) {
					norm = getNumArg(&Arg, &args, len + offset, "normalization_weight");
				} else if(cmdcmp(arg, "proximity") == 0) {
					proxi = getNumArg(&Arg, &args, len + offset, "proximity");
				} else if(cmdcmp(arg, "flag") == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(cmdcmp(arg, "flag_help") == 0) {
					flag = -1;
				} else if(cmdcmp(arg, "distance") == 0) {
					method = getArgDie(&Arg, &args, len + offset, "distance");
				} else if(cmdcmp(arg, "distance_help") == 0) {
					method = 0;
				} else if(cmdcmp(arg, "significance_lvl") == 0) {
					alpha = getdArg(&Arg, &args, len + offset, "significance_lvl");
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
				} else if(cmdcmp(arg, "threads") == 0) {
					t = getNumArg(&Arg, &args, len + offset, "threads");
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
						filenames = getArgListDie(&Arg, &args, len, "i");
						numFile = getArgListLen(&Arg, &args);
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'n') {
						noutputfilename = getArgDie(&Arg, &args, len, "n");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'x') {
						precision = getNumArg(&Arg, &args, len, "x");
						opt = 0;
					} else if(opt == 'y') {
						methfilename = getArgDie(&Arg, &args, len, "y");
						opt = 0;
					} else if(opt == 'V') {
						diffilename = getArgDie(&Arg, &args, len, "V");
						opt = 0;
					} else if(opt == 'r') {
						targetTemplate = getArgDie(&Arg, &args, len, "r");
						opt = 0;
					} else if(opt == 'a') {
						addfilename = getArgDie(&Arg, &args, len, "a");
						opt = 0;
					} else if(opt == 'E') {
						minDepth = getdArg(&Arg, &args, len, "E");
						opt = 0;
					} else if(opt == 'C') {
						minCov = getdArg(&Arg, &args, len, "C") / 100;
						opt = 0;
					} else if(opt == 'L') {
						minLength = getNumArg(&Arg, &args, len, "L");
						opt = 0;
					} else if(opt == 'W') {
						norm = getNumArg(&Arg, &args, len, "W");
						opt = 0;
					} else if(opt == 'P') {
						proxi = getNumArg(&Arg, &args, len, "P");
						opt = 0;
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
					} else if(opt == 'd') {
						method = getArgDie(&Arg, &args, len, "d");
						opt = 0;
					} else if(opt == 'D') {
						method = 0;
					} else if(opt == 'l') {
						alpha = getdArg(&Arg, &args, len, "l");
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
					} else if(opt == 't') {
						t = getNumArg(&Arg, &args, len, "t");
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
		filenames = Arg;
		numFile = args;
	}
	
	/* verify input */
	if(minCov < 0 || 1 < minCov) {
		invaArg("\"--min_cov\"");
	}
	if(ByteScale == 0) {
		if(size == sizeof(short unsigned)) {
			invaArg("\"--short_precision\"");
		} else {
			invaArg("\"--byte_precision\"");
		}
	}
	if(alpha < 0) {
		invaArg("\"--significance_lvl\"");
	}
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "# Format flags output, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "#   1:\tRelaxed Phylip\n");
		fprintf(stdout, "#   2:\tDistances are pairwise, always true on *.mat files\n");
		fprintf(stdout, "#   4:\tInclude template name in phylip file\n");
		fprintf(stdout, "#   8:\tInclude insignificant bases in distance calculation, only affects fasta input\n");
		fprintf(stdout, "#  16:\tDistances based on fasta input\n");
		fprintf(stdout, "#  32:\tDo not include insignificant bases in pruning\n");
		fprintf(stdout, "#\n");
		return 0;
	}
	
	/* distance method */
	if(method == 0) {
		fprintf(stdout, "# Distance calculation methods:\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "# cos:\tCalculate distance between positions as the angle between the count vectors.\n");
		fprintf(stdout, "# z:\tMake consensus comparison if vectors passes a McNemar test\n");
		fprintf(stdout, "# chi2:\tCalculate the chi square distance\n");
		fprintf(stdout, "# nchi2:\tCalculate the normalized chi square distance\n");
		fprintf(stdout, "# c:\tCalculate the Clausen distance between the count vectors. d(A,B) = (||A-B||_1 / sum(max{Ai, Bi}))\n");
		fprintf(stdout, "# nc:\tCalculate the normalized Clausen distance between the count vectors.\n");
		fprintf(stdout, "# bc:\tCalculate the Bray-Curtis dissimilarity between the count vectors.\n");
		fprintf(stdout, "# nbc:\tCalculate the normalized Bray-Curtis dissimilarity between the count vectors.\n");
		fprintf(stdout, "# ln:\tCalculate distance between positions as the n-norm distance between the count vectors. Replace \"n\" with the waned norm\n");
		fprintf(stdout, "# linf:\tCalculate distance between positions as the l_infinity distance between the count vectors.\n");
		fprintf(stdout, "# nln:\tCalculate distance between positions as the normalized n-norm distance between the count vectors. Replace last \"n\" with the waned norm\n");
		fprintf(stdout, "# nlinf:\tCalculate distance between positions as the normalized l_infinity distance between the count vectors.\n");
		fprintf(stdout, "#\n");
		return 0;
	} else if(strcmp(method, "cos") == 0) {
		veccmp = &coscmp;
	} else if(strcmp(method, "z") == 0) {
		veccmp = &zcmp;
	} else if(strcmp(method, "chi2") == 0) {
		veccmp = &chi2cmp;
	} else if(strcmp(method, "nchi2") == 0) {
		veccmp = &nchi2cmp;
	} else if(strcmp(method, "nc") == 0) {
		veccmp = &nccmp;
	} else if(strcmp(method, "c") == 0) {
		veccmp = &ccmp;
	} else if(strcmp(method, "np") == 0) {
		veccmp = &npcmp;
	} else if(strcmp(method, "p") == 0) {
		veccmp = &pcmp;
	} else if(strcmp(method, "nbc") == 0) {
		veccmp = &nbccmp;
	} else if(strcmp(method, "bc") == 0) {
		veccmp = &bccmp;
	} else if(strcmp(method, "nl1") == 0) {
		veccmp = &nl1cmp;
	} else if(strcmp(method, "nl2") == 0) {
		veccmp = &nl2cmp;
	} else if(strcmp(method, "nlinf") == 0) {
		veccmp = &nlinfcmp;
	} else if(strcmp(method, "l1") == 0) {
		veccmp = &l1cmp;
	} else if(strcmp(method, "l2") == 0) {
		veccmp = &l2cmp;
	} else if(strcmp(method, "linf") == 0) {
		veccmp = &linfcmp;
	} else if(*method == 'l') {
		veccmp = &lncmp;
		n = strtoul(method + 1, &errorMsg, 10);
		if(*errorMsg != 0) {
			invaArg("\"-d ln\"");
		}
		veccmp(0, (short unsigned *)(&n), 0, 0);
	} else if(strncmp(method, "nl", 2) == 0) {
		veccmp = &nlncmp;
		n = strtoul(method + 2, &errorMsg, 10);
		if(*errorMsg != 0) {
			invaArg("\"-d nln\"");
		}
		veccmp(0, (short unsigned *)(&n), 0, 0);
	} else {
		invaArg("\"-d\"");
	}
	
	/* set print precision */
	setPrecisionPhy(precision);
	
	/* tmp dir */
	if(tmp) {
		tmpF(tmp);
	}
	
	/* set precision */
	ltdMatrixInit(-size);
	ltdMatrixMinit(-size);
	
	/* set function variables */
	zcmp(0, (short unsigned *)(&alpha), 0, 0);
	if(flag & 32) {
		getIncPosPtr = &getIncPosInsigPrune;
	} else if(flag & 8) {
		getIncPosPtr = &getIncPosInsig;
	}
	
	/* check for required input */
	if(!numFile && targetTemplate) {
		numFile = 1;
	}
	
	if(addfilename && filenames) {
		return add2Matrix(*filenames, addfilename, outputfilename, noutputfilename, diffilename, targetTemplate, minCov, norm, minDepth, minLength, proxi, flag, sep, veccmp, t);
	} else {
		makeMatrix(numFile, filenames, outputfilename, noutputfilename, diffilename, targetTemplate, minCov, alpha, norm, minDepth, minLength, proxi, flag, veccmp, methfilename, t);
	}
	
	return 0;
}
