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
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static void makeMatrix(unsigned numFile, char **filenames, char *outputfilename, char *noutputfilename, char *diffilename, char *targetTemplate, double minCov, double alpha, unsigned norm, unsigned minDepth, unsigned minLength, unsigned proxi, unsigned flag, double (*veccmp)(short unsigned*, short unsigned*, int, int), char *methfilename, int tnum) {
	
	int i, informat, cSize;
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
	
	/* open output */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	if(noutputfilename) {
		if(strcmp(noutputfilename, outputfilename) == 0) {
			noutfile = outfile;
		} else if(*noutputfilename == '-' && noutputfilename[1] == '-' && noutputfilename[2] == 0) {
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
		} else if(*diffilename == '-' && diffilename[1] == '-' && diffilename[2] == 0) {
			diffile = stdout;
		} else {
			diffile = sfopen(diffilename, "wb");
		}
	} else {
		diffile = 0;
	}
	motif = 0;
	
	/* determine format of input */
	infile = setFileBuff(1048576);
	if(flag & 16) {
		informat = '>';
	} else if(numFile && numFile != 1) {
		informat = fileExist(infile, *filenames);
	} else {
		informat = '#';
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
	} else if(numFile < 2) {
		/* create many matrices */
		/* requires *.union */
		unionfile = setFileBuff(524288);
		if(filenames) {
			openAndDetermine(unionfile, *filenames);
		} else {
			openAndDetermine(unionfile, "--");
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

static int add2Matrix(char *path, char *addfilename, char *outputfilename, char *noutputfilename, char *diffilename, char *targetTemplate, double minCov, unsigned norm, unsigned minDepth, unsigned minLength, unsigned proxi, unsigned flag, double (*veccmp)(short unsigned*, short unsigned*, int, int), int tnum) {
	
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
	if(!(names = getFilenamesPhy(path, n, infile))) {
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-n", "Output number of nucleotides included", "False/None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-m", "Mask methylation motifs from <file>", "False/None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nv", "Output nucleotide variations", "False/None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-r", "Target reference", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-a", "Add file to existing matrix", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-md", "Minimum depth", "15");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mc", "Minimum coverage", "50.0%");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum overlapping length", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nm", "Normalization", "1000000");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-pr", "Minimum proximity between SNPs", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "cos");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-p", "Minimum lvl. of signifiacnce", "0.05");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Number of threads", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_dist(int argc, char *argv[]) {
	
	int size;
	unsigned args, numFile, flag, norm, minDepth, minLength, proxi, n, t;
	char *arg, *targetTemplate, **filenames, *addfilename, *errorMsg;
	char *outputfilename, *noutputfilename, *methfilename, *diffilename;
	double minCov, alpha;
	double (*veccmp)(short unsigned*, short unsigned*, int, int);
	
	/* set defaults */
	size = sizeof(double);
	numFile = 0;
	flag = 1;
	norm = 1000000;
	minDepth = 15;
	minLength = 1;
	proxi = 0;
	t = 1;
	targetTemplate = 0;
	filenames = 0;
	addfilename = 0;
	outputfilename = "--";
	noutputfilename = 0;
	methfilename = 0;
	diffilename = 0;
	minCov = 0.5;
	alpha = 0.05;
	veccmp = &coscmp;
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					filenames = argv + args;
					numFile = 1;
					while(++args < argc && (*argv[args] != '-' || (argv[args][1] == '-' && argv[args][2] == 0))) {
						++numFile;
					}
					if(numFile == 0) {
						missArg("\"-i\"");
					}
					--args;
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "n") == 0) {
				if(++args < argc) {
					noutputfilename = argv[args];
				} else {
					noutputfilename = "--";
				}
			} else if(strcmp(arg, "m") == 0) {
				if(++args < argc) {
					methfilename = argv[args];
				} else {
					missArg("\"-m\"");
				}
			} else if(strcmp(arg, "nv") == 0) {
				if(++args < argc) {
					diffilename = argv[args];
				} else {
					diffilename = "--";
				}
			} else if(strcmp(arg, "r") == 0) {
				if(++args < argc) {
					targetTemplate = argv[args];
				} else {
					missArg("\"-r\"");
				}
			} else if(strcmp(arg, "a") == 0) {
				if(++args < argc) {
					addfilename = argv[args];
				} else {
					missArg("\"-a\"");
				}
			} else if(strcmp(arg, "md") == 0) {
				if(++args < argc) {
					if((minDepth = strtoul(argv[args], &errorMsg, 10)) == 0) {
						minDepth = 1;
					}
					if(*errorMsg != 0) {
						invaArg("\"-md\"");
					}
				} else {
					missArg("\"-md\"");
				}
			} else if(strcmp(arg, "ml") == 0) {
				if(++args < argc) {
					if((minLength = strtoul(argv[args], &errorMsg, 10)) == 0) {
						minLength = 1;
					}
					if(*errorMsg != 0) {
						invaArg("\"-ml\"");
					}
				} else {
					missArg("\"-ml\"");
				}
			} else if(strcmp(arg, "mc") == 0) {
				if(++args < argc) {
					minCov = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || minCov < 0 || 100 < minCov) {
						invaArg("\"-mc\"");
					}
					minCov /= 100.0;
				} else {
					missArg("\"-mc\"");
				}
			} else if(strcmp(arg, "nm") == 0) {
				if(++args < argc) {
					norm = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-nm\"");
					}
				} else {
					missArg("\"-nm\"");
				}
			} else if(strcmp(arg, "pr") == 0) {
				if(++args < argc) {
					proxi = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-pr\"");
					}
				} else {
					missArg("\"-pr\"");
				}
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					flag = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
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
				fprintf(stdout, "#  16:\tDistances based on fasta input\n");
				fprintf(stdout, "#  32:\tDo not include insignificant bases in pruning\n");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "cos") == 0) {
						veccmp = &coscmp;
					} else if(strcmp(arg, "z") == 0) {
						veccmp = &zcmp;
					} else if(strcmp(arg, "chi2") == 0) {
						veccmp = &chi2cmp;
					} else if(strcmp(arg, "nchi2") == 0) {
						veccmp = &nchi2cmp;
					} else if(strcmp(arg, "nc") == 0) {
						veccmp = &nccmp;
					} else if(strcmp(arg, "c") == 0) {
						veccmp = &ccmp;
					} else if(strcmp(arg, "np") == 0) {
						veccmp = &npcmp;
					} else if(strcmp(arg, "p") == 0) {
						veccmp = &pcmp;
					} else if(strcmp(arg, "nbc") == 0) {
						veccmp = &nbccmp;
					} else if(strcmp(arg, "bc") == 0) {
						veccmp = &bccmp;
					} else if(strcmp(arg, "nl1") == 0) {
						veccmp = &nl1cmp;
					} else if(strcmp(arg, "nl2") == 0) {
						veccmp = &nl2cmp;
					} else if(strcmp(arg, "nlinf") == 0) {
						veccmp = &nlinfcmp;
					} else if(strcmp(arg, "l1") == 0) {
						veccmp = &l1cmp;
					} else if(strcmp(arg, "l2") == 0) {
						veccmp = &l2cmp;
					} else if(strcmp(arg, "linf") == 0) {
						veccmp = &linfcmp;
					} else if(*arg == 'l') {
						veccmp = &lncmp;
						n = strtoul(arg + 1, &errorMsg, 10);
						if(*errorMsg != 0) {
							invaArg("\"-d ln\"");
						}
						veccmp(0, (short unsigned *)(&n), 0, 0);
					} else if(strncmp(arg, "nl", 2) == 0) {
						veccmp = &nlncmp;
						n = strtoul(arg + 2, &errorMsg, 10);
						if(*errorMsg != 0) {
							invaArg("\"-d nln\"");
						}
						veccmp(0, (short unsigned *)(&n), 0, 0);
					} else {
						invaArg("\"-d\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "dh") == 0) {
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
			} else if(strcmp(arg, "fp") == 0) {
				size = sizeof(float);
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
			} else if(strcmp(arg, "t") == 0) {
				if(++args < argc) {
					t = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-t\"");
					}
				} else {
					missArg("\"-t\"");
				}
			} else if(strcmp(arg, "p") == 0) {
				if(++args < argc) {
					alpha = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || alpha < 0) {
						invaArg("\"-p\"");
					}
				} else {
					missArg("\"-p\"");
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
		return add2Matrix(*filenames, addfilename, outputfilename, noutputfilename, diffilename, targetTemplate, minCov, norm, minDepth, minLength, proxi, flag, veccmp, t);
	} else {
		makeMatrix(numFile, filenames, outputfilename, noutputfilename, diffilename, targetTemplate, minCov, alpha, norm, minDepth, minLength, proxi, flag, veccmp, methfilename, t);
	}
	
	return 0;
}
