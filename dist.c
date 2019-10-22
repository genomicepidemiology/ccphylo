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
#include "dist.h"
#include "ltdmatrix.h"
#include "matcmp.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "str.h"
#include "unionparse.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static void makeMatrix(unsigned numFile, char **filenames, char *outputfilename, char *noutputfilename, char *targetTemplate, double minCov, double alpha, unsigned norm, unsigned minDepth, unsigned minLength, unsigned format, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	unsigned n, pos, len;
	unsigned char *include;
	FILE *outfile, *noutfile;
	FileBuff *infile, *unionfile;
	Matrix *distMat, *nDest;
	MatrixCounts *mat1;
	NucCount *mat2;
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
	
	/* init */
	mat1 = initMat(1048576, 128);
	mat2 = initNucCount(128);
	infile = setFileBuff(1048576);
	
	if(targetTemplate) {
		/* init */
		distMat = ltdMatrix_init(numFile);
		nDest = ltdMatrix_init(numFile);
		if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
			ERROR();
		}
		include = smalloc(numFile);
		memset(include, 1, numFile);
		
		/* make ltd matrix */
		ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, targetTemplate, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
		/* print ltd matrix */
		if(1 < distMat->n) {
			printphy(outfile, distMat, filenames, include, format);
			if(noutputfilename) {
				printphy(noutfile, nDest, filenames, include, format);
			}
		}
	} else if(numFile == 1) {
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
		
		/* get filenames */
		n = numFile;
		while(n--) {
			/* truncate */
			len = strlen(filenames[n]);
			if((pos = rstrpos(filenames[n], '.', len)) != -1) {
				filenames[n][pos] = 0;
				len = pos;
			}
			
			/* add .mat.gz */
			strcpy(filenames[n] + len, ".mat.gz");
			len += 7;
			
			/* check if file exists */
			if(!fileExist(filenames[n])) {
				len -= 3;
				filenames[n][len] = 0;
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
			ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, entry->target, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
			
			/* print ltd matrix */
			if(1 < distMat->n) {
				printphy(outfile, distMat, filenames, include, format);
				if(noutputfilename) {
					printphy(noutfile, nDest, filenames, include, format);
				}
			}
		}
		closeFileBuff(unionfile);
		destroyFileBuff(unionfile);
		UnionEntry_destroy(entry);
	} else {
		fprintf(stderr, "Missing input.\n");
		exit(1);
	}
	
	Matrix_destroy(distMat);
	Matrix_destroy(nDest);
	destroyMat(mat1);
	destroyNucCount(mat2);
	destroyFileBuff(infile);
	free(targetStamps);
	free(include);
	
	/* close */
	fclose(outfile);
	if(noutputfilename && strcmp(noutputfilename, outputfilename) != 0) {
		fclose(noutfile);
	}
}

static int add2Matrix(char *path, char *addfilename, char *outputfilename, char *noutputfilename, char *targetTemplate, double minCov, double alpha, unsigned norm, unsigned minDepth, unsigned minLength, unsigned format, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	int n, pos;
	double *D, *N;
	char *ptr;
	FILE *outfile, *noutfile;
	FileBuff *infile;
	MatrixCounts *mat1;
	NucCount *mat2;
	Qseqs **names;
	
	/* init */
	mat1 = initMat(1048576, 128);
	mat2 = initNucCount(128);
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
	
	/* open output and init new row(s) */
	outfile = sfopen(outputfilename, "rb+");
	if(noutputfilename) {
		noutfile = sfopen(noutputfilename, "rb+");
	} else {
		noutfile = 0;
	}
	
	/* calculate new row */
	if(ltdRow_get(D, N, mat1, mat2, infile, targetTemplate, addfilename, names, n, norm, minDepth, minLength, minCov, veccmp)) {
		fprintf(stderr, "Distance measures failed and thus the matrix was not updated.\n");
		return 1;
	}
	
	/* print updated matrix */
	printphyUpdate(outfile, ++n, addfilename, D, format);
	fclose(outfile);
	if(noutfile) {
		printphyUpdate(noutfile, n, addfilename, N, format);
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
	destroyMat(mat1);
	destroyNucCount(mat2);
	destroyFileBuff(infile);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo dist calculates distances between samples based on overlaps between nucleotide count matrices created by e.g. KMA.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-n", "Output number of nucleotides included", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-r", "Target reference", "None");
	/* here */
	//fprintf(out, "# %16s\t%-32s\t%s\n", "-a", "Add file to existing matrix", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-md", "Minimum depth", "15");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mc", "Minimum coverage", "50.0%");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum overlapping length", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nm", "Normalization", "1000000");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output format", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "cos");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-p", "Minimum lvl. of signifiacnce", "0.05");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_dist(int argc, char *argv[]) {
	
	unsigned args, numFile, format, norm, minDepth, minLength, n;
	char *arg, *targetTemplate, **filenames, *addfilename, *errorMsg;
	char *outputfilename, *noutputfilename;
	double minCov, alpha;
	double (*veccmp)(short unsigned*, short unsigned*, int, int);
	
	/* set defaults */
	numFile = 0;
	format = 1;
	norm = 1000000;
	minDepth = 15;
	minLength = 1;
	targetTemplate = 0;
	filenames = 0;
	addfilename = 0;
	outputfilename = "--";
	noutputfilename = 0;
	minCov = 0.5;
	alpha = 0.05;
	veccmp = &coscmp;
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				filenames = argv + ++args;
				numFile = 1;
				while(++args < argc && (*argv[args] != '-' || (argv[args][1] == '-' && argv[args][2] == 0))) {
					++numFile;
				}
				if(numFile == 0) {
					missArg("\"-i\"");
				}
				--args;
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
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					format = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
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
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "cos") == 0) {
						veccmp = &coscmp;
					} else if(strcmp(arg, "z") == 0) {
						veccmp = &zcmp;
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
				fprintf(stdout, "Distance calculation methods:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "# cos:\tCalculate distance between positions as the angle between the count vectors.\n");
				fprintf(stdout, "# z:\tMake consensus comparison if vectors passes a McNemar test\n");
				fprintf(stdout, "# nbc:\tCalculate the normalized Bray-Curtis dissimilarity between the count vectors.\n");
				fprintf(stdout, "# bc:\tCalculate the Bray-Curtis dissimilarity between the count vectors.\n");
				fprintf(stdout, "# nln:\tCalculate distance between positions as the normalized n-norm distance between the count vectors. Replace last \"n\" with the waned norm\n");
				fprintf(stdout, "# nlinf:\tCalculate distance between positions as the normalized l_infinity distance between the count vectors.\n");
				fprintf(stdout, "# ln:\tCalculate distance between positions as the n-norm distance between the count vectors. Replace \"n\" with the waned norm\n");
				fprintf(stdout, "# linf:\tCalculate distance between positions as the l_infinity distance between the count vectors.\n");
				fprintf(stdout, "#\n");
				return 0;
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
	
	/* set function variables */
	zcmp(0, (short unsigned *)(&alpha), 0, 0);
	
	/* check for required input */
	if(!numFile && targetTemplate) {
		numFile = 1;
	}
	
	if(addfilename && filenames) {
		return add2Matrix(*filenames, addfilename, outputfilename, noutputfilename, targetTemplate, minCov, alpha, norm, minDepth, minLength, format, veccmp);
	} else {
		makeMatrix(numFile, filenames, outputfilename, noutputfilename, targetTemplate, minCov, alpha, norm, minDepth, minLength, format, veccmp);
	}
	
	return 0;
}
