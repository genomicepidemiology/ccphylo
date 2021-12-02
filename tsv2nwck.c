/* Philip T.L.C. Clausen Jul 2021 plan@dtu.dk */

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
#include "dat.h"
#include "datclust.h"
#include "distcmp.h"
#include "filebuff.h"
#include "hclust.h"
#include "nwck.h"
#include "pherror.h"
#include "qseqs.h"
#include "tmp.h"
#include "tsv.h"
#include "tsv2nwck.h"
#include "vector.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

int tsv2nwck(char *inputfilename, char *outputfilename, unsigned char sep) {
	
	int i, *P;
	FILE *outfile;
	FileBuff *infile;
	Dat *Dmat;
	Qseqs **names;
	Vector *Q;
	
	/* init */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	
	/* load tsv */
	if(!(Dmat = loadTsv(infile, sep))) {
		fprintf(stderr, "Input matrix contained zero rows.\n");
		return 0;
	}
	 closeFileBuff(infile);
	 destroyFileBuff(infile);
	
	/* Get Q */
	Q = vector_init(Dmat->m);
	P = smalloc(Dmat->m * sizeof(int));
	initQ_Dmat(Dmat, Q, P);
	Dat_destroy(Dmat);
	
	/* Init newick clustering format */
	names = smalloc(Dmat->m * sizeof(Qseqs *));
	names += Dmat->m;
	i = Dmat->m;
	while(i--) {
		*--names = setQseqs(10);
		(*names)->len = sprintf((char *)((*names)->seq), "%d", i);
	}
	
	/* Cluster Q */
	tclust(Q, P, names);
	fprintf(outfile, "%s;\n", (*names)->seq);
	
	/* clean up */
	fclose(outfile);
	vector_destroy(Q);
	free(P);
	names += Dmat->m;
	i = Dmat->m;
	while(i--) {
		destroyQseqs(*--names);
	}
	free(names);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo tsv2phy converts tsv files to phylip distance files.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-s", "Separator", "\\t");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "cos");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on matrix", "double");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-sp", "Short precision on matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tsv2nwck(int argc, char *argv[]) {
	
	int args, size;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	unsigned char sep;
	double exponent;
	
	/* init */
	size = sizeof(double);
	inputfilename = "--";
	outputfilename = "--";
	sep = '\t';
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					inputfilename = argv[args];
				} else {
					missArg("\"-i\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "s") == 0) {
				if(argc <= ++args) {
					missArg("\"-s\"");
				} else if((sep = argv[args][0]) == 0 || argv[args][1] != 0) {
					invaArg("\"-s\"");
				}
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "cos") == 0) {
						distcmp_d = &coscmp_d;
						distcmp_f = &coscmp_f;
						distcmp_b = &coscmp_b;
					} else if(strcmp(arg, "chi2") == 0) {
						distcmp_d = &chi2cmp_d;
						distcmp_f = &chi2cmp_f;
						distcmp_b = &chi2cmp_b;
					} else if(strcmp(arg, "bc") == 0) {
						distcmp_d = &bccmp_d;
						distcmp_f = &bccmp_f;
						distcmp_b = &bccmp_b;
					} else if(strcmp(arg, "l1") == 0) {
						distcmp_d = &l1cmp_d;
						distcmp_f = &l1cmp_f;
						distcmp_b = &l1cmp_b;
					} else if(strcmp(arg, "l2") == 0) {
						distcmp_d = &l2cmp_d;
						distcmp_f = &l2cmp_f;
						distcmp_b = &l2cmp_b;
					} else if(strcmp(arg, "linf") == 0) {
						distcmp_d = &linfcmp_d;
						distcmp_f = &linfcmp_f;
						distcmp_b = &linfcmp_b;
					} else if(*arg == 'l') {
						distcmp_d = &lncmp_d;
						distcmp_f = &lncmp_f;
						distcmp_b = &lncmp_b;
						exponent = strtod(arg + 1, &errorMsg);
						if(*errorMsg != 0) {
							invaArg("\"-d ln\"");
						}
						distcmp_d(0, &exponent, 0);
						distcmp_f(0, (float *)(&exponent), 0);
						distcmp_b(0, (unsigned char *)(&exponent), 0);
					} else if(strcmp(arg, "p") == 0) {
						distcmp_d = &pearcmp_d;
						distcmp_f = &pearcmp_f;
						distcmp_b = &pearcmp_b;
					} else {
						invaArg("\"-d\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "dh") == 0) {
				fprintf(stdout, "# Distance calculation methods:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "# cos:\tCalculate cosine distance between vectors.\n");
				fprintf(stdout, "# chi2:\tCalculate the chi square distance\n");
				fprintf(stdout, "# bc:\tCalculate the Bray-Curtis dissimilarity between vectors.\n");
				fprintf(stdout, "# ln:\tCalculate distance between vectors as the n-norm distance between the count vectors. Replace \"n\" with the waned norm\n");
				fprintf(stdout, "# linf:\tCalculate distance between vectors as the l_infinity distance between the count vectors.\n");
				fprintf(stdout, "# p:\tCalculate Pearsons correlation between vectors.\n");
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
				Dat_init = &DatMinit;
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
	DatInit(-size, 0);
	DatMinit(-size, 0);
	
	/* set ptrs */
	pairQ = &minQ;
	formLastNodePtr = &formLastBiNode;
	
	return tsv2nwck(inputfilename, outputfilename, sep);
}
