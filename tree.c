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
#include "matrix.h"
#include "nj.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tree.h"
#include "vector.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

void formTree(char *inputfilename, char *outputfilename) {
	
	/* here */
	/* generalize to output several trees */
	unsigned i, *N;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D, *Q;
	Qseqs **names;
	Vector *sD;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	names = smalloc(32 * sizeof(Qseqs *));
	i = 33;
	names += 32;
	while(--i) {
		*--names = setQseqs(32);
	}
	D = ltdMatrix_init(32);
	Q = ltdMatrix_init(32);
	sD = vector_init(32);
	N = smalloc(32 * sizeof(unsigned));
	
	/* set */
	openAndDetermine(infile, inputfilename);
	names = loadPhy(D, names, infile);
	
	if(1 < D->n) {
		/* make tree */
		N = nj(D, Q, sD, N, names);
		
		/* output tree */
		fprintf(outfile, ">%s%s;\n", "nj", (*names)->seq);
	}
	
	/* clean */
	fclose(outfile);
	closeFileBuff(infile);
	free(N);
	Matrix_destroy(D);
	Matrix_destroy(Q);
	vector_destroy(sD);
	free(names);
	destroyFileBuff(infile);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo forms tree(s) in newick format given a set of phylip distance matrices.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file.", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-m", "Tree construction method.", "nj");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mh", "Help on option \"-m\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tree(int argc, char *argv[]) {
	
	unsigned args, method;
	char *arg, *inputfilename, *outputfilename;
	
	/* set defaults */
	method = 1;
	inputfilename = "--";
	outputfilename = "--";
	
	if((args = 1) == argc) {
		fprintf(stderr, "Missing arguments, printing helpmessage.\n");
		return helpMessage(stderr);
	}
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
			} else if(strcmp(arg, "m") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "nj") == 0) {
						method = 1;
					} else {
						invaArg("\"-m\"");
					}
				} else {
					missArg("\"-m\"");
				}
			} else if(strcmp(arg, "mh") == 0) {
				fprintf(stdout, "Tree construction methods:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "# nj:\tNeighbour-Joining.\n");
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
	
	if(method == 1) {
		formTree(inputfilename, outputfilename);
	}
	
	return 0;
}
