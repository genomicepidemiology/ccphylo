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
#include "tmp.h"
#include "tree.h"
#include "vector.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

void formTree(char *inputfilename, char *outputfilename, int thread_num) {
	
	unsigned i, *N;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	Vector *sD;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	i = 32;
	D = ltdMatrix_init(i);
	sD = vector_init(i);
	N = smalloc(i * sizeof(unsigned));
	names = smalloc(i * sizeof(Qseqs *));
	names += i;
	++i;
	while(--i) {
		*--names = setQseqs(64);
	}
	
	/* set */
	openAndDetermine(infile, inputfilename);
	
	/* generate trees */
	while((names = loadPhy(D, names, header, infile)) && D->n) {
		if(2 < D->n) {
			/* make tree */
			N = nj_thread(D, sD, N, names, thread_num);
			
			/* output tree */
			if(header->len) {
				fprintf(outfile, ">%s%s;\n", header->seq, (*names)->seq);
			} else {
				fprintf(outfile, "%s;\n", (*names)->seq);
			}
		} else if(D->n == 2) {
			/* output tree */
			if(header->len) {
				fprintf(outfile, ">%s(%s,%s:%.2f);\n", header->seq, (*names)->seq, names[1]->seq, **(D->mat));
			} else {
				fprintf(outfile, "(%s,%s:%.2f);\n", (*names)->seq, names[1]->seq, **(D->mat));
			}
		}
	}
	
	/* clean */
	fclose(outfile);
	closeFileBuff(infile);
	free(N);
	Matrix_destroy(D);
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tree(int argc, char *argv[]) {
	
	unsigned args, thread_num;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	
	/* set defaults */
	thread_num = 1;
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
			} else if(strcmp(arg, "m") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "nj") == 0) {
						updateDptr = &updateD;
						minDist = &initQ;
						minDist_thread = &initQ_thread;
					} else if(strcmp(arg, "upgma") == 0) {
						updateDptr = &updateD_UPGMA;
						minDist = &minD;
						minDist_thread = &minD_thread;
					} else if(strcmp(arg, "cf") == 0) {
						updateDptr = &updateD_CF;
						minDist = &minD;
						minDist_thread = &minD_thread;
					} else if(strcmp(arg, "ff") == 0) {
						updateDptr = &updateD_FF;
						minDist = &minD;
						minDist_thread = &minD_thread;
					} else {
						invaArg("\"-m\"");
					}
				} else {
					missArg("\"-m\"");
				}
			} else if(strcmp(arg, "mh") == 0) {
				fprintf(stdout, "# Tree construction methods:\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "# %-8s\t%s\n", "nj", "Neighbour-Joining");
				fprintf(stdout, "# %-8s\t%s\n", "upgma", "UPGMA");
				fprintf(stdout, "# %-8s\t%s\n", "cf", "K-means Closest First");
				fprintf(stdout, "# %-8s\t%s\n", "ff", "K-means Furthest First");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "t") == 0) {
				if(++args < argc) {
					thread_num = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-t\"");
					}
				} else {
					missArg("\"-t\"");
				}
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
	
	/* make tree */
	formTree(inputfilename, outputfilename, thread_num);
	
	return 0;
}
