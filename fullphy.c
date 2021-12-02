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
#include <time.h>
#include "bytescale.h"
#include "filebuff.h"
#include "fullphy.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

void formFullPhy(char *inputfilename, char *outputfilename, int flag) {
	
	int i;
	unsigned char **sNames;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	time_t t0, t1;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	i = 32;
	D = ltdMatrix_init(i);
	names = smalloc(i * sizeof(Qseqs *));
	names += i;
	++i;
	while(--i) {
		*--names = setQseqs(4);
	}
	sNames = smalloc(i * sizeof(unsigned char *));
	
	/* set */
	openAndDetermine(infile, inputfilename);
	
	/* generate trees */
	t0 = clock();
	while((names = loadPhy(D, names, header, infile)) && D->n) {
		t1 = clock();
		fprintf(stderr, "# Total time used loading matrix: %.2f s.\n", difftime(t1, t0) / 1000000);
		t0 = t1;
		sNames = realloc(sNames, D->size * sizeof(unsigned char *));
		if(!sNames) {
			ERROR();
		} else {
			i = D->n;
			while(i--) {
				sNames[i] = names[i]->seq;
			}
		}
		printfullphy(outfile, D, (char **) sNames, flag);
		t1 = clock();
		fprintf(stderr, "# Total time outputting full matrix: %.2f s.\n", difftime(t1, t0) / 1000000);
		t0 = t1;
	}
	
	/* clean */
	fclose(outfile);
	closeFileBuff(infile);
	Matrix_destroy(D);
	free(names);
	destroyFileBuff(infile);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo forms tree(s) in newick format given a set of phylip distance matrices.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file.", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-sp", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_fullphy(int argc, char *argv[]) {
	
	int size;
	unsigned args, flag;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	
	/* set defaults */
	size = sizeof(double);
	flag = 1;
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
	
	/* make full phy */
	formFullPhy(inputfilename, outputfilename, flag);
	
	return 0;
}
