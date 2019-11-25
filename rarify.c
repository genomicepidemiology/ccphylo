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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "filebuff.h"
#include "matparse.h"
#include "pherror.h"
#include "rarify.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

int rarify(char *inputfilename, char *outputfilename, long unsigned nf, long unsigned rf) {
	
	unsigned i, pos;
	long unsigned count, remainder;
	short unsigned *counts;
	FILE *outfile;
	FileBuff *infile;
	NucCount *mat;
	
	/* init */
	mat = initNucCount(128);
	infile = setFileBuff(1048576);
	
	/* open files */
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	
	/* rarify matrix */
	remainder = 0;
	while(FileBuffGetRow(infile, mat)) {
		if(mat->ref) {
			/* rarify counts */
			counts = mat->counts + 6;
			pos = 0;
			i = 7;
			while(--i) {
				if((count = *--counts)) {
					/* rarify */
					count *= rf;
					remainder += count % nf;
					count /= nf;
					
					/* check remainder */
					if(rf <= remainder) {
						count += remainder / rf;
						remainder %= rf;
					}
					
					/* update */
					*counts = count;
				}
			}
			
			/* output new counts */
			fprintf(outfile, "%c\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", mat->ref, *counts, counts[1], counts[2], counts[3], counts[4], counts[5]);	
		} else if(*(mat->name)) {
			fprintf(outfile, "#%s\n", mat->name);
		} else {
			fprintf(outfile, "\n");	
		}
	}
	closeFileBuff(infile);
	destroyFileBuff(infile);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo rarify rarifies an KMA matrix.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nf", "Total number of fragments in sample", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-rf", "Rarification factor", "10000000");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_rarify(int argc, char *argv[]) {
	
	int args;
	long unsigned nf, rf;
	char *arg, *inputfilename, *outputfilename;
	
	nf = 0;
	rf = 10000000;
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
					missArg("\"-i\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "nf") == 0) {
				if(++args < argc) {
					nf = strtoul(argv[args], &arg, 10);
					if(*arg != 0) {
						invaArg("\"-nf\"");
					}
				} else {
					missArg("\"-nf\"");
				}
			} else if(strcmp(arg, "rf") == 0) {
				if(++args < argc) {
					rf = strtoul(argv[args], &arg, 10);
					if(*arg != 0) {
						invaArg("\"-rf\"");
					}
				} else {
					missArg("\"-rf\"");
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
	
	/* insuffient input */
	if(!nf) {
		fprintf(stderr, "Missing argument:\t\"-nf\"\n");
		return helpMessage(stderr);
	}
	
	/* rarify */
	rarify(inputfilename, outputfilename, nf, rf);
	
	return 0;
}
