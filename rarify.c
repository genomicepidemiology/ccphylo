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
#include "cmdline.h"
#include "filebuff.h"
#include "matparse.h"
#include "pherror.h"
#include "rarify.h"

int rarify(char *inputfilename, char *outputfilename, long unsigned nf, long unsigned rf) {
	
	unsigned i;
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
	if(*outputfilename == '-' && outputfilename[1] == 0) {
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
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'A', "fragment_amount", "Total number of fragments in sample", "0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'R', "rarification_factor", "Rarification factor", "10000000");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	nf	A	fragment_amount
	rf	R	rarification_factor
	h	h	help
	*/
}

int main_rarify(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, len, offset;
	long unsigned nf, rf;
	char **Arg, *arg, *inputfilename, *outputfilename, opt;
	
	nf = 0;
	rf = 10000000;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	
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
					inputfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "fragment_amount") == 0) {
					nf = getNumArg(&Arg, &args, len + offset, "fragment_amount");
				} else if(cmdcmp(arg, "rarification_factor") == 0) {
					rf = getNumArg(&Arg, &args, len + offset, "rarification_factor");
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
						inputfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'A') {
						nf = getNumArg(&Arg, &args, len, "A");
						opt = 0;
					} else if(opt == 'R') {
						rf = getNumArg(&Arg, &args, len, "R");
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
		inputfilename = *Arg;
		if(--args) {
			nonOptError();
		}
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
