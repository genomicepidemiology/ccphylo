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
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo calculates distances between samples based on overlaps between nucleotide count matrices created by e.g. KMA.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s).", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Target template(s)", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-md", "Minimum depth", "15");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-nm", "Normalization", "1000000");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output format", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "cos");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_dist(int argc, char *argv[]) {
	
	unsigned args, numFile, format, norm, minDepth, n;
	char *arg, *targetTemplate, **filenames, *outputfilename, *errorMsg;
	double alpha, (*veccmp)(short unsigned*, short unsigned*, int, int);
	FILE *outfile;
	Matrix *distMat;
	
	/* set defaults */
	numFile = 0;
	format = 0;
	norm = 1000000;
	minDepth = 15;
	targetTemplate = 0;
	filenames = 0;
	outputfilename = "--";
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
			} else if(strcmp(arg, "t") == 0) {
				if(++args < argc) {
					targetTemplate = argv[args];
				} else {
					missArg("\"-t\"");
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
	if(!numFile || !targetTemplate) {
		fprintf(stderr, "Missing arguments, printing helpmessage.\n");
		return helpMessage(stderr);
	}
	
	/* make ltd matrix */
	distMat = ltdMatrix_get(targetTemplate, filenames, numFile, norm, minDepth, veccmp);
	
	/* print ltd matrix */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	printphy(outfile, distMat, filenames, format);
	fclose(outfile);
	
	return 0;
}
