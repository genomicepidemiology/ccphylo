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
#include "dbparse.h"
#include "filebuff.h"
#include "hashmapstr.h"
#include "pherror.h"
#include "resparse.h"
#include "union.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

HashMapStr * unionRes(char **filenames, int numFile, char *outputfilename, double minCov, double minDepth, unsigned minLength) {
	
	unsigned n;
	FileBuff *infile;
	HashMapStr *entries;
	ResEntry *entry;
	
	/* init */
	entries = HashMapStr_init(128);
	entry = ResEntry_init(128);
	infile = setFileBuff(1048576);
	minLength *= 100;
	
	/* iterate files */
	for(n = 0; n < numFile; ++n) {
		openAndDetermine(infile, *filenames++);
		if(FileBuffValidateHeader(infile)) {
			fprintf(stderr, "Malformed res file:\t%s\n", *--filenames);
			exit(1);
		}
		while(FileBuffGetEntry(infile, entry)) {
			if(minCov <= entry->Template_Coverage && minDepth <= entry->Depth && minLength <= entry->Template_length * entry->Template_Coverage) {
				HashMapStr_add(entries, entry->Template->seq, n);
			}
		}
		closeFileBuff(infile);
	}
	
	/* clean up */
	destroyFileBuff(infile);
	
	return entries;
}

int unionResPrint(char **filenames, int numFile, char *outputfilename, double minCov, double minDepth, unsigned minLength) {
	
	FILE *outfile;
	HashMapStr *entries;
	
	/* init */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	
	/* get union */
	entries = unionRes(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	
	/* print results */
	numFile = HashMapStr_print(entries, outfile);
	
	/* clean up */
	fclose(outfile);
	HashMapStr_destroy(entries);
	
	return numFile;
}

int unionResOrderPrint(char **filenames, int numFile, char *outputfilename, char *dbfilename, double minCov, double minDepth, unsigned minLength) {
	
	int nc, num;
	unsigned *ptr;
	FILE *outfile, *templatefile;
	BucketStr *node;
	HashMapStr *entries;
	Qseqs *templatename;
	
	/* init */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	templatename = setQseqs(strlen(dbfilename) + 64);
	templatename->len = sprintf((char *) templatename->seq, "%s.name", dbfilename);
	templatefile = sfopen((char *) templatename->seq, "rb");
	
	/* get union */
	entries = unionRes(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	
	/* print db ordered results */
	nc = 0;
	while(entries->n && nameLoad(templatename, templatefile)) {
		if((node = HashMapStr_get(entries, templatename->seq)) && 0 < (num = node->num)) {
			nc += fprintf(outfile, "%s\t%d", templatename->seq, ++num);
			++num;
			ptr = node->uList - 1;
			while(--num) {
				nc += fprintf(outfile, "\t%s", filenames[*++ptr]);
			}
			nc += fprintf(outfile, "\n");
			/* destroy node */
			free(node->str);
			free(node->uList);
			free(node);
		}
	}
	
	/* clean */
	fclose(outfile);
	fclose(templatefile);
	HashMapStr_destroy(entries);
	destroyQseqs(templatename);
	
	return nc;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo union finds the union between templates in res files created by e.g. KMA.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s).", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Print ordered wrt. template DB filename", "none");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-md", "Minimum depth", "1.0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mc", "Minimum coverage", "50.0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum consensus length", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_union(int argc, char *argv[]) {
	
	unsigned args, numFile, minLength;
	double minCov, minDepth;
	char *arg, **filenames, *outputfilename, *templatefilename, *errorMsg;
	
	numFile = 0;
	filenames = 0;
	outputfilename = "--";
	templatefilename = 0;
	minDepth = 1;
	minCov = 50.0;
	minLength = 1;
	
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
					templatefilename = argv[args];
				} else {
					missArg("\"-t\"");
				}
			} else if(strcmp(arg, "md") == 0) {
				if(++args < argc) {
					minDepth = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0) {
						invaArg("\"-md\"");
					}
				} else {
					missArg("\"-md\"");
				}
			} else if(strcmp(arg, "mc") == 0) {
				if(++args < argc) {
					minCov = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0) {
						invaArg("\"-mc\"");
					}
				} else {
					missArg("\"-mc\"");
				}
			} else if(strcmp(arg, "ml") == 0) {
				if(++args < argc) {
					minLength = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-mc\"");
					}
				} else {
					missArg("\"-mc\"");
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
	
	/* check for required input */
	if(!numFile || !filenames) {
		fprintf(stderr, "Missing arguments, printing helpmessage.\n");
		return helpMessage(stderr);
	}
	
	if(templatefilename) {
		unionResOrderPrint(filenames, numFile, outputfilename, templatefilename, minCov, minDepth, minLength);
	} else {
		unionResPrint(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	}
	
	return 0;
}
