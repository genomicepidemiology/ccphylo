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
#include "cmdline.h"
#include "dbparse.h"
#include "filebuff.h"
#include "hashmapstr.h"
#include "pherror.h"
#include "resparse.h"
#include "seq2fasta.h"
#include "union.h"

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
	
	int nc;
	FILE *outfile;
	HashMapStr *entries;
	
	/* init */
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	
	/* get union */
	entries = unionRes(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	
	/* print tested filenames */
	fprintf(outfile, "%d", numFile);
	nc = numFile + 1;
	while(--nc) {
		fprintf(outfile, "\t%s", *filenames++);
	}
	fprintf(outfile, "\n");
	
	/* print results */
	nc = HashMapStr_print(entries, outfile);
	
	/* clean up */
	fclose(outfile);
	HashMapStr_destroy(entries);
	
	return nc;
}

int unionResOrderPrint(char **filenames, int numFile, char *outputfilename, char *dbfilename, char *reffilename, double minCov, double minDepth, unsigned minLength) {
	
	int nc, num, tnum, *template_lengths, seqlist[2];
	unsigned *ptr;
	char *templatefilename;
	FILE *outfile, *templatefile, *reffile;
	BucketStr *node;
	HashMapStr *entries;
	Qseqs *templatename;
	
	/* init */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	templatename = setQseqs(64);
	tnum = strlen(dbfilename);
	templatefilename = smalloc(tnum + 64);
	sprintf(templatefilename, "%s.name", dbfilename);
	templatefile = sfopen(templatefilename, "rb");
	templatefilename[tnum] = 0;
	reffile = 0;
	tnum = 1;
	
	/* get union */
	entries = unionRes(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	
	/* print tested filenames */
	if(reffilename) {
		reffile = sfopen(reffilename, "wb");
		fprintf(outfile, "%d\t%s", numFile + 1, reffilename);
		seqlist[0] = 1;
		template_lengths = getLengths(templatefilename);
	} else {
		fprintf(outfile, "%d", numFile);
		seqlist[0] = 0;
		template_lengths = 0;
	}
	nc = numFile + 1;
	while(--nc) {
		fprintf(outfile, "\t%s", *filenames++);
	}
	fprintf(outfile, "\n");
	
	/* print db ordered results */
	nc = 0;
	while(entries->n && nameLoad(templatename, templatefile)) {
		if((node = HashMapStr_get(entries, templatename->seq)) && 0 < (num = node->num)) {
			if(reffile) {
				seqlist[1] = tnum;
				printFastaList(reffile, templatefilename, template_lengths, seqlist);
				fflush(reffile);
				
				nc += fprintf(outfile, "%s\t%d\t%d", templatename->seq, (num += 2), 0);
				ptr = node->uList - 1;
				while(--num) {
					nc += fprintf(outfile, "\t%d", *++ptr + 1);
				}
			} else {
				nc += fprintf(outfile, "%s\t%d", templatename->seq, ++num);
				++num;
				ptr = node->uList - 1;
				while(--num) {
					nc += fprintf(outfile, "\t%d", *++ptr);
				}
			}
			nc += fprintf(outfile, "\n");
			/* destroy node */
			free(node->str);
			free(node->uList);
			free(node);
		}
		++tnum;
	}
	/* clean */
	fclose(outfile);
	fclose(templatefile);
	HashMapStr_destroy(entries);
	destroyQseqs(templatename);
	
	/* get references */
	if(reffile) {
		fclose(reffile);
		free(template_lengths);
	}
	
	return nc;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo union finds the union between templates in res files created by e.g. KMA.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file(s)", "None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'B', "database", "Print ordered wrt. template DB filename", "None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'r', "reference_file", "Create reference fasta file", "None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'E', "min_depth", "Minimum depth", "15");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'C', "min_cov", "Minimum coverage", "50.0%");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'L', "min_len", "Minimum overlapping length", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	db	B	database
	r	r	reference_file
	md	E	min_depth
	mc	C	min_cov
	ml	L	min_len
	h	h	help
	
	*/
}

int main_union(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, numFile, minLength, len, offset;
	double minCov, minDepth;
	char **Arg, *arg, **filenames, *outputfilename, *templatefilename;
	char *reffilename, opt;
	
	numFile = 0;
	filenames = 0;
	outputfilename = (char *)(stdstream);
	templatefilename = 0;
	reffilename = 0;
	minDepth = 1;
	minCov = 50.0;
	minLength = 1;
	
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
					filenames = getArgListDie(&Arg, &args, len + offset, "input");
					numFile = getArgListLen(&Arg, &args);
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "database") == 0) {
					templatefilename = getArgDie(&Arg, &args, len + offset, "database");
				} else if(cmdcmp(arg, "reference_file") == 0) {
					reffilename = getArgDie(&Arg, &args, len + offset, "reference_file");
				} else if(cmdcmp(arg, "min_depth") == 0) {
					minDepth = getdArg(&Arg, &args, len + offset, "min_depth");
				} else if(cmdcmp(arg, "min_cov") == 0) {
					minCov = getdArg(&Arg, &args, len + offset, "min_cov");
				} else if(cmdcmp(arg, "min_len") == 0) {
					minLength = getNumArg(&Arg, &args, len + offset, "min_len");
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
						filenames = getArgListDie(&Arg, &args, len, "i");
						numFile = getArgListLen(&Arg, &args);
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'B') {
						templatefilename = getArgDie(&Arg, &args, len, "B");
						opt = 0;
					} else if(opt == 'r') {
						reffilename = getArgDie(&Arg, &args, len, "r");
						opt = 0;
					} else if(opt == 'E') {
						minDepth = getdArg(&Arg, &args, len, "E");
						opt = 0;
					} else if(opt == 'C') {
						minCov = getdArg(&Arg, &args, len, "C");
						opt = 0;
					} else if(opt == 'L') {
						minLength = getNumArg(&Arg, &args, len, "L");
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
		filenames = Arg;
		numFile = args;
	}
	
	/* check for required input */
	if(!numFile || !filenames) {
		fprintf(stderr, "Missing arguments, printing helpmessage.\n");
		return helpMessage(stderr);
	}
	
	if(templatefilename) {
		unionResOrderPrint(filenames, numFile, outputfilename, templatefilename, reffilename, minCov, minDepth, minLength);
	} else if(reffilename) {
		fprintf(stderr, "Database is needed in order to reconstruct the reference(s).\n");
		exit(1);
	} else {
		unionResPrint(filenames, numFile, outputfilename, minCov, minDepth, minLength);
	}
	
	return 0;
}
