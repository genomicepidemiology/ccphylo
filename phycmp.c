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
#include "bytescale.h"
#include "cmdline.h"
#include "distcmp.h"
#include "filebuff.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"

int entriescmp(Qseqs **names1, Qseqs **names2, int n) {
	
	int cmp;
	
	if(!n) {
		return 0;
	} else if(!names1 && !names2) {
		return 0;
	} else if(!names1 || !names2) {
		return 1;
	}
	
	--names1;
	--names2;
	++n;
	while(--n) {
		cmp = strcmp((char*)((*++names1)->seq), (char*)((*++names2)->seq));
		if(cmp) {
			return cmp;
		}
	}
	
	return 0;
}

void phyfilecmp(char **inputfilenames, int filenum, char *outputfilename, int flag, char sep, char quotes) {
	
	int i;
	long n;
	double d;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D1, *D2;
	Qseqs **names1, **names2, *header1, *header2;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	i = 32;
	D1 = ltdMatrix_init(i);
	D2 = ltdMatrix_init(i);
	header1 = setQseqs(64);
	header2 = setQseqs(64);
	names1 = smalloc(i * sizeof(Qseqs *));
	names1 += i;
	names2 = smalloc(i * sizeof(Qseqs *));
	names2 += i;
	++i;
	while(--i) {
		*--names1 = setQseqs(4);
		*--names2 = setQseqs(4);
	}
	
	/* load matrices */
	openAndDetermine(infile, *inputfilenames);
	names1 = loadPhy(D1, names1, header1, infile, sep, quotes);
	if(filenum != 1) {
		closeFileBuff(infile);
		openAndDetermine(infile, inputfilenames[1]);
	}
	names2 = loadPhy(D2, names2, header2, infile, sep, quotes);
	closeFileBuff(infile);
	destroyFileBuff(infile);
	
	/* validate matrices for comparison */
	if(!D1->n || !D2->n) {
		fprintf(stderr, "Missing matrix\n");
		exit(1);
	} else if(D1->n != D2->n) {
		fprintf(stderr, "Matrices differ in size.\n");
		exit(1);
	} else if(entriescmp(names1, names2, D1->n)) {
		fprintf(stderr, "Matrices has different entries.\n");
		exit(1);
	}
	
	/* get comparisons */
	n = D1->n;
	n *= (n - 1);
	n >>= 1;
	if(flag & 1) { /* cos */
		d = D1->mat ? coscmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? coscmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? coscmp_s(*(D1->smat), *(D2->smat), n) : coscmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "cos:\t%f\n", d);
	}
	if(flag & 2) { /* chi2 */
		d = D1->mat ? chi2cmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? chi2cmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? chi2cmp_s(*(D1->smat), *(D2->smat), n) : chi2cmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "chi2:\t%f\n", d);	
	}
	if(flag & 4) { /* bray-curtis */
		d = D1->mat ? bccmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? bccmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? bccmp_s(*(D1->smat), *(D2->smat), n) : bccmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "bc:\t%f\n", d);		
	}
	if(flag & 8) { /* l1 */
		d = D1->mat ? l1cmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? l1cmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? l1cmp_s(*(D1->smat), *(D2->smat), n) : l1cmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "l1:\t%f\n", d);			
	}
	if(flag & 16) { /* l2 */
		d = D1->mat ? l2cmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? l2cmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? l2cmp_s(*(D1->smat), *(D2->smat), n) : l2cmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "l2:\t%f\n", d);				
	}
	if(flag & 32) { /* linf */
		d = D1->mat ? linfcmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? linfcmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? linfcmp_s(*(D1->smat), *(D2->smat), n) : linfcmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "linf:\t%f\n", d);				
	}
	if(flag & 64) { /* p */
		d = D1->mat ? pearcmp_d(*(D1->mat), *(D2->mat), n) : D1->fmat ? pearcmp_f(*(D1->fmat), *(D2->fmat), n) : D1->smat ? pearcmp_s(*(D1->smat), *(D2->smat), n) : pearcmp_b(*(D1->bmat), *(D2->bmat), n);
		fprintf(outfile, "p:\t%f\n", d);				
	}
	
	/* clean up */
	fclose(outfile);
	i = D1->n;
	names1 += i;
	names2 += i;
	++i;
	while(--i) {
		destroyQseqs(*--names1);
		destroyQseqs(*--names2);
	}
	free(names1);
	free(names2);
	destroyQseqs(header1);
	destroyQseqs(header2);
	Matrix_destroy(D1);
	Matrix_destroy(D2);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "# CCPhylo phycmp compares two distance matrices in phylip format.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file(s)", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'p', "float_precision", "Float precision on distance matrix", "False / double");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 's', "short_precision", "Short precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "byte_precision", "Byte precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'H', "mmap", "Allocate matrix on the disk", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tmp", "Set directory for temporary files", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_phycmp(int argc, char **argv) {
	
	const char *stdstream = "-";
	int size, len, offset, args, flag, filenum;
	char **Arg, *arg, **filenames, *outputfilename, *tmp, opt, sep, quotes;
	
	/* set defaults */
	size = sizeof(double);
	flag = 1;
	filenum = 0;
	filenames = 0;
	outputfilename = (char *)(stdstream);
	tmp = 0;
	sep = '\t';
	quotes = '\0';
	
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
					filenum = getArgListLen(&Arg, &args);
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "flag") == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(cmdcmp(arg, "flag_help") == 0) {
					flag = -1;
				} else if(cmdcmp(arg, "float_precision") == 0) {
					size = sizeof(float);
				} else if(cmdcmp(arg, "short_precision") == 0) {
					size = sizeof(short unsigned);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "short_precision");
				} else if(cmdcmp(arg, "byte_precision") == 0) {
					size = sizeof(unsigned char);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "byte_precision");
				} else if(cmdcmp(arg, "mmap") == 0) {
					ltdMatrix_init = &ltdMatrixMinit;
				} else if(cmdcmp(arg, "tmp") == 0) {
					tmp = getArgDie(&Arg, &args, len + offset, "tmp");
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
						filenum = getArgListLen(&Arg, &args);
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
					} else if(opt == 'p') {
						size = sizeof(float);
					} else if(opt == 's') {
						size = sizeof(short unsigned);
						ByteScale = getdDefArg(&Arg, &args, len, ByteScale, "s");
						opt = 0;
					} else if(opt == 'b') {
						size = sizeof(unsigned char);
						ByteScale = getdDefArg(&Arg, &args, len, ByteScale, "b");
						opt = 0;
					} else if(opt == 'H') {
						ltdMatrix_init = &ltdMatrixMinit;
					} else if(opt == 'T') {
						tmp = getArgDie(&Arg, &args, len, "T");
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
		filenum = args;
	}
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "# Distance calculation methods:\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "# 1\tcos: Calculate cosine distance between vectors.\n");
		fprintf(stdout, "# 2\tchi2: Calculate the chi square distance\n");
		fprintf(stdout, "# 4\tbc: Calculate the Bray-Curtis dissimilarity between vectors.\n");
		fprintf(stdout, "# 8\tl1: Calculate distance between vectors as the 1-norm distance between the count vectors.\n");
		fprintf(stdout, "# 16\tl2: Calculate distance between vectors as the 2-norm distance between the count vectors.\n");
		fprintf(stdout, "# 32\tlinf: Calculate distance between vectors as the l_infinity distance between the count vectors.\n");
		fprintf(stdout, "# 64\tp: Calculate Pearsons correlation between vectors.\n");
		fprintf(stdout, "#\n");
		return 0;
	}
	
	/* tmp dir */
	if(tmp) {
		tmpF(tmp);
	}
	
	/* set precision */
	ltdMatrixInit(-size);
	ltdMatrixMinit(-size);
	
	/* make full phy */
	phyfilecmp(filenames, filenum, outputfilename, flag, sep, quotes);
	
	return 0;
}
