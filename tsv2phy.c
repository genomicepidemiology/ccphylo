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
#include "cmdline.h"
#include "dat.h"
#include "distcmp.h"
#include "filebuff.h"
#include "pherror.h"
#include "tmp.h"
#include "tsv.h"
#include "tsv2phy.h"
#define dfprintf(fPtr, str, d) if(fprintf(fPtr, str, d) < 0) {if(errno) {ERROR();} else {fprintf(stderr, "Could not extend file.\n"); exit(1);}}
#define ddfprintf(fPtr, str, precision, d) if(fprintf(fPtr, str, precision, d) < 0) {if(errno) {ERROR();} else {fprintf(stderr, "Could not extend file.\n"); exit(1);}}

int tsv2phy(char *inputfilename, char *outputfilename, unsigned format, unsigned char sep, int precision) {
	
	int m, n, i, j;
	double **Dptr, *Di, **Dj;
	float **Dfptr, *Dfi, **Dfj;
	short unsigned **Dsptr, *Dsi, **Dsj;
	unsigned char **Dbptr, *Dbi, **Dbj;
	FILE *outfile;
	FileBuff *infile;
	Dat *Dmat;
	
	/* init */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == 0) {
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
	
	/* print matrix */
	dfprintf(outfile, "%10d", Dmat->m);
	m = Dmat->m;
	n = Dmat->n;
	Dptr = Dmat->mat;
	Dfptr = Dmat->fmat;
	Dsptr = Dmat->smat;
	Dbptr = Dmat->bmat;
	i = -1;
	while(++i < m) {
		/* print entry name */
		if(format & 1) {
			dfprintf(outfile, "\n%d", i);
		} else {
			dfprintf(outfile, "\n%-10d", i);
		}
		
		/* print distance */
		j = i + 1;
		if(Dptr) {
			Di = Dptr[i];
			Dj = Dptr - 1;
			while(--j) {
				ddfprintf(outfile, "\t%.*g", precision, distcmp_d(Di, *++Dj, n));
			}
		} else if(Dfptr) {
			Dfi = Dfptr[i];
			Dfj = Dfptr - 1;
			while(--j) {
				ddfprintf(outfile, "\t%.*g", precision, distcmp_f(Dfi, *++Dfj, n));
			}
		} else if(Dsptr) {
			Dsi = Dsptr[i];
			Dsj = Dsptr - 1;
			while(--j) {
				ddfprintf(outfile, "\t%.*g", precision, distcmp_s(Dsi, *++Dsj, n));
			}
		} else {
			Dbi = Dbptr[i];
			Dbj = Dbptr - 1;
			while(--j) {
				ddfprintf(outfile, "\t%.*g", precision, distcmp_b(Dbi, *++Dbj, n));
			}
		}
	}
	dfprintf(outfile, "%c", '\n');
	
	/* clean up */
	sfclose(outfile);
	Dat_destroy(Dmat);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo tsv2phy converts tsv files to phylip distance files.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'x', "print_precision", "Floating point print precision", "9");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'd', "distance", "Distance method", "cos");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'D', "distance_help", "Help on option \"-d\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'p', "float_precision", "Float precision on distance matrix", "False / double");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 's', "short_precision", "Short precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "byte_precision", "Byte precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'H', "mmap", "Allocate matrix on the disk", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tmp", "Set directory for temporary files", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	s	S	separator
	d	d	distance
	dh	D	distance_help
	f	f	flag
	fh	F	flag_help
	fp	p	float_precision
	sp	s	short_precision
	bp	b	byte_precision
	mm	H	mmap
	tmp	T	tmp
	h	h	help
	*/
}

int main_tsv2phy(int argc, char **argv) {
	
	const char *stdstream = "-";
	int size, len, offset, args, flag, precision;
	char **Arg, *arg, *inputfilename, *outputfilename, *tmp, *method;
	char *errorMsg, opt;
	unsigned char sep;
	double exponent;
	
	/* init */
	size = sizeof(double);
	flag = 1;
	precision = 9;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	sep = '\t';
	method = "cos";
	tmp = 0;
	
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
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "print_precision") == 0) {
					precision = getNumArg(&Arg, &args, len + offset, "print_precision");
				} else if(cmdcmp(arg, "distance") == 0) {
					method = getArgDie(&Arg, &args, len + offset, "distance");
				} else if(cmdcmp(arg, "distance_help") == 0) {
					method = 0;
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
					Dat_init = &DatMinit;
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
						inputfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'x') {
						precision = getNumArg(&Arg, &args, len, "x");
						opt = 0;
					} else if(opt == 'd') {
						method = getArgDie(&Arg, &args, len, "d");
						opt = 0;
					} else if(opt == 'D') {
						method = 0;
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
						Dat_init = &DatMinit;
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
		inputfilename = *Arg;
		if(--args) {
			nonOptError();
		}
	}
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "Format flags output format, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "# 1:\tRelaxed Phylip\n");
		fprintf(stdout, "#\n");
		return 0;
	}
	
	/* distance method */
	if(method == 0) {
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
	}if(strcmp(method, "cos") == 0) {
		distcmp_d = &coscmp_d;
		distcmp_f = &coscmp_f;
		distcmp_b = &coscmp_b;
	} else if(strcmp(method, "chi2") == 0) {
		distcmp_d = &chi2cmp_d;
		distcmp_f = &chi2cmp_f;
		distcmp_b = &chi2cmp_b;
	} else if(strcmp(method, "bc") == 0) {
		distcmp_d = &bccmp_d;
		distcmp_f = &bccmp_f;
		distcmp_b = &bccmp_b;
	} else if(strcmp(method, "l1") == 0) {
		distcmp_d = &l1cmp_d;
		distcmp_f = &l1cmp_f;
		distcmp_b = &l1cmp_b;
	} else if(strcmp(method, "l2") == 0) {
		distcmp_d = &l2cmp_d;
		distcmp_f = &l2cmp_f;
		distcmp_b = &l2cmp_b;
	} else if(strcmp(method, "linf") == 0) {
		distcmp_d = &linfcmp_d;
		distcmp_f = &linfcmp_f;
		distcmp_b = &linfcmp_b;
	} else if(*method == 'l') {
		distcmp_d = &lncmp_d;
		distcmp_f = &lncmp_f;
		distcmp_b = &lncmp_b;
		exponent = strtod(method + 1, &errorMsg);
		if(*errorMsg != 0) {
			invaArg("\"--distance ln\"");
		}
		distcmp_d(0, &exponent, 0);
		distcmp_f(0, (float *)(&exponent), 0);
		distcmp_b(0, (unsigned char *)(&exponent), 0);
	} else if(strcmp(method, "p") == 0) {
		distcmp_d = &pearcmp_d;
		distcmp_f = &pearcmp_f;
		distcmp_b = &pearcmp_b;
	} else {
		invaArg("\"--distance\"");
	}
	
	/* tmp dir */
	if(tmp) {
		tmpF(tmp);
	}
	
	/* set precision */
	DatInit(-size, 0);
	DatMinit(-size, 0);
	
	return tsv2phy(inputfilename, outputfilename, flag, sep, precision);
}
