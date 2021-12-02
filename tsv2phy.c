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
#include "dat.h"
#include "distcmp.h"
#include "filebuff.h"
#include "pherror.h"
#include "tmp.h"
#include "tsv.h"
#include "tsv2phy.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);
#define dfprintf(fPtr, str, d) if(fprintf(fPtr, str, d) < 0) {if(errno) {ERROR();} else {fprintf(stderr, "Could not extend file.\n"); exit(1);}}

int tsv2phy(char *inputfilename, char *outputfilename, unsigned format, unsigned char sep) {
	
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
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
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
				dfprintf(outfile, "\t%f", distcmp_d(Di, *++Dj, n));
			}
		} else if(Dfptr) {
			Dfi = Dfptr[i];
			Dfj = Dfptr - 1;
			while(--j) {
				dfprintf(outfile, "\t%f", distcmp_f(Dfi, *++Dfj, n));
			}
		} else if(Dsptr) {
			Dsi = Dsptr[i];
			Dsj = Dsptr - 1;
			while(--j) {
				dfprintf(outfile, "\t%f", distcmp_s(Dsi, *++Dsj, n));
			}
		} else {
			Dbi = Dbptr[i];
			Dbj = Dbptr - 1;
			while(--j) {
				dfprintf(outfile, "\t%f", distcmp_b(Dbi, *++Dbj, n));
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-s", "Separator", "\\t");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-d", "Distance method", "cos");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-dh", "Help on option \"-d\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	//fprintf(out, "# %16s\t%-32s\t%s\n", "-sp", "Short precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tsv2phy(int argc, char *argv[]) {
	
	int size;
	unsigned args, format;
	char *arg, *inputfilename, *outputfilename, *errorMsg;
	unsigned char sep;
	double exponent;
	
	/* init */
	size = sizeof(double);
	format = 1;
	inputfilename = "--";
	outputfilename = "--";
	sep = '\t';
	
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
			} else if(strcmp(arg, "s") == 0) {
				if(argc <= ++args) {
					missArg("\"-s\"");
				} else if((sep = argv[args][0]) == 0 || argv[args][1] != 0) {
					invaArg("\"-s\"");
				}
			} else if(strcmp(arg, "d") == 0) {
				if(++args < argc) {
					arg = argv[args];
					if(strcmp(arg, "cos") == 0) {
						distcmp_d = &coscmp_d;
						distcmp_f = &coscmp_f;
						distcmp_b = &coscmp_b;
					} else if(strcmp(arg, "chi2") == 0) {
						distcmp_d = &chi2cmp_d;
						distcmp_f = &chi2cmp_f;
						distcmp_b = &chi2cmp_b;
					} else if(strcmp(arg, "bc") == 0) {
						distcmp_d = &bccmp_d;
						distcmp_f = &bccmp_f;
						distcmp_b = &bccmp_b;
					} else if(strcmp(arg, "l1") == 0) {
						distcmp_d = &l1cmp_d;
						distcmp_f = &l1cmp_f;
						distcmp_b = &l1cmp_b;
					} else if(strcmp(arg, "l2") == 0) {
						distcmp_d = &l2cmp_d;
						distcmp_f = &l2cmp_f;
						distcmp_b = &l2cmp_b;
					} else if(strcmp(arg, "linf") == 0) {
						distcmp_d = &linfcmp_d;
						distcmp_f = &linfcmp_f;
						distcmp_b = &linfcmp_b;
					} else if(*arg == 'l') {
						distcmp_d = &lncmp_d;
						distcmp_f = &lncmp_f;
						distcmp_b = &lncmp_b;
						exponent = strtod(arg + 1, &errorMsg);
						if(*errorMsg != 0) {
							invaArg("\"-d ln\"");
						}
						distcmp_d(0, &exponent, 0);
						distcmp_f(0, (float *)(&exponent), 0);
						distcmp_b(0, (unsigned char *)(&exponent), 0);
					} else if(strcmp(arg, "p") == 0) {
						distcmp_d = &pearcmp_d;
						distcmp_f = &pearcmp_f;
						distcmp_b = &pearcmp_b;
					} else {
						invaArg("\"-d\"");
					}
				} else {
					missArg("\"-d\"");
				}
			} else if(strcmp(arg, "dh") == 0) {
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
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					format = strtoul(argv[args], &arg, 10);
					if(*arg != 0) {
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
				Dat_init = &DatMinit;
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
	DatInit(-size, 0);
	DatMinit(-size, 0);
	
	return tsv2phy(inputfilename, outputfilename, format, sep);
}
