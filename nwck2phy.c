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
#include "filebuff.h"
#include "matrix.h"
#include "nwck.h"
#include "nwck2phy.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

int newick2phy(char *inputfilename, char *outputfilename, unsigned format) {
	
	int i, j, n, org_i, namesize;
	unsigned char **sNames;
	double Li, Lj, d, **Dmat, *D_i, *ptr;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs *name, **names, *header;
	
	/* init */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "rb");
	}
	D = ltdMatrix_init(128);
	name = setQseqs(1024);
	header = setQseqs(64);
	names = smalloc(128 * sizeof(Qseqs *));
	i = 128;
	while(i--) {
		names[i] = setQseqs(32);
	}
	namesize = 128;
	sNames = smalloc(128 * sizeof(char *));
	
	Dmat = D->mat;
	while(getNwck(infile, name, header)) {
		n = getSizeNwck(name);
		
		/* alloc */
		if(D->size < n) {
			if(!(names = realloc(names, n * sizeof(Qseqs *)))) {
				ERROR();
			}
			for(i = D->size; i < n; ++i) {
				names[i] = smalloc(sizeof(Qseqs));
			}
			ltdMatrix_realloc(D, n);
		}
		
		/* get distances */
		(*names)->seq = name->seq;
		(*names)->len = name->len;
		(*names)->size = name->size;
		org_i = 0;
		D->n = 1;
		while(D->n != n) {
			/* slpit */
			if(splitNwck(names[org_i], names[D->n], &Li, &Lj)) {
				/* 
				update last / new row,
				new dist is: dist to previous node + limblength
				*/
				ptr = Dmat[D->n] - 1;
				if(Lj < 0) {
					j = D->n + 2;
					while(--j) {
						*++ptr = Lj;
					}
				} else {
					i = org_i;
					D_i = Dmat[i] - 1;
					j = i + 1;
					while(--j) {
						*++ptr = *++D_i < 0 ? -1.0 : Lj + *D_i;
					}
					*++ptr = Lj + Li;
					for(j = i + 1; j <= D->n; ++j) {
						*++ptr = (d = Dmat[j][i]) < 0 ? -1.0 : Lj + d;
					}
				}
				
				/* update originating row */
				if(Li < 0) {
					i = org_i;
					D_i = Dmat[i] - 1;
					j = i + 1;
					while(--j) {
						*++D_i = Li;
					}
					
					/* update originating column */
					for(i = org_i + 1; i < D->n; ++i) {
						Dmat[i][org_i] = Li;
					}
				} else {
					i = org_i;
					D_i = Dmat[i] - 1;
					j = i + 1;
					while(--j) {
						if(0 <= *++D_i) {
							*D_i += Li;
						}
					}
					
					/* update originating column */
					for(i = org_i + 1; i < D->n; ++i) {
						if(0 <= *(D_i = Dmat[i] + org_i)) {
							*D_i += Li;
						}
					}
				}
				
				/* update matrix size */
				++D->n;
			} else {
				++org_i;
			}
		}
		
		/* convert Qseqs* to char* */
		if(namesize < D->n) {
			namesize = D->n;
			free(sNames);
			sNames = smalloc(namesize * sizeof(char *));
		}
		i = D->n;
		while(i--) {
			sNames[i] = names[i]->seq;
		}
		
		/* print new phylip matrix */
		printphy(outfile, D, (char **) sNames, 0, (char *) header->seq, format);
		
	}
	
	closeFileBuff(infile);
	destroyFileBuff(infile);
	destroyQseqs(name);
	free(names);
	free(sNames);
	
	return 0;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo nwck2phy converts newick files to phylip distance files.\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_nwck2phy(int argc, char *argv[]) {
	
	unsigned args, format;
	char *arg, *inputfilename, *outputfilename;
	
	format = 1;
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
	
	newick2phy(inputfilename, outputfilename, format);
	
	return 0;
}
