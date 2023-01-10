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
#include "bytescale.h"
#include "cmdline.h"
#include "filebuff.h"
#include "matrix.h"
#include "nwck.h"
#include "nwck2phy.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"

int newick2phy(char *inputfilename, char *outputfilename, unsigned format) {
	
	int i, j, n, org_i, namesize;
	unsigned char **sNames;
	double Li, Lj, d, **Dmat, *D_i, *ptr;
	float **Dfmat, *Df_i, *fptr;
	short unsigned **Dsmat, *Ds_i, *sptr;
	unsigned char **Dbmat, *Db_i, *bptr;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs *name, **names, *header;
	
	/* init */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	D = ltdMatrix_init(128);
	name = setQseqs(1024);
	header = setQseqs(64);
	names = smalloc(D->size * sizeof(Qseqs *));
	i = D->size;
	while(i--) {
		names[i] = smalloc(sizeof(Qseqs));
	}
	namesize = D->size;
	sNames = smalloc(D->size * sizeof(char *));
	
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
		Dmat = 0;
		Dfmat = 0;
		Dsmat = 0;
		Dbmat = 0;
		if(D->mat) {
			Dmat = D->mat;
		} else if(D->fmat) {
			Dfmat = D->fmat;
		} else if(D->smat) {
			Dsmat = D->smat;
		} else {
			Dbmat = D->bmat;
		}
		
		/* get distances */
		(*names)->seq = name->seq;
		(*names)->len = name->len;
		(*names)->size = name->size;
		org_i = 0;
		D->n = 1;
		if(Dmat) {
			while(D->n != n) {
				/* split */
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
		} else if(Dfmat) {
			while(D->n != n) {
				/* slpit */
				if(splitNwck(names[org_i], names[D->n], &Li, &Lj)) {
					/* 
					update last / new row,
					new dist is: dist to previous node + limblength
					*/
					fptr = Dfmat[D->n] - 1;
					if(Lj < 0) {
						j = D->n + 2;
						while(--j) {
							*++fptr = Lj;
						}
					} else {
						i = org_i;
						Df_i = Dfmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++fptr = *++Df_i < 0 ? -1.0 : Lj + *Df_i;
						}
						*++fptr = Lj + Li;
						for(j = i + 1; j <= D->n; ++j) {
							*++fptr = (d = Dfmat[j][i]) < 0 ? -1.0 : Lj + d;
						}
					}
					
					/* update originating row */
					if(Li < 0) {
						i = org_i;
						Df_i = Dfmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++Df_i = Li;
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							Dfmat[i][org_i] = Li;
						}
					} else {
						i = org_i;
						Df_i = Dfmat[i] - 1;
						j = i + 1;
						while(--j) {
							if(0 <= *++Df_i) {
								*Df_i += Li;
							}
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							if(0 <= *(Df_i = Dfmat[i] + org_i)) {
								*Df_i += Li;
							}
						}
					}
					
					/* update matrix size */
					++D->n;
				} else {
					++org_i;
				}
			}
		} else if(Dsmat) {
			while(D->n != n) {
				/* slpit */
				if(splitNwck(names[org_i], names[D->n], &Li, &Lj)) {
					/* 
					update last / new row,
					new dist is: dist to previous node + limblength
					*/
					sptr = Dsmat[D->n] - 1;
					if(Lj < 0) {
						j = D->n + 2;
						while(--j) {
							*++sptr = dtouc(Lj, 0);
						}
					} else {
						i = org_i;
						Ds_i = Dsmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++sptr = uctod(*++Ds_i) < 0 ? dtouc(-1.0, 0) : dtouc(Lj, 0) + *Ds_i;
						}
						*++sptr = dtouc(Lj + Li, 0);
						for(j = i + 1; j <= D->n; ++j) {
							d = uctod(Dsmat[j][i]);
							*++sptr = d < 0 ? dtouc(-1.0, 0) : dtouc((Lj + d), 0);
						}
					}
					
					/* update originating row */
					if(Li < 0) {
						i = org_i;
						Ds_i = Dsmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++Ds_i = dtouc(Li, 0);
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							Dsmat[i][org_i] = dtouc(Li, 0);
						}
					} else {
						i = org_i;
						Ds_i = Dsmat[i] - 1;
						j = i + 1;
						while(--j) {
							if(0 <= uctod(*++Ds_i)) {
								*Ds_i += dtouc(Li, 0);
							}
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							if(0 <= uctod(*(Ds_i = Dsmat[i] + org_i))) {
								*Ds_i += dtouc(Li, 0);
							}
						}
					}
					
					/* update matrix size */
					++D->n;
				} else {
					++org_i;
				}
			}
		} else {
			while(D->n != n) {
				/* slpit */
				if(splitNwck(names[org_i], names[D->n], &Li, &Lj)) {
					/* 
					update last / new row,
					new dist is: dist to previous node + limblength
					*/
					bptr = Dbmat[D->n] - 1;
					if(Lj < 0) {
						j = D->n + 2;
						while(--j) {
							*++bptr = dtouc(Lj, 0);
						}
					} else {
						i = org_i;
						Db_i = Dbmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++bptr = uctod(*++Db_i) < 0 ? dtouc(-1.0, 0) : dtouc(Lj, 0) + *Db_i;
						}
						*++bptr = dtouc(Lj + Li, 0);
						for(j = i + 1; j <= D->n; ++j) {
							d = uctod(Dbmat[j][i]);
							*++bptr = d < 0 ? dtouc(-1.0, 0) : dtouc((Lj + d), 0);
						}
					}
					
					/* update originating row */
					if(Li < 0) {
						i = org_i;
						Db_i = Dbmat[i] - 1;
						j = i + 1;
						while(--j) {
							*++Db_i = dtouc(Li, 0);
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							Dbmat[i][org_i] = dtouc(Li, 0);
						}
					} else {
						i = org_i;
						Db_i = Dbmat[i] - 1;
						j = i + 1;
						while(--j) {
							if(0 <= uctod(*++Db_i)) {
								*Db_i += dtouc(Li, 0);
							}
						}
						
						/* update originating column */
						for(i = org_i + 1; i < D->n; ++i) {
							if(0 <= uctod(*(Db_i = Dbmat[i] + org_i))) {
								*Db_i += dtouc(Li, 0);
							}
						}
					}
					
					/* update matrix size */
					++D->n;
				} else {
					++org_i;
				}
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
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'x', "print_precision", "Floating point print precision", "9");
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

int main_nwck2phy(int argc, char **argv) {
	
	const char *stdstream = "-";
	int size, len, offset, args, flag, precision;
	char **Arg, *arg, *inputfilename, *outputfilename, *tmp, opt;
	
	/* init */
	size = sizeof(double);
	flag = 1;
	precision = 9;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	stripEntry = &noStripDir;
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
				} else if(cmdcmp(arg, "print_precision") == 0) {
					precision = getNumArg(&Arg, &args, len + offset, "print_precision");
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
						inputfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'x') {
						precision = getNumArg(&Arg, &args, len, "x");
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
		inputfilename = *Arg;
		if(--args) {
			nonOptError();
		}
	}
	
	/* set print precision */
	setPrecisionPhy(precision);
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "Format flags output format, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "# 1:\tRelaxed Phylip\n");
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
	
	newick2phy(inputfilename, outputfilename, flag);
	
	return 0;
}
