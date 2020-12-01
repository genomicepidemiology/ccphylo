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

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trim.h"
#include "filebuff.h"
#include "fsacmp.h"
#include "matrix.h"
#include "meth.h"
#include "methparse.h"
#include "pherror.h"
#include "phy.h"
#include "seqparse.h"
#define exchange(src1, src2, tmp) tmp = src1; src1 = src2; src2 = tmp;
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

void printTrimFsa(char *filename, unsigned char *seq, int len, unsigned *includes, int flag) {
	
	const unsigned char bases[16] = "ACGTN-RYSWKMBDHV";
	int i, shifter;
	unsigned char *seqPtr;
	
	seqPtr = seq;
	shifter = 32;
	for(i = 0; i < len; ++i, ++seqPtr) {
		if((*includes >> --shifter) & 1) {
			*seqPtr = bases[*seqPtr];
		} else if(flag & 1) {
			*seqPtr = 'N';
		} else {
			*seqPtr = tolower(bases[*seqPtr]);
		}
		if(!shifter) {
			shifter = 32;
			++includes;
		}
	}
	*seqPtr = '\n';
	
	fprintf(stdout, ">%s\n", stripDir(filename));
	fwrite(seq, 1, len + 1, stdout);
}

void fsaTrim(int numFile, char *targetTemplate, char **filenames, unsigned minLength, double minCov, unsigned flag, unsigned proxi, char *methfilename) {
	
	int i, pair, inludeN, len, cSize;
	unsigned *includes;
	long unsigned *nibbleSeq;
	unsigned char *trans, **seqs, **seqsPtr;
	FileBuff *infile;
	MethMotif *motif;
	Qseqs *header, *ref, *seq;
	
	/* init */
	includes = 0;
	--filenames;
	pair = flag & 2 ? 1 : 0;
	inludeN = numFile;
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	ref = 0;
	seq = setQseqs(1048576);
	if(!(seqs = calloc(numFile, sizeof(unsigned char *)))) {
		ERROR();
	}
	seqsPtr = seqs;
	trans = getIupacBitTable(flag);
	nibbleSeq = 0;
	len = 0;
	if(methfilename) {
		openAndDetermine(infile, methfilename);
		motif = getMethMotifs(infile, seq);
		closeFileBuff(infile);
	} else {
		motif = 0;
	}
	
	/* load sequences and trim positions */
	i = numFile + 1;
	while(--i) {
		/* open and validate file */
		if(openAndDetermine(infile, *++filenames) != '>') {
			fprintf(stderr, "\"%s\" is not fasta.\n", *filenames);
			exit(1);
		}
		
		/* find the correct entry */
		while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
		
		/* get sequence */
		if(!header->len) {
			fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
			--inludeN;
		} else if(FileBuffgetFsaSeq(infile, seq, trans)) {
			if(ref) {
				if(seq->len != ref->len) {
					fprintf(stderr, "Sequences does not match: %s\n", *filenames);
					exit(1);
				}
				
				/* make / update inclusion array(s) */
				if(pair) {
					initIncPos(includes, len);
					qseq2nibble(seq, nibbleSeq);
					maskMotifs(nibbleSeq, includes, len, motif);
					getIncPosPtr(includes, seq, seq, proxi);
					if(getNpos(includes, len) < minLength) {
						fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filenames);
						--inludeN;
					}
				} else {
					qseq2nibble(seq, nibbleSeq);
					maskMotifs(nibbleSeq, includes, len, motif);
					getIncPosPtr(includes, seq, ref, proxi);
					*seqsPtr = smalloc(len + 1);
					memcpy(*seqsPtr, seq->seq, len + 1);
				}
			} else {
				len = seq->len;
				minLength = minLength < minCov * len ? minCov * len : minLength;
				cSize = len / 32 + 1;
				if(!(nibbleSeq = calloc(cSize, sizeof(long unsigned)))) {
					ERROR();
				}
				includes = smalloc(cSize * sizeof(unsigned));
				
				initIncPos(includes, len);
				qseq2nibble(seq, nibbleSeq);
				maskMotifs(nibbleSeq, includes, len, motif);
				getIncPos(includes, seq, seq, proxi);
				
				if(getNpos(includes, len) < minLength) {
					free(nibbleSeq);
					free(includes);
					includes = 0;
					--inludeN;
				} else if(!pair) {
					*seqsPtr = smalloc(len + 1);
					memcpy(*seqsPtr, seq->seq, len + 1);
					ref = seq;
					seq = setQseqs(len + 1);
				}
			}
		} else {
			fprintf(stderr, "Missing template sequence (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
			--inludeN;
		}
		
		if(pair) {
			printTrimFsa(*filenames, seq->seq, len, includes, flag);
		} else {
			++seqsPtr;
		}
		closeFileBuff(infile);
	}
	
	/* make final output */
	if(!inludeN) {
		errno = 1;
		fprintf(stderr, "All sequences were trimmed away.\n");
	} if(!pair) {
		if(minLength <= getNpos(includes, len)) {
			/* print all sequences */
			i = numFile + 1;
			while(--i) {
				if(*--seqsPtr) {
					printTrimFsa(*filenames, *seqsPtr, len, includes, flag);
				}
				--filenames;
			}
		} else {
			errno = 1;
			fprintf(stderr, "No sufficient overlap was found.\n");
		}
	}
	
	/* clean */
	destroyFileBuff(infile);
	destroyQseqs(header);
	destroyQseqs(seq);
	if(ref) {
		destroyQseqs(ref);
	}
	i = numFile + 1;
	seqsPtr = seqs - 1;
	while(--i) {
		if(*++seqsPtr) {
			free(*seqsPtr);
		}
	}
	free(seqs);
	if(nibbleSeq) {
		free(nibbleSeq);
	}
	destroyMethMotifs(motif);
	
}


static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo trims multiple alignments from different files, and merge them into one\n");
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file(s)", "stdin");
	//fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-m", "Mask methylation motifs from <file>", "False/None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-r", "Target reference", "None");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mc", "Minimum coverage", "50.0%");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-ml", "Minimum overlapping length", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-pr", "Minimum proximity between SNPs", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_trim(int argc, char *argv[]) {
	
	unsigned args, numFile, flag, minLength, proxi;
	char *arg, *targetTemplate, **filenames, *errorMsg, *methfilename;
	double minCov;
	
	/* set defaults */
	numFile = 0;
	flag = 0;
	minLength = 1;
	proxi = 0;
	targetTemplate = 0;
	filenames = 0;
	methfilename = 0;
	minCov = 0.5;
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					filenames = argv + args;
					numFile = 1;
					while(++args < argc && (*argv[args] != '-' || (argv[args][1] == '-' && argv[args][2] == 0))) {
						++numFile;
					}
					if(numFile == 0) {
						missArg("\"-i\"");
					}
					--args;
				}
			} else if(strcmp(arg, "m") == 0) {
				if(++args < argc) {
					methfilename = argv[args];
				} else {
					missArg("\"-m\"");
				}
			} else if(strcmp(arg, "r") == 0) {
				if(++args < argc) {
					targetTemplate = argv[args];
				} else {
					missArg("\"-r\"");
				}
			} else if(strcmp(arg, "ml") == 0) {
				if(++args < argc) {
					if((minLength = strtoul(argv[args], &errorMsg, 10)) == 0) {
						minLength = 1;
					}
					if(*errorMsg != 0) {
						invaArg("\"-ml\"");
					}
				} else {
					missArg("\"-ml\"");
				}
			} else if(strcmp(arg, "mc") == 0) {
				if(++args < argc) {
					minCov = strtod(argv[args], &errorMsg);
					if(*errorMsg != 0 || minCov < 0 || 100 < minCov) {
						invaArg("\"-mc\"");
					}
					minCov /= 100.0;
				} else {
					missArg("\"-mc\"");
				}
			} else if(strcmp(arg, "pr") == 0) {
				if(++args < argc) {
					proxi = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-pr\"");
					}
				} else {
					missArg("\"-pr\"");
				}
			} else if(strcmp(arg, "f") == 0) {
				if(++args < argc) {
					flag = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-f\"");
				}
			} else if(strcmp(arg, "fh") == 0) {
				fprintf(stdout, "# Format flags output, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#   1:\tHard mask\n");
				fprintf(stdout, "#   2:\tPairwise comparison\n");
				fprintf(stdout, "#   8:\tInclude insignificant bases in distance calculation, only affects fasta input\n");
				fprintf(stdout, "#  32:\tDo not include insignificant bases in pruning\n");
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
	
	/* set pointers */
	if(flag & 32) {
		getIncPosPtr = &getIncPosInsigPrune;
	} else if(flag & 8) {
		getIncPosPtr = &getIncPosInsig;
	}
	
	/* check for required input */
	if(!numFile && targetTemplate) {
		numFile = 1;
	}
	if(!targetTemplate) {
		fprintf(stderr, "Target template is required.\n");
		exit(1);
	}
	fsaTrim(numFile, targetTemplate, filenames, minLength, minCov, flag, proxi, methfilename);
	
	return 0;
}
