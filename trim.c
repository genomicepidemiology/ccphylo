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
#include "cmdline.h"
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

void printTrimFsa(FILE *out, char *filename, unsigned char *seq, int len, unsigned *includes, int flag) {
	
	const unsigned char bases[16] = "ACGTN-RYSWKMBDHV";
	int i, shifter;
	unsigned char *seqPtr;
	
	
	fprintf(out, ">%s\n", stripDir(filename));
	seqPtr = seq;
	shifter = 32;
	if((flag & 18) == 16) { /* not 2 && 16 */
		for(i = 0; i < len; ++i, ++seqPtr) {
			if((*includes >> --shifter) & 1) {
				putc(bases[*seqPtr], out);
			}
			if(!shifter) {
				shifter = 32;
				++includes;
			}
		}
		putc('\n', out);
	} else {
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
		cfwrite(seq, 1, len + 1, out);
	}
}

void fsaTrim(int numFile, char *targetTemplate, char **filenames, char *outputfilename, unsigned minLength, double minCov, unsigned flag, unsigned proxi, char *methfilename) {
	
	int i, pair, inludeN, len, cSize, numSeqs, maxSeqs;
	unsigned *includes;
	long unsigned *nibbleSeq;
	unsigned char *trans, **seqs, **seqsPtr, **seqnames;
	FILE *out;
	FileBuff *infile;
	MethMotif *motif;
	Qseqs *header, *ref, *seq;
	
	/* here */
	/*
	add "seqnames" as "filenames" when no reference entry is given
	loop over files when no reference is added
	*/
	
	/* init */
	numSeqs = 0;
	includes = 0;
	--filenames;
	pair = flag & 2 ? 1 : 0;
	inludeN = 0;
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	ref = 0;
	seq = setQseqs(1048576);
	if(!(seqs = calloc(numFile, sizeof(unsigned char *)))) {
		ERROR();
	}
	seqsPtr = seqs;
	trans = getIupacBitTable(flag);
	maxSeqs = numFile;
	if(!pair && !targetTemplate) {
		seqnames = smalloc(maxSeqs * sizeof(unsigned char *));
	} else {
		seqnames = 0;
	}
	nibbleSeq = 0;
	len = 0;
	if(methfilename) {
		openAndDetermine(infile, methfilename);
		motif = getMethMotifs(infile, seq);
		closeFileBuff(infile);
	} else {
		motif = 0;
	}
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		out = stdout;
	} else {
		out = sfopen(outputfilename, "wb");
	}
	
	/* load sequences and trim positions */
	i = numFile + 1;
	while(--i) {
		/* open and validate file */
		if(openAndDetermine(infile, *++filenames) != '>') {
			fprintf(stderr, "\"%s\" is not fasta.\n", *filenames);
			exit(1);
		}
		
		/* go over entries, when no reference has been parsed */
		do {
			if(numSeqs == maxSeqs) {
				if(!(seqs = realloc(seqs, (maxSeqs <<= 1) * sizeof(unsigned char *)))) {
					ERROR();
				} else if(!targetTemplate && !(seqnames = realloc(seqnames, maxSeqs * sizeof(unsigned char *)))) {
					ERROR();
				} else {
					seqsPtr = seqs + numSeqs;
					memset(seqsPtr, 0, numSeqs * sizeof(unsigned char *));
				}
			}
			
			/* find the correct entry */
			while(FileBuffgetFsaHeader(infile, header) && (targetTemplate && strcmp((char *) header->seq, targetTemplate)));
			
			/* get sequence */
			if(FileBuffgetFsaSeq(infile, seq, trans)) {
				if(ref) {
					if(seq->len != ref->len) {
						fprintf(stderr, "Sequences does not match: %s %s\n", header->seq, *filenames);
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
						} else {
							++inludeN;
						}
					} else {
						qseq2nibble(seq, nibbleSeq);
						maskMotifs(nibbleSeq, includes, len, motif);
						getIncPosPtr(includes, seq, ref, proxi);
						*seqsPtr = smalloc(len + 1);
						memcpy(*seqsPtr, seq->seq, len + 1);
						++numSeqs;
						++inludeN;
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
					} else if(!pair) {
						*seqsPtr = smalloc(len + 1);
						memcpy(*seqsPtr, seq->seq, len + 1);
						ref = seq;
						seq = setQseqs(len + 1);
						seq->len = ref->len;
						++numSeqs;
					}
					++inludeN;
				}
				if(!pair && !targetTemplate) {
					seqnames[numSeqs - 1] = smalloc(header->len + 1);
					memcpy(seqnames[numSeqs - 1], header->seq, header->len + 1);
				}
				if(pair) {
					printTrimFsa(out, targetTemplate ? *filenames : (char *)(header->seq), seq->seq, len, includes, flag);
				} else {
					++seqsPtr;
				}
			} else if(targetTemplate && !pair) {
				++seqsPtr;
			}
		} while(targetTemplate == 0 && header->len);
		
		if(targetTemplate && (!header->len || !seq->len)) {
			fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
		}
		
		closeFileBuff(infile);
	}
	
	/* make final output */
	if(!inludeN) {
		errno = 1;
		fprintf(stderr, "All sequences were trimmed away.\n");
	} if(!pair) {
		if(minLength <= getNpos(includes, len)) {
			/* create pseudo inclusion criteria */
			if(flag & 16) {
				/* here */
				pseudoAlnPrune(includes, seqs, len, (targetTemplate ? numFile : numSeqs));
			}
			/* print all sequences */
			i = (targetTemplate ? numFile : numSeqs) + 1;
			while(--i) {
				if(*--seqsPtr) {
					printTrimFsa(out, targetTemplate ? *filenames : (char *)(seqnames[i-1]), *seqsPtr, len, includes, flag);
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
	i = numSeqs + 1;
	seqsPtr = seqs - 1;
	while(--i) {
		if(*++seqsPtr) {
			free(*seqsPtr);
		}
	}
	free(seqs);
	if(seqnames) {
		i = numSeqs + 1;
		seqsPtr = seqnames - 1;
		while(--i) {
			if(*++seqsPtr) {
				free(*seqsPtr);
			}
		}
		free(seqnames);
	}
	if(nibbleSeq) {
		free(nibbleSeq);
	}
	destroyMethMotifs(motif);
	fclose(out);
}

void fsaTrimOld(int numFile, char *targetTemplate, char **filenames, char *outputfilename, unsigned minLength, double minCov, unsigned flag, unsigned proxi, char *methfilename) {
	
	int i, pair, inludeN, len, cSize;
	unsigned *includes;
	long unsigned *nibbleSeq;
	unsigned char *trans, **seqs, **seqsPtr;
	FILE *out;
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
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		out = stdout;
	} else {
		out = sfopen(outputfilename, "wb");
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
		while(FileBuffgetFsaHeader(infile, header) && (targetTemplate && strcmp((char *) header->seq, targetTemplate)));
		
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
			printTrimFsa(out, *filenames, seq->seq, len, includes, flag);
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
					printTrimFsa(out, *filenames, *seqsPtr, len, includes, flag);
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
	fclose(out);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo trims multiple alignments from different files, and merge them into one\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file(s)", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'y', "methylation_motifs", "Mask methylation motifs from <file>", "False/None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'r', "reference", "Target reference identifier", "None");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'C', "min_cov", "Minimum coverage", "50.0%");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'L', "min_len", "Minimum overlapping length", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'P', "proximity", "Minimum proximity between SNPs", "0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
	
	/*
	i	i	input
	o	o	output
	m	y	methylation_motifs
	r	r	reference
	mc	C	min_cov
	ml	L	min_len
	pr	P	proximity
	f	f	flag
	fh	F	flag_help
	h	h	help
	*/
}

int main_trim(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, len, offset, numFile, flag, minLength, proxi;
	char **Arg, *arg, *targetTemplate, **filenames, *methfilename;
	char *outputfilename, opt;
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
					filenames = getArgListDie(&Arg, &args, len + offset, "input");
					numFile = getArgListLen(&Arg, &args);
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "methylation_motifs") == 0) {
					methfilename = getArgDie(&Arg, &args, len + offset, "methylation_motifs");
				} else if(cmdcmp(arg, "reference") == 0) {
					targetTemplate = getArgDie(&Arg, &args, len + offset, "reference");
				} else if(cmdcmp(arg, "min_cov") == 0) {
					minCov = getdArg(&Arg, &args, len + offset, "min_cov") / 100;
				} else if(cmdcmp(arg, "min_len") == 0) {
					minLength = getNumArg(&Arg, &args, len + offset, "min_len");
				} else if(cmdcmp(arg, "proximity") == 0) {
					proxi = getNumArg(&Arg, &args, len + offset, "proximity");
				} else if(cmdcmp(arg, "flag") == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(cmdcmp(arg, "flag_help") == 0) {
					flag = -1;
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
					} else if(opt == 'y') {
						methfilename = getArgDie(&Arg, &args, len, "y");
						opt = 0;
					} else if(opt == 'r') {
						targetTemplate = getArgDie(&Arg, &args, len, "r");
						opt = 0;
					} else if(opt == 'C') {
						minCov = getdArg(&Arg, &args, len, "C") / 100;
						opt = 0;
					} else if(opt == 'L') {
						minLength = getNumArg(&Arg, &args, len, "L");
						opt = 0;
					} else if(opt == 'P') {
						proxi = getNumArg(&Arg, &args, len, "P");
						opt = 0;
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
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
	
	/* set pointers */
	if(flag == -1) {
		fprintf(stdout, "# Format flags output, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "#   1:\tHard mask\n");
		fprintf(stdout, "#   2:\tPairwise comparison\n");
		fprintf(stdout, "#   8:\tInclude insignificant bases in distance calculation, only affects fasta input\n");
		fprintf(stdout, "#  16:\tCreate pseudo alignment, not compatible with pairwise comparison\n");
		fprintf(stdout, "#  32:\tDo not include insignificant bases in pruning\n");
		fprintf(stdout, "#\n");
		return 0;
	} else if(flag & 32) {
		getIncPosPtr = &getIncPosInsigPrune;
	} else if(flag & 8) {
		getIncPosPtr = &getIncPosInsig;
	}
	
	/* check for required input */
	if(!numFile && targetTemplate) {
		numFile = 1;
	}
	
	fsaTrim(numFile, targetTemplate, filenames, outputfilename, minLength, minCov, flag, proxi, methfilename);
	
	return 0;
}
