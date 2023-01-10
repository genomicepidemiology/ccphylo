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
#include <time.h>
#include "bytescale.h"
#include "cmdline.h"
#include "dnj.h"
#include "filebuff.h"
#include "hclust.h"
#include "matrix.h"
#include "nj.h"
#include "nwck.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"
#include "tmp.h"
#include "tree.h"
#include "vector.h"

void formTree(char *inputfilename, char *outputfilename, int flag, char sep, char quotes, char m, int thread_num) {
	
	int i, *N;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	Vector *sD, *Q;
	time_t t0, t1;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == 0) ? stdout : sfopen(outputfilename, "wb");
	infile = setFileBuff(1048576);
	header = setQseqs(64);
	i = 32;
	D = ltdMatrix_init(i);
	sD = vector_init(i);
	if(m == 'e') {
		Q = 0;
		N = smalloc(i * sizeof(int));
	} else {
		Q = vector_init(i);
		N = smalloc(2 * i * sizeof(int));
	}
	names = smalloc(i * sizeof(Qseqs *));
	names += i;
	++i;
	while(--i) {
		*--names = setQseqs(4);
	}
	/* set ptr according ot flag */
	if(flag) {
		if(flag & 1) {
			formLastNodePtr = &formLastBiNode;
		}
		if(flag & 2) {
			limbLengthPtr = &limbLengthNeg;
		}
	}
	
	/* set */
	openAndDetermine(infile, inputfilename);
	
	/* generate trees */
	t0 = clock();
	while((names = loadPhy(D, names, header, infile, sep, quotes)) && D->n) {
		t1 = clock();
		fprintf(stderr, "# Total time used loading matrix: %.2f s.\n", difftime(t1, t0) / 1000000);
		t0 = t1;
		if(2 < D->n) {
			/* make tree */
			if(m == 'd') {
				N = dnj_thread(D, sD, Q, N, names, thread_num);
			} else if(m == 'e') {
				N = nj_thread(D, sD, N, names, thread_num);
			} else { /* m == 'h' */
				N = hclust(D, sD, Q, N, names);
			}
		} else if(D->n == 2) {
			/* form tree */
			formLastBiNode(*names, names[1], (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(**(D->bmat))));
		}
		
		/* output tree */
		if(header->len) {
			fprintf(outfile, ">%s%s;\n", header->seq, (*names)->seq);
		} else {
			fprintf(outfile, "%s;\n", (*names)->seq);
		}
		
		t1 = clock();
		fprintf(stderr, "# Total time used Constructing tree: %.2f s.\n", difftime(t1, t0) / 1000000);
		t0 = t1;
	}
	
	/* clean */
	fclose(outfile);
	closeFileBuff(infile);
	free(N);
	Matrix_destroy(D);
	vector_destroy(sD);
	free(names);
	destroyFileBuff(infile);
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo forms tree(s) in newick format given a set of phylip distance matrices.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'q', "quotes", "Quote taxa", "\\0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'x', "print_precision", "Floating point print precision", "9");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'm', "method", "Tree construction method.", "dnj");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'M', "method_help", "Help on option \"-m\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'p', "float_precision", "Float precision on distance matrix", "False / double");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 's', "short_precision", "Short precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "byte_precision", "Byte precision on distance matrix", "False / double / 1e0");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'g', "free", "Gradually free up D", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'H', "mmap", "Allocate matrix on the disk", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tmp", "Set directory for temporary files", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 't', "threads", "Number of threads", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tree(int argc, char **argv) {
	
	const char *stdstream = "-";
	int size, len, offset, args, flag, thread_num, precision;
	char **Arg, *arg, *inputfilename, *outputfilename, *method, *tmp;
	char m, opt, sep, quotes;
	
	/* set defaults */
	size = sizeof(double);
	flag = 0;
	thread_num = 1;
	precision = 9;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	method = "dnj";
	m = 'd';
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
					inputfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "quotes") == 0) {
					quotes = getcArgDie(&Arg, &args, len + offset, "quotes");
				} else if(cmdcmp(arg, "print_precision") == 0) {
					precision = getNumArg(&Arg, &args, len + offset, "print_precision");
				} else if(cmdcmp(arg, "method") == 0) {
					method = getArgDie(&Arg, &args, len + offset, "method");
				} else if(cmdcmp(arg, "method_help") == 0) {
					method = "mh";
				} else if(cmdcmp(arg, "flag") == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(cmdcmp(arg, "flag_help") == 0) {
					flag = -1;
				} else if(cmdcmp(arg, "threads") == 0) {
					thread_num = getNumArg(&Arg, &args, len + offset, "threads");
				} else if(cmdcmp(arg, "float_precision") == 0) {
					size = sizeof(float);
				} else if(cmdcmp(arg, "short_precision") == 0) {
					size = sizeof(short unsigned);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "short_precision");
				} else if(cmdcmp(arg, "byte_precision") == 0) {
					size = sizeof(unsigned char);
					ByteScale = getdDefArg(&Arg, &args, len + offset, ByteScale, "byte_precision");
				} else if(cmdcmp(arg, "free") == 0) {
					ltdMatrixShrink = &ltdMatrix_shrink;
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
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'q') {
						quotes = getcArgDie(&Arg, &args, len, "q");
						opt = 0;
					} else if(opt == 'x') {
						precision = getNumArg(&Arg, &args, len, "x");
						opt = 0;
					} else if(opt == 'm') {
						method = getArgDie(&Arg, &args, len, "m");
						opt = 0;
					} else if(opt == 'M') {
						method = "mh";
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
					} else if(opt == 't') {
						thread_num = getNumArg(&Arg, &args, len, "t");
						opt = 0;
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
					} else if(opt == 'g') {
						ltdMatrixShrink = &ltdMatrix_shrink;
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
	setPrecisionNwck(precision);
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "# Format flags output, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "#   1:\tStrictly bifurcate the root\n");
		fprintf(stdout, "#   2:\tAllow negative branchlengths\n");
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
	
	/* set ptrs according to method */
	if(strcmp(method, "nj") == 0) {
		/* neighbor-joining */
		m = 'e';
		updateDptr = &updateD;
		minDist = &initQ;
		minDist_thread = &initQ_thread;
		initQchunkPtr = &initQchunk;
	} else if(strcmp(method, "upgma") == 0) {
		/* upgma */
		m = 'd';
		Qpair = &UPGMApair;
		minDist_thread = &minQ_thread;
		Qrow = &UPGMArow;
		nextQrow = &nextQminRow;
		Qbool = &minQbool;
		initDsDQN = &initDmin;
		pairQ = &minQ;
		updateDsDQNPtr = &updateUPGMA;
		popArrangePtr = &UPGMA_popArrange;
		qPos = &minPos;
		/*
		m = 'e';
		updateDptr = &updateD_UPGMA;
		minDist = &minD;
		minDist_thread = &minD_thread;
		minDchunkPtr = &minDchunk;
		*/
	} else if(strcmp(method, "cf") == 0) {
		/* closest first */
		m = 'h';
		initDsDQN = &initDmin;
		pairQ = &minQ;
		updateDsDQNPtr = &updateCF;
		popArrangePtr = &UPGMA_popArrange;
		/*
		m = 'e';
		updateDptr = &updateD_CF;
		minDist = &minD;
		minDist_thread = &minD_thread;
		minDchunkPtr = &minDchunk;
		*/
	} else if(strcmp(method, "ff") == 0) {
		/* furthest first */
		m = 'd';
		Qpair = &UPGMApair;
		minDist_thread = &minQ_thread;
		Qrow = &UPGMArow;
		nextQrow = &nextQminRow;
		Qbool = &minQbool;
		initDsDQN = &initDmin;
		pairQ = &minQ;
		updateDsDQNPtr = &updateFF;
		popArrangePtr = &UPGMA_popArrange;
		qPos = &minPos;
		/*
		m = 'e';
		updateDptr = &updateD_FF;
		minDist = &minD;
		minDist_thread = &minD_thread;
		minDchunkPtr = &minDchunk;
		*/
	} else if(strcmp(method, "mn") == 0) {
		/* minimum neighbors */
		m = 'e';
		updateDptr = &updateD;
		minDist = &initQ_MN;
		minDist_thread = &initQ_thread;
		initQchunkPtr = &initQ_MNchunk;
	} else if(strcmp(method, "frank") == 0) {
		/* not frank, and add to -mh */
		/* wrong assumptions, should be moved to dnj for threading and similar workflow */
		/*
		m = 'h';
		initDsDQN = &initDmin;
		pairQ = &maxQ;
		updateDsDQNPtr = &updateCF;
		popArrangePtr = &UPGMA_popArrange;
		*/
		m = 'e';
		updateDptr = &updateD_CF;
		minDist = &maxD;
		minDist_thread = &minD_thread;
		minDchunkPtr = &maxDchunk;
	} else if(strcmp(method, "hnj") == 0) {
		/* heuristic neighbor-joining */
		m = 'h';
		initDsDQN = &initHNJ;
		pairQ = &minQ;
		updateDsDQNPtr = &updateHNJ;
		popArrangePtr = &HNJ_popArrange;
	} else if(0 && strcmp(method, "hmn") == 0) {
		/* wrong, will never work */
		m = 'h';
		initDsDQN = &initHMN;
		pairQ = &maxQ;
		updateDsDQNPtr = &updateHMN;
		popArrangePtr = &HMN_popArrange;
	} else if(strcmp(method, "dnj") == 0) {
		/* dynamic neighbor-joining */
		m = 'd';
		Qpair = &minQpair;
		minDist_thread = &minQ_thread;
		Qrow = &minQrow;
		nextQrow = &nextQminRow;
		Qbool = &minQbool;
		initDsDQN = &initHNJ;
		pairQ = &minQ;
		updateDsDQNPtr = &updateDNJ;
		popArrangePtr = &DNJ_popArrange;
		qPos = &minPos;
	} else if(0 && strcmp(method, "dmn") == 0) {
		/* wrong, will never work */
		m = 'd';
		Qpair = &maxQpair;
		minDist_thread = &minQ_thread;
		Qrow = &maxQrow;
		nextQrow = &nextQmaxRow;
		Qbool = &maxQbool;
		initDsDQN = &initHMN;
		pairQ = &maxQ;
		updateDsDQNPtr = &updateDMN;
		popArrangePtr = &HMN_popArrange;
		qPos = &maxPos;
	} else if(strcmp(method, "mh") == 0) {
		fprintf(stdout, "# Tree construction methods:\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "# %-8s\t%s\n", "nj", "Neighbor-Joining");
		fprintf(stdout, "# %-8s\t%s\n", "upgma", "UPGMA");
		fprintf(stdout, "# %-8s\t%s\n", "cf", "K-means Closest First");
		fprintf(stdout, "# %-8s\t%s\n", "ff", "K-means Furthest First");
		fprintf(stdout, "# %-8s\t%s\n", "mn", "Minimum Neighbors");
		fprintf(stdout, "# %-8s\t%s\n", "hnj", "Heuristic Neighbor-Joining");
		//fprintf(stdout, "# %-8s\t%s\n", "hmn", "Heuristic Minimum Neighbors");
		fprintf(stdout, "# %-8s\t%s\n", "dnj", "Dynamic Neighbor-Joining");
		//fprintf(stdout, "# %-8s\t%s\n", "dmn", "Dynamic Minimum Neighbors");
		fprintf(stdout, "#\n");
		return 0;
	} else {
		invaArg("\"-m\"");
	}
	
	/* make tree */
	formTree(inputfilename, outputfilename, flag, sep, quotes, m, thread_num);
	
	return 0;
}
