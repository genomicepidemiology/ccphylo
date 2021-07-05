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
#define missArg(opt) fprintf(stderr, "Missing argument at %s.", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

void formTree(char *inputfilename, char *outputfilename, int flag, char m, int thread_num) {
	
	int i, *N;
	FILE *outfile;
	FileBuff *infile;
	Matrix *D;
	Qseqs **names, *header;
	Vector *sD, *Q;
	time_t t0, t1;
	
	/* init */
	outfile = (*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) ? stdout : sfopen(outputfilename, "wb");
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
	while((names = loadPhy(D, names, header, infile)) && D->n) {
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
			
			/* output tree */
			if(header->len) {
				fprintf(outfile, ">%s%s;\n", header->seq, (*names)->seq);
			} else {
				fprintf(outfile, "%s;\n", (*names)->seq);
			}
		} else if(D->n == 2) {
			/* output tree */
			if(header->len) {
				fprintf(outfile, ">%s(%s,%s:%.2f);\n", header->seq, (*names)->seq, names[1]->seq, **(D->mat));
			} else {
				fprintf(outfile, "(%s,%s:%.2f);\n", (*names)->seq, names[1]->seq, **(D->mat));
			}
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
	fprintf(out, "# %16s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-i", "Input file.", "stdin");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-o", "Output file", "stdout");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-f", "Output flags", "0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fh", "Help on option \"-f\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-m", "Tree construction method.", "dnj");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mh", "Help on option \"-m\"", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-fp", "Float precision on distance matrix", "double");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-bp", "Byte precision on distance matrix", "double / 1e0");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-gf", "Gradually free up D", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-mm", "Allocate matrix on the disk", "False");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-tmp", "Set directory for temporary files", "");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-t", "Number of threads", "1");
	fprintf(out, "# %16s\t%-32s\t%s\n", "-h", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_tree(int argc, char *argv[]) {
	
	int size;
	unsigned args, flag, thread_num;
	char *arg, *inputfilename, *outputfilename, *method, *errorMsg, m;
	
	/* set defaults */
	size = sizeof(double);
	flag = 0;
	thread_num = 1;
	inputfilename = "--";
	outputfilename = "--";
	method = "dnj";
	m = 'd';
	
	args = 1;
	while(args < argc) {
		arg = argv[args];
		if(*arg++ == '-') {
			if(strcmp(arg, "i") == 0) {
				if(++args < argc) {
					inputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "o") == 0) {
				if(++args < argc) {
					outputfilename = argv[args];
				} else {
					missArg("\"-o\"");
				}
			} else if(strcmp(arg, "m") == 0) {
				if(++args < argc) {
					method = argv[args];
				} else {
					missArg("\"-m\"");
				}
			} else if(strcmp(arg, "mh") == 0) {
				method = "mh";
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
				fprintf(stdout, "#   1:\tStrictly bifurcate the root\n");
				fprintf(stdout, "#   2:\tAllow negative branchlengths\n");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(arg, "t") == 0) {
				if(++args < argc) {
					thread_num = strtoul(argv[args], &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-t\"");
					}
				} else {
					missArg("\"-t\"");
				}
			} else if(strcmp(arg, "fp") == 0) {
				size = sizeof(float);
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
			} else if(strcmp(arg, "gf") == 0) {
				ltdMatrixShrink = &ltdMatrix_shrink;
			} else if(strcmp(arg, "mm") == 0) {
				ltdMatrix_init = &ltdMatrixMinit;
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
		/* minimum neighbours */
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
	formTree(inputfilename, outputfilename, flag, m, thread_num);
	
	return 0;
}
