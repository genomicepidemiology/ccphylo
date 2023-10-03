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
#include "matrix.h"
#include "qseqs.h"
#define destroyTable(trans)(free(trans - 128))

extern void (*getIncPosPtr)(unsigned*, Qseqs*, Qseqs*, unsigned);

unsigned char * get2BitTable(unsigned flag);
unsigned char * getIupacBitTable(unsigned flag);
void initIncPos(unsigned *include, int len);
void getIncPos(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi);
void getIncPosInsigPrune(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi);
void getIncPosInsig(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi);
void maskProxi(unsigned *include, unsigned *include1, unsigned *include2, long unsigned *seq1, long unsigned *seq2, unsigned len, unsigned proxi);
int getNpos(unsigned *include, int len);
void pseudoAlnPrune(unsigned *include, unsigned char **seqs, int len, int n);
unsigned fsacmp(long unsigned *seq1, long unsigned *seq2, unsigned *include, int len);
long unsigned fsacmpair(long unsigned *seq1, long unsigned *seq2, unsigned *include, int len);
void printDiff(FILE *outfile, int samplei, int samplej, int nuci, int pos, int nucj);
unsigned fsacmprint(FILE *outfile, int samplei, int samplej, long unsigned *seq1, long unsigned *seq2, unsigned *include, int len);
long unsigned fsacmpairint(FILE *outfile, int samplei, int samplej, long unsigned *seq1, long unsigned *seq2, unsigned *include, int len);
unsigned cmpFsa(Matrix *D, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, FILE *diffile);
void cmpairFsa(Matrix *D, Matrix *N, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, unsigned minLength, double minCov, unsigned proxi, FILE *diffile);
