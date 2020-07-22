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

#include <pthread.h>
#include "matrix.h"
#include "qseqs.h"
#include "vector.h"

#ifndef NJ
typedef struct njThread NJthread;
struct njThread {
	pthread_t id;
	Matrix *D;
	Vector *sD;
	unsigned *N;
	double min;
	int mi;
	int mj;
	unsigned num;
	struct njThread *next;
};
#define NJ 1
#endif

extern void (*updateDptr)(Matrix *, Vector *, unsigned *, unsigned, unsigned, double, double);
extern long unsigned (*minDist)(Matrix *, Vector *, unsigned *);
extern void * (*minDist_thread)(void *);

void limbLength(double *Li, double *Lj, unsigned i, unsigned j, Vector *sD, unsigned *N, double D_ij);
unsigned * initSummaD(Vector *sD, Matrix *D, unsigned *N);
long unsigned initQ(Matrix *D, Vector *sD, unsigned *N);
double initQchunk(Matrix *D, double *sD, unsigned *N, double min, int *mi, int *mj, int i, int j);
void * initQ_thread(void *arg);
long unsigned minD(Matrix *D, Vector *sD, unsigned *N);
double minDchunk(Matrix *D, unsigned *N, double min, int *mi, int *mj, int i, int j);
void * minD_thread(void *arg);
long unsigned minPair(Matrix *Q);
void updateD(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj);
void updateD_UPGMA(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj);
void updateD_FF(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj);
void updateD_CF(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj);
unsigned * nj(Matrix *D, Vector *sD, unsigned *N, Qseqs **names);
unsigned * nj_thread(Matrix *D, Vector *sD, unsigned *N, Qseqs **names, int thread_num);
