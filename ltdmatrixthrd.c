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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fbseek.h"
#include "filebuff.h"
#include "ltdmatrixthrd.h"
#include "matcmp.h"
#include "matparse.h"
#include "matrix.h"
#include "pherror.h"
#include "threader.h"
#define strcmp2(str1, str2) (*(str1) == *(str2) && strcmp(str1 + 1, str2 + 1) == 0)

static CmpMatArg * formThread(CmpMatArg *thread, Matrix *D, Matrix *N, MatrixCounts *mat1, Qseqs **filenames, unsigned n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	CmpMatArg *node;
	
	/* get input */
	node = smalloc(sizeof(CmpMatArg));
	node->D = D;
	node->N = N;
	node->mat1 = mat1;
	node->filenames = filenames;
	node->n = n;
	node->norm = norm;
	node->minDepth = minDepth;
	node->minLength = minLength;
	node->minCov = minCov;
	node->veccmp = veccmp;
	
	/* link to remaining threads */
	node->next = thread;
	
	return node;
}

static void joinThreads(CmpMatArg *src) {
	
	CmpMatArg *src_next;
	
	while(src) {
		src_next = src->next;
		if((errno = pthread_join(src->id, NULL))) {
			ERROR();
		}
		free(src);
		src = src_next;
	}
}

void matCmpThreadOut(int tnum, void * (*func)(void*), Matrix *D, Matrix *N, MatrixCounts *mat1, Qseqs **filenames, unsigned n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	int i;
	CmpMatArg *thread;
	
	/* thread out */
	i = tnum;
	thread = 0;
	do {
		/* make thread */
		thread = formThread(thread, D, N, mat1, filenames, n, norm, minDepth, minLength, minCov, veccmp);
		
		/* start thread */
	} while(--i && !(errno = pthread_create(&thread->id, NULL, func, thread)));
	
	/* check if total number if threads were created */
	if(i) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		fprintf(stderr, "Will continue with %d threads.\n", tnum - i);
	}
	
	/* run function on main thread too */
	func(thread);
	
	/* join threads */
	joinThreads(thread->next);
	free(thread);
}

void * cmpMatRowThrd(void *arg) {
	
	static volatile int lock[1] = {0};
	static int pj = 0; /* sample, position in matrix */
	CmpMatArg *thread = arg;
	char *targetTemplate, *filename;
	unsigned j, n, norm, minDepth, minLength;
	double *D, *N, minCov, dist;
	double (*veccmp)(short unsigned*, short unsigned*, int, int);
	FileBuff *infile;
	MatrixCounts *mat1;
	NucCount *mat2;
	Qseqs **filenames;
	
	/* get input */
	D = (double *) thread->D;
	N = (double *) thread->N;
	mat1 = thread->mat1;
	filenames = thread->filenames;
	n = thread->n;
	norm = thread->norm;
	minDepth = thread->minDepth;
	minLength = thread->minLength;
	minCov = thread->minCov;
	veccmp = thread->veccmp;
	
	/* init */
	targetTemplate = (char *) mat1->name->seq;
	infile = setFileBuff(1048576);
	mat2 = initNucCount(128);
	
	/* calculate distances */
	while(pj < n) {
		lock(lock);
		j = pj++;
		unlock(lock);
		
		if(j < n) {
			filename = (char *) filenames[j]->seq;
			/* open matrix file, and find target */
			openAndDetermine(infile, filename);
			while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			
			/* get distance between the matrices */
			dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
			if(dist < 0) {
				if(dist == -1.0) {
					fprintf(stderr, "No sufficient overlap with sample:\t%s\n", filename);
				} else if(dist == -2.0) {
					fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filename);
					exit(1);
				} else {
					fprintf(stderr, "Failed to produce a distance metric for sample:\t%s\n", filename);
					exit(1);
				}
			}
			
			D[j] = dist;
			if(N) {
				N[j] = mat2->total;
			}
			
			/* close mtrix file */
			closeFileBuff(infile);
		}
	}
	
	/* clean */
	destroyFileBuff(infile);
	destroyNucCount(mat2);
	
	return NULL;
}

int ltdRowThrd(double *D, double *N, char *targetTemplate, char *addfilename, Qseqs **filenames, int n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int), int tnum) {
	
	FileBuff *infile;
	MatrixCounts *mat1;
	NucCount *mat2;
	
	/* init */
	infile = setFileBuff(1048576);
	mat1 = initMat(1048576, 128);
	mat2 = initNucCount(128);
	
	/* load new sample matrix into memory */
	/* open matrix file, and find target */
	openAndDetermine(infile, addfilename);
	while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
	
	/* initialize mat1 */
	setMatName(mat1, mat2);
	if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
		fprintf(stderr, "Malformed matrix in:\t%s\n", addfilename);
		exit(1);
	}
	closeFileBuff(infile);
	
	/* validate matrix */
	if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
		fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, addfilename);
		return 1;
	}
	
	/* strip matrix for insersions */
	stripMat(mat1);
	
	/* clean */
	destroyFileBuff(infile);
	destroyNucCount(mat2);
	
	/* thread out */
	if(n < tnum) {
		fprintf(stderr, "Adjustning number of nodes to %d, to conform with the matrix size.\n", (tnum = n));
	}
	matCmpThreadOut(tnum, &cmpMatRowThrd, (Matrix *) D, (Matrix *) N, mat1, filenames, n, norm, minDepth, minLength, minCov, veccmp);
	
	/* clean */
	destroyMat(mat1);
	
	return 0;
}
