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

#define _XOPEN_SOURCE 600
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bytescale.h"
#include "fbseek.h"
#include "filebuff.h"
#include "ltdmatrixthrd.h"
#include "matcmp.h"
#include "matparse.h"
#include "matrix.h"
#include "pherror.h"
#include "threader.h"
#define strcmp2(str1, str2) (*(str1) == *(str2) && strcmp(str1 + 1, str2 + 1) == 0)

static CmpMatArg * formThread(CmpMatArg *thread, Matrix *D, Matrix *N, MatrixCounts *mat1, Qseqs **filenames, unsigned n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int), unsigned char *include, TimeStamp **targetStamps, int srtd) {
	
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
	node->infile = setFileBuff(1048576);
	node->mat2 = initNucCount(128);
	node->include = include;
	node->targetStamps = targetStamps;
	node->srtd = srtd;
	
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
		destroyFileBuff(src->infile);
		destroyNucCount(src->mat2);
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
		thread = formThread(thread, D, N, mat1, filenames, n, norm, minDepth, minLength, minCov, veccmp, 0, 0, 1);
		
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
	destroyFileBuff(thread->infile);
	destroyNucCount(thread->mat2);
	free(thread);
}

void * cmpMatRowThrd(void *arg) {
	
	static volatile int Lock = 0;
	static int pj = 0; /* sample, position in matrix */
	volatile int *lock = &Lock;
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
	infile = thread->infile;
	mat2 = thread->mat2;
	
	/* init */
	targetTemplate = (char *) mat1->name->seq;
	
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
	
	return NULL;
}

void * cmpMatThrd(void *arg) {
	
	static volatile int pi = 0, cellFinish = 0, purge = 1, Lock = 0;
	static int pj = 0, sj = 0, srtd = 1; /* sample, position in matrix */
	volatile int *lock = &Lock;
	CmpMatArg *thread = arg;
	char *targetTemplate, *filename, **filenames;
	unsigned char *include;
	unsigned i, j, n, s, norm, minDepth, minLength;
	double *Dptr, *Nptr, minCov, dist;
	double (*veccmp)(short unsigned*, short unsigned*, int, int);
	float *Dfptr, *Nfptr;
	short unsigned *Dsptr, *Nsptr;
	unsigned char *Dbptr, *Nbptr;
	FileBuff *infile;
	Matrix *D, *N;
	MatrixCounts *mat1;
	NucCount *mat2;
	TimeStamp **targetStamps, **targetStamp;
	
	if(!thread) {
		pi = 0;
		cellFinish = 0;
		pj = 0;
		sj = 0;
		srtd = 1;
		purge = 1;
		return NULL;
	}
	
	/* get input */
	purge = 0;
	D = thread->D;
	N = thread->N;
	mat1 = thread->mat1;
	filenames = (char **) thread->filenames;
	n = thread->n;
	norm = thread->norm;
	minDepth = thread->minDepth;
	minLength = thread->minLength;
	minCov = thread->minCov;
	veccmp = thread->veccmp;
	infile = thread->infile;
	mat2 = thread->mat2;
	include = thread->include;
	targetStamps = thread->targetStamps;
	
	/* master */
	if(!n) {
		lock(lock);
		pj = 0;
		sj = 0;
		n = ++pi;
		cellFinish = n;
		srtd = thread->srtd;
		unlock(lock);
	}
	
	do {
		/* wait for next row to be ready */
		while(pi <= pj) {
			while(pi < pj) {
				usleep(256);
			}
			
			if(purge) {
				return NULL;
			}
		}
		
		while(pi <= pj) {
			usleep(256);
		}
		
		/* get next cell position */
		lock(lock);
		i = pi;
		j = pj++;
		s = sj++;
		if(!include[s]) {
			while(!include[++s]);
			sj = s + 1;
		}
		unlock(lock);
		
		/* init */
		Dptr = 0;
		Nptr = 0;
		Dfptr = 0;
		Nfptr = 0;
		Dsptr = 0;
		Nsptr = 0;
		Dbptr = 0;
		Nbptr = 0;
		if(D->mat) {
			Dptr = D->mat[i];
			Nptr = N ? N->mat[i] : 0;
		} else if(D->fmat) {
			Dfptr = D->fmat[i];
			Nfptr = N ? N->fmat[i] : 0;
		} else if(D->smat) {
			Dsptr = D->smat[i];
			Nsptr = N ? N->smat[i] : 0;
		} else {
			Dbptr = D->bmat[i];
			Nbptr = N ? N->bmat[i] : 0;
		}
		targetTemplate = (char *) mat1->name->seq;
		
		/* calculate distances */
		while(j < i) {
			
			targetStamp = targetStamps + s;
			filename = (char *) filenames[s];
			/* open matrix file, and find target */
			openAndDetermine(infile, filename);
			if(*targetStamp && srtd) {
				seekFileBiff(infile, *targetStamp);
				if(include[j] == 1) {
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				} else {
					/* set name */
					memcpy(mat2->name, mat1->name->seq, mat1->name->len);
				}
			} else {
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			}
			
			/* make timestamp */
			if(include[j] == 1) {
				*targetStamp = timeStampFileBuff(infile, *targetStamp);
				include[j] = 2;
			}
			
			/* get distance between the matrices */
			dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
			if(dist < 0) {
				if(dist == -1.0) {
					fprintf(stderr, "No sufficient overlap between samples:\t%s\t%s\n", filenames[i], filename);
				} else if(dist == -2.0) {
					fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filename);
					exit(1);
				} else {
					fprintf(stderr, "Failed to produce a distance metric for sample:\t%s\n", filename);
					exit(1);
				}
			}
			
			if(Dptr) {
				Dptr[j] = dist;
				if(N) {
					Nptr[j] = mat2->total;
				}
			} else if(Dfptr) {
				Dfptr[j] = dist;
				if(N) {
					Nfptr[j] = mat2->total;
				}
			} else if(Dsptr) {
				Dsptr[j] = dtouc(dist, 0.5);
				if(N) {
					Nsptr[j] = dtouc(mat2->total, 0.5);
				}
			} else {
				Dbptr[j] = dtouc(dist, 0.5);
				if(N) {
					Nbptr[j] = dtouc(mat2->total, 0.5);
				}
			}
			
			/* close mtrix file */
			closeFileBuff(infile);
			
			/* signal that cell is done */
			lock(lock);
			--cellFinish;
			j = pj++;
			s = sj++;
			if(!include[s]) {
				while(!include[++s]);
				sj = s + 1;
			}
			unlock(lock);
		}
	} while(i != n);
	
	/* wait for remaining cells to complete */
	while(cellFinish) {
		usleep(256);
	}
	
	return NULL;
}

void ltdMatrixThrd(Matrix *D, Matrix *N, MatrixCounts *mat1, TimeStamp **targetStamps, unsigned char *include, char *targetTemplate, char **filenames, int numFile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int), int tnum) {
	
	int i, n, srtd, len;
	CmpMatArg *thread;
	FileBuff *infile;
	NucCount *mat2;
	TimeStamp **targetStamp;
	
	/* thread out */
	n = -1;
	include += numFile;
	i = numFile + 1;
	while(--i) {
		if(*--include) {
			++n;
		}
	}
	if(n < tnum) {
		tnum = n;
		if(n <= 0) {
			return;
		}
	}
	i = tnum;
	thread = 0;
	do {
		/* make thread */
		thread = formThread(thread, D, N, mat1, (Qseqs **) filenames, n, norm, minDepth, minLength, minCov, veccmp, include, targetStamps, 1);
		
		/* start thread */
	} while(--i && !(errno = pthread_create(&thread->id, NULL, &cmpMatThrd, thread)));
	
	/* check if total number if threads were created */
	if(i) {
		fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
		fprintf(stderr, "Will continue with %d threads.\n", tnum - i);
	}
	thread->n = 0;
	infile = thread->infile;
	mat2 = thread->mat2;
	
	/* find first valid matrix */
	i = 0;
	do {
		if(include[i]) {
			/* open matrix file, and find target */
			openAndDetermine(infile, filenames[i]);
			if(*(targetStamp = targetStamps + i) && srtd) {
				seekFileBiff(infile, *targetStamp);
				/* first time template -> seek to new template */
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				if(!(strcmp2(targetTemplate, (char *) mat2->name))) {
					/* reset strm */
					sfseek(infile->file, 0, SEEK_SET);
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
					srtd = 0;
				}
			} else {
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			}
			
			if((strcmp2(targetTemplate, (char *) mat2->name))) {
				/* make timestamp */
				*targetStamp = timeStampFileBuff(infile, *targetStamp);
				include[i] = 2;
				
				n = 0;
				len = 0;
				while(FileBuffGetRow(infile, mat2) && mat2->ref) {
					if(mat2->ref != '-') {
						++len;
						/* validate position */
						if(minDepth <= mat2->total) {
							++n;
						}
					}
				}
				
				if(n < minLength || n < minCov * len) {
					fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filenames[i]);
					include[i] = 0;
				}
			} else {
				fprintf(stderr, "Template (\"%s\") is not included in:\t%s\n", targetTemplate, filenames[i]);
				include[i] = 0;
			}
			
			closeFileBuff(infile);
		}
	} while(!include[i++] && i < numFile);
	
	/* get distances */
	D->n = 0;
	while(i < numFile) {
		if(include[i]) {
			/* open matrix file, and find target */
			openAndDetermine(infile, filenames[i]);
			if(*(targetStamp = targetStamps + i) && srtd) {
				seekFileBiff(infile, *targetStamp);
				/* first time template -> seek to new template */
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				if(!(strcmp2(targetTemplate, (char *) mat2->name))) {
					/* reset strm */
					sfseek(infile->file, 0, SEEK_SET);
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
					srtd = 0;
				}
			} else {
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			}
			
			if((strcmp2(targetTemplate, (char *) mat2->name))) {
				/* make timestamp */
				if((*targetStamp = timeStampFileBuff(infile, *targetStamp))) {
					include[i] = 2;
				}
				
				/* initialize mat1 */
				setMatName(mat1, mat2);
				if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
					fprintf(stderr, "Input is not DB sorted.\n");
					/* reset strm */
					sfseek(infile->file, 0, SEEK_SET);
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
					
					if(strcmp2(targetTemplate, (char *) mat2->name)) {
						/* try load again */
						setMatName(mat1, mat2);
						if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
							fprintf(stderr, "Malformed matrix in:\t%s\n", filenames[i]);
							exit(1);
						}
					} else {
						fprintf(stderr, "Template (\"%s\") was not found in sample:\t%s\n", targetTemplate, filenames[i]);
						mat1->nNucs = 0;
						mat1->len = 0;
						include[i] = 0;
					}
					/* mark run as unsorted */
					srtd = 0;
				}
			} else {
				fprintf(stderr, "Template (\"%s\") is not included in:\t%s\n", targetTemplate, filenames[i]);
				include[i] = 0;
			}
			closeFileBuff(infile);
			
			/* validate matrix */
			if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
				fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filenames[i]);
				include[i] = 0;
			} else {
				/* strip matrix for insersions */
				stripMat(mat1);
			}
			
			if(include[i]) {
				thread->srtd = srtd;
				cmpMatThrd(thread);
			}
		}
		++i;
	}
	
	/* get number of include templates */
	++numFile;
	--include;
	n = 0;
	while(--numFile) {
		if(*++include) {
			++n;
		}
	}
	D->n = n;
	if(N) {
		N->n = n;
	}
	
	/* reset function */
	cmpMatThrd(0);
	
	/* join threads */
	joinThreads(thread->next);
	destroyFileBuff(thread->infile);
	destroyNucCount(thread->mat2);
	free(thread);
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
