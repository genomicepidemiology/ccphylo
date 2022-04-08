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
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "bytescale.h"
#include "filebuff.h"
#include "fsacmp.h"
#include "fsacmpthrd.h"
#include "matrix.h"
#include "pherror.h"
#include "threader.h"
#include "seqparse.h"

static CmpFsaArg * formThread(CmpFsaArg *thread, Matrix *D, Matrix *N, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, unsigned minLength, double minCov, FILE *diffile, char *targetTemplate, Qseqs *ref, Qseqs **filenames, unsigned proxi) {
	
	CmpFsaArg *node;
	
	/* make thread */
	node = smalloc(sizeof(CmpFsaArg));
	node->D = D;
	node->N = N;
	node->n = n;
	node->len = len;
	node->seqs = seqs;
	node->include = include;
	node->includes = includes;
	node->norm = norm;
	node->minLength = minLength;
	node->minCov = minCov;
	node->diffile = diffile;
	node->targetTemplate = targetTemplate;
	node->ref = ref;
	node->filenames = filenames;
	node->proxi = proxi;
	
	/* link to remaining threads */
	node->next = thread;
	
	return node;
}

static void joinThreads(CmpFsaArg *src) {
	
	CmpFsaArg *src_next;
	
	while(src) {
		src_next = src->next;
		if((errno = pthread_join(src->id, NULL))) {
			ERROR();
		}
		free(src);
		src = src_next;
	}
}

void fsaCmpThreadOut(int tnum, void * (*func)(void*), Matrix *D, Matrix *N, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, unsigned minLength, double minCov, FILE *diffile, char *targetTemplate, Qseqs *ref, Qseqs **filenames, unsigned proxi) {
	
	int i;
	CmpFsaArg *thread;
	
	/* thread out */
	i = tnum;
	thread = 0;
	do {
		/* make thread */
		thread = formThread(thread, D, N, n, len, seqs, include, includes, norm, minLength, minCov, diffile, targetTemplate, ref, filenames, proxi);
		
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
	
	/* reset function */
	func(0);
}

void * cmpFsaThrd(void *arg) {
	
	/* this is the threaded version of cmpFsa */
	static volatile int Lock = 0;
	static int si = 1, sj = 0, pi = 0, pj = 0, inc = 0; /* sample, position in matrix */
	volatile int *lock = &Lock;
	CmpFsaArg *thread = arg;
	int n, len, Dn;
	unsigned i, j, norm, *includes;
	long unsigned **seqs, *seqi, *seqj, dist;
	unsigned char *include;
	double *Dptr, nFactor;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	FILE *diffile;
	Matrix *D;
	
	/* get input */
	if(thread) {
		D = thread->D;
		n = thread->n;
		len = thread->len;
		seqs = thread->seqs;
		include = thread->include;
		includes = *(thread->includes);
		norm = thread->norm;
		diffile = thread->diffile;
	} else {
		/* reset function */
		si = 1;
		sj = 0;
		pi = 0;
		pj = 0;
		inc = 0;
		return NULL;
	}
	
	/* set stuff that only needs to be set once */
	lock(lock);
	if(!pi) {
		/* seek first included sample */
		i = 0;
		while(i != n && *include == 0) {
			++include;
			++i;
		}
		si = i + 1;
		
		/* get number of included samples */
		Dn = *include;
		while(++i != n) {
			Dn += *++include;
		}
		include -= (n - 1);
		D->n = Dn;
		inc = getNpos(includes, len);
		fprintf(stderr, "# %d / %d bases included in distance matrix.\n", inc, len);
		pi = 1;
	}
	unlock(lock);
	
	/* get constant stuff */
	if(norm) {
		nFactor = norm;
		nFactor /= inc;
	} else {
		nFactor = 1.0;
	}
	Dn = D->n;
	if(Dn < 2) {
		return NULL;
	}
	
	/* calculate distances */
	while(pi != Dn) {
		lock(lock);
		/* last distance been found or is being computed */
		if(pi == Dn) {
			unlock(lock);
			return NULL;
		}
		
		/* get next valid sample combo */
		i = si;
		j = sj;
		while(include[i] == 0 && include[j] == 0) {
			while(i < n && include[i] == 0) {
				++i;
			}
			while(j < i && include[j] == 0) {
				++j;
			}
			
			if(i == j) {
				++i;
				j = 0;
			}
		}
		
		/* get sequences */
		seqi = seqs[i];
		seqj = seqs[j];
		
		/* update si, sj */
		si = i;
		sj = j + 1;
		if(si == sj) {
			++si;
			sj = 0;
		}
		
		/* get the corresponding position in distance matrix */
		Dptr = 0;
		Dfptr = 0;
		Dsptr = 0;
		Dbptr = 0;
		if(D->mat) {
			Dptr = D->mat[pi] + pj;
		} else if(D->fmat) {
			Dfptr = D->fmat[pi] + pj;
		} else if(D->smat) {
			Dsptr = D->smat[pi] + pj;
		} else {
			Dbptr = D->bmat[pi] + pj;
		}
		
		/* update next position in D */
		if(pi == ++pj) {
			++pi;
			pj = 0;
		}
		unlock(lock);
		
		if(diffile) {
			dist = fsacmprint(diffile, i, j, seqi, seqj, includes, len);
		} else {
			dist = fsacmp(seqi, seqj, includes, len);
		}
		if(Dptr) {
			*Dptr = nFactor * dist;
		} else if(Dfptr) {
			*Dfptr = nFactor * dist;
		} else if(Dsptr) {
			*Dsptr = dtouc(nFactor * dist, 0.5);
		} else {
			*Dbptr = dtouc(nFactor * dist, 0.5);
		}
	}
	
	return NULL;
}

void * cmpairFsaThrd(void *arg) {
	
	/* this is the threaded version of cmpairFsa */
	static volatile int Lock = 0;
	static int si = 1, sj = 0, pi = 0, pj = 0, base = 0; /* sample, position in matrix */
	volatile int *lock = &Lock;
	CmpFsaArg *thread = arg;
	int n, len, Dn;
	unsigned i, j, inc, norm, minLength, proxi;
	unsigned **includes, *includesi, *includesj, *includesPair;
	long unsigned **seqs, *seqi, *seqj, dist;
	unsigned char *include, *Dbptr, *Nbptr;
	double minCov, *Dptr, *Nptr;
	float *Dfptr, *Nfptr;
	short unsigned *Dsptr, *Nsptr;
	FILE *diffile;
	Matrix *D, *N;
	
	/* get input */
	if(thread) {
		D = thread->D;
		N = thread->N;
		n = thread->n;
		len = thread->len;
		seqs = thread->seqs;
		include = thread->include;
		includes = thread->includes;
		norm = thread->norm;
		minLength = thread->minLength;
		minCov = thread->minCov;
		diffile = thread->diffile;
		minLength = minLength < minCov * len ? minCov * len : minLength;
		proxi = thread->proxi;
	} else {
		/* reset function */
		si = 1;
		sj = 0;
		pi = 0;
		pj = 0;
		base = 0;
		return NULL;
	}
	
	/* set stuff that only needs to be set once */
	lock(lock);
	if(!pi) {
		/* seek first included sample */
		i = 0;
		while(i != n && *include == 0) {
			++include;
			++i;
		}
		base = i;
		si = i + 1;
		sj = i;
		
		/* get number of included samples */
		Dn = *include;
		while(++i != n) {
			Dn += *++include;
		}
		include -= (n - 1);
		D->n = Dn;
		if(N) {
			N->n = Dn;
		}
		pi = 1;
	}
	unlock(lock);
	
	/* get actual number of samples in D */
	Dn = D->n;
	if(Dn < 2) {
		return NULL;
	}
	includesPair = smalloc((len / 32 + 1) * sizeof(unsigned));
	
	/* calculate distances */
	while(pi != Dn) {
		lock(lock);
		/* last distance been found or is being computed */
		if(pi == Dn) {
			unlock(lock);
			free(includesPair);
			return NULL;
		}
		
		/* get next valid sample combo */
		i = si;
		j = sj;
		while(include[i] == 0 || include[j] == 0) {
			while(i < n && include[i] == 0) {
				++i;
			}
			while(j < i && include[j] == 0) {
				++j;
			}
			
			if(i == j) {
				++i;
				j = base;
			}
		}
		
		/* get sequences */
		seqi = seqs[i];
		seqj = seqs[j];
		includesi = includes[i];
		includesj = includes[j];
		
		/* update si, sj */
		si = i;
		sj = j + 1;
		if(si == sj) {
			++si;
			sj = base;
		}
		
		/* get the corresponding position in distance matrix */
		Dptr = 0;
		Nptr = 0;
		Dfptr = 0;
		Nfptr = 0;
		Dsptr = 0;
		Nsptr = 0;
		Dbptr = 0;
		Nbptr = 0;
		if(D->mat) {
			Dptr = D->mat[pi] + pj;
			Nptr = N ? N->mat[pi] + pj : 0;
		} else if(D->fmat) {
			Dfptr = D->fmat[pi] + pj;
			Nfptr = N ? N->fmat[pi] + pj : 0;
		} else if(D->smat) {
			Dsptr = D->smat[pi] + pj;
			Nsptr = N ? N->smat[pi] + pj : 0;
		} else {
			Dbptr = D->bmat[pi] + pj;
			Nbptr = N ? N->bmat[pi] + pj : 0;
		}
		
		/* update next position in D */
		if(pi == ++pj) {
			++pi;
			pj = 0;
		}
		unlock(lock);
		
		/* mask out proximity SNPs */
		maskProxi(includesPair, includesi, includesj, seqi, seqj, len, proxi);
		
		/* calculate distance of pair */
		if(diffile) {
			dist = fsacmpairint(diffile, i, j, seqi, seqj, includesPair, len);
		} else {
			dist = fsacmpair(seqi, seqj, includesPair, len);
		}
		
		/* separate distance and included bases */
		inc = dist & UINT_MAX;
		if(Dptr) {
			if(minLength <= inc) {
				if(norm) {
					*Dptr = (dist >> 32) * norm;
					*Dptr /= inc;
				} else {
					*Dptr = (dist >> 32);
				}
			} else {
				*Dptr = -1.0;
			}
			if(N) {
				*Nptr = inc;
			}
		} else if(Dfptr) {
			if(minLength <= inc) {
				if(norm) {
					*Dfptr = (dist >> 32) * norm;
					*Dfptr /= inc;
				} else {
					*Dfptr = (dist >> 32);
				}
			} else {
				*Dfptr = -1.0;
			}
			if(N) {
				*Nfptr = inc;
			}
		} else if(Dsptr) {
			if(minLength <= inc) {
				if(norm) {
					*Dsptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
				} else {
					*Dsptr = dtouc((dist >> 32), 0.5);
				}
			} else {
				*Dsptr = dtouc(-1.0, 0);
			}
			if(N) {
				*Nsptr = dtouc(inc, 0.5);
			}
		} else {
			if(minLength <= inc) {
				if(norm) {
					*Dbptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
				} else {
					*Dbptr = dtouc((dist >> 32), 0.5);
				}
			} else {
				*Dbptr = dtouc(-1.0, 0);
			}
			if(N) {
				*Nbptr = dtouc(inc, 0.5);
			}
		}
	}
	
	free(includesPair);
	return NULL;
}

void * cmpFsaRowThrd(void *arg) {
	
	static volatile int Lock = 0;
	static int pj = 0; /* sample, position in matrix */
	volatile int *lock = &Lock;
	CmpFsaArg *thread = arg;
	int n, len;
	unsigned j, minLength, inc, proxi;
	unsigned *includeadd, *includeseq;
	long unsigned *addL, *seqL, dist;
	char *targetTemplate, *filename;
	unsigned char *trans;
	double minCov, norm, *D, *N;
	FILE *diffile;
	FileBuff *infile;
	Qseqs *ref, *header, *seq, **filenames;
	
	/* get input */
	if(thread) {
		D = (double *) thread->D;
		N = (double *) thread->N;
		n = thread->n;
		len = thread->len;
		addL = *(thread->seqs);
		trans = thread->include;
		includeadd = *(thread->includes);
		norm = thread->norm;
		minLength = thread->minLength;
		minCov = thread->minCov;
		diffile = thread->diffile;
		minLength = minLength < minCov * len ? minCov * len : minLength;
		targetTemplate = thread->targetTemplate;
		ref = thread->ref;
		filenames = thread->filenames;
		proxi = thread->proxi;
	} else {
		/* reset function */
		pj = 0;
		return NULL;
	}
	
	/* init */
	infile = setFileBuff(1048576);
	header = setQseqs(32);
	seq = setQseqs(ref->size);
	includeseq = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	seqL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	
	/* calculate distances */
	while(pj < n) {
		lock(lock);
		j = pj++;
		unlock(lock);
		
		if(j < n) {
			filename = (char *) filenames[j]->seq;
			/* open matrix file, and find target */
			openAndDetermine(infile, filename);
			while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
			FileBuffgetFsaSeq(infile, seq, trans);
			closeFileBuff(infile);
			
			/* get included positions */
			memcpy(includeseq, includeadd, ((len >> 5) + 1) * sizeof(unsigned));
			getIncPosPtr(includeseq, seq, ref, proxi);
			
			/* convert seq to nibbles */
			qseq2nibble(seq, seqL);
			
			/* get distance */
			if(diffile) {
				dist = fsacmpairint(diffile, n, j, addL, seqL, includeseq, len);
			} else {
				dist = fsacmpair(addL, seqL, includeseq, len);
			}
			
			/* separate distance and included bases */
			if(minLength <= (inc = dist & UINT_MAX)) {
				D[j] = norm ? ((dist >> 32) * norm / inc) : (dist >> 32);
			} else {
				D[j] = -1.0;
				inc = 0;
				fprintf(stderr, "No sufficient overlap with sample:\t%s\n", filename);
			}
			if(N) {
				N[j] = inc;
			}
		}
	}
	
	/* clean */
	destroyFileBuff(infile);
	destroyQseqs(header);
	destroyQseqs(seq);
	free(includeseq);
	free(seqL);
	
	return NULL;
}

int ltdFsaRowThrd(double *D, double *N, char *targetTemplate, char *addfilename, char *diffilename, Qseqs **filenames, int n, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, int tnum) {
	
	/* this is the threaded version of ltdFsaRow_get */
	unsigned char *trans;
	unsigned *includeadd, len;
	long unsigned *addL;
	FILE *diffile;
	FileBuff *infile;
	Qseqs *header, *seq;
	
	/* init  */
	len = 0;
	trans = get2BitTable(flag);
	infile = setFileBuff(1048576);
	header = setQseqs(32);
	seq = setQseqs(1048576);
	
	/* load new sample into memory */
	/* open file, and find target */
	openAndDetermine(infile, addfilename);
	
	/* find the correct entry */
	while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
	
	/* get sequence */
	if(!header->len) {
		fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, addfilename);
		exit(1);
	} else if(FileBuffgetFsaSeq(infile, seq, trans)) {
		closeFileBuff(infile);
	} else {
		fprintf(stderr, "\"%s\" is not fasta.\n", addfilename);
		exit(1);
	}
	if(len != 0 && len != seq->len) {
		fprintf(stderr, "New sequence does not match the existing sequences.\n");
		exit(1);
	} else {
		len = seq->len;
		minLength = minLength < minCov * len ? minCov * len : minLength;
	}
	
	/* validate new seq */
	includeadd = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	addL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	initIncPos(includeadd, len);
	getIncPosPtr(includeadd, seq, seq, proxi);
	if(getNpos(includeadd, len) < minLength) {
		fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, addfilename);
		return 1;
	} else if(diffilename) {
		if(*diffilename == '-' && diffilename[1] == 0) {
			diffile = stdout;
		} else {
			diffile = sfopen(diffilename, "ab");
		}
	} else {
		diffile = 0;
	}
	
	/* convert seq to nibbles */
	qseq2nibble(seq, addL);
	
	/* clean */
	destroyQseqs(header);
	destroyFileBuff(infile);
	
	/* thread out */
	if(n < tnum) {
		fprintf(stderr, "Adjustning number of nodes to %d, to conform with the matrix size.\n", (tnum = n));
	}
	fsaCmpThreadOut(tnum, &cmpFsaRowThrd, (Matrix *) D, (Matrix *) N, n, len, &addL, trans, &includeadd, norm, minLength, minCov, diffile, targetTemplate, seq, filenames, proxi);
	
	/* update size of matrix */
	if(diffile && diffile != stdout) {
		fclose(diffile);
	}
	
	/* clean */
	destroyQseqs(seq);
	destroyTable(trans);
	free(includeadd);
	free(addL);
	
	return 0;
}
