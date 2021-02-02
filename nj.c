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
#include <stdlib.h>
#include "matrix.h"
#include "nj.h"
#include "nwck.h"
#include "qseqs.h"
#include "pherror.h"
#include "str.h"
#include "threader.h"
#include "vector.h"

void (*limbLengthPtr)(double*, double*, unsigned, unsigned, Vector *, unsigned*, double) = &limbLength;
double (*initQchunkPtr)(Matrix *, double*, unsigned*, double, int*, int*, int, int) = &initQchunk;
double (*minDchunkPtr)(Matrix *, unsigned*, double, int*, int*, int, int) = &minDchunk;
void (*updateDptr)(Matrix *, Vector *, unsigned *, unsigned, unsigned, double, double) = &updateD;
long unsigned (*minDist)(Matrix *, Vector *, unsigned *) = &initQ;
void * (*minDist_thread)(void *) = &initQ_thread;

void limbLength(double *Li, double *Lj, unsigned i, unsigned j, Vector *sD, unsigned *N, double D_ij) {
	
	unsigned Ni, Nj;
	double delta_ij;
	
	Ni = N[i] - 2;
	Nj = N[j] - 2;
	if(0 < Ni && 0 < Nj) {
		/* semi-standard */
		delta_ij = ((sD->vec[i] - D_ij) / Ni) - ((sD->vec[j] - D_ij) / Nj);
		*Li = (D_ij + delta_ij) / 2;
		*Lj = (D_ij - delta_ij) / 2;
		
		/* org:
		delta_ij = (sD->vec[i] - sD->vec[j]) / N = (sD->vec[i] - D_ij) / N) - ((sD->vec[j] - D_ij) / N;
		*/
		
		/* avoid negative branch lengths */
		if(*Li < 0) {
			*Lj = D_ij;
			*Li = 0;
		} else if(*Lj < 0) {
			*Li = D_ij;
			*Lj = 0;
		}
	} else if(0 < Ni) {
		/* j only shares similarity with i */
		*Li = 0;
		*Lj = D_ij;
	} else if(0 < Nj) {
		/* i only shares similarity with j */
		*Li = D_ij;
		*Lj = 0;
	} else {
		/* i and j only share similarity with eachother */
		*Li = (*Lj = D_ij / 2);
	}
}

void limbLengthNeg(double *Li, double *Lj, unsigned i, unsigned j, Vector *sD, unsigned *N, double D_ij) {
	
	unsigned Ni, Nj;
	double delta_ij;
	
	Ni = N[i] - 2;
	Nj = N[j] - 2;
	if(0 < Ni && 0 < Nj) {
		/* semi-standard */
		delta_ij = ((sD->vec[i] - D_ij) / Ni) - ((sD->vec[j] - D_ij) / Nj);
		*Li = (D_ij + delta_ij) / 2;
		*Lj = (D_ij - delta_ij) / 2;
		
		/* org:
		delta_ij = (sD->vec[i] - sD->vec[j]) / N = (sD->vec[i] - D_ij) / N) - ((sD->vec[j] - D_ij) / N;
		*/
	} else if(0 < Ni) {
		/* j only shares similarity with i */
		*Li = 0;
		*Lj = D_ij;
	} else if(0 < Nj) {
		/* i only shares similarity with j */
		*Li = D_ij;
		*Lj = 0;
	} else {
		/* i and j only share similarity with eachother */
		*Li = (*Lj = D_ij / 2);
	}
}

unsigned * initSummaD(Vector *sD, Matrix *D, unsigned *N) {
	
	int i, j;
	unsigned *Nptr;
	double dist, *sDptr, *sDvec;
	
	/*
	sD(i) = summa(d(i,k));
	*/
	
	/* alloc */
	if(sD->size < (sD->n = D->n)) {
		vector_realloc(sD, sD->n);
		if(!(N = realloc(N, sD->n * sizeof(unsigned)))) {
				ERROR();
		}
	}
	
	/* init to zero */
	i = sD->n + 1;
	sDptr = sD->vec - 1;
	Nptr = N - 1;
	while(--i) {
		*++sDptr = 0;
		*++Nptr = 1;
	}
	
	/* compute sum's */
	sDptr = *(D->mat) - 1;
	sDvec = sD->vec;
	for(i = 1; i < sD->n; ++i) {
		for(j = 0; j < i; j++) {
			if(0 <= (dist = *++sDptr)) {
				sDvec[i] += dist;
				sDvec[j] += dist;
				++N[i];
				++N[j];
			}
		}
	}
	
	return N;
}

long unsigned initQ(Matrix *D, Vector *sD, unsigned *N) {
	
	int i, j, mi, mj;
	unsigned N_i, *Ni, *Nj;
	long unsigned pos;
	double dist, min, Q, sD_i, *sDi, *sDj, *Dptr;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	mi = 0;
	mj = 0;
	min = 1;
	Dptr = *(D->mat) - 1;
	sDi = sD->vec;
	Ni = N;
	i = 0;
	while(++i < D->n) {
		N_i = *++Ni;
		Nj = N - 1;
		sD_i = *++sDi;
		sDj = sD->vec - 1;
		j = -1;
		while(++j < i) {
			if(0 <= (dist = *++Dptr)) {
				dist *= ((N_i + *++Nj - 4) >> 1);
				if((Q = dist - sD_i - *++sDj) < min) {
					min = Q;
					mi = i;
					mj = j;
				}
			} /* else {
				// missing
				Q = 1;
			} */
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double initQchunk(Matrix *D, double *sD, unsigned *N, double min, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end;
	unsigned N_i;
	double sD_i, dist, Q, *Dptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	N_i = N[i];
	N += --j;
	sD_i = sD[i];
	sD += j;
	Dptr = D->mat[i] + j;
	
	/* get chunk */
	while(++j < end) {
		if(0 <= (dist = *++Dptr)) {
			dist *= ((N_i + *++N - 4) >> 1);
			if((Q = dist - sD_i - *++sD) < min) {
				min = Q;
				*mi = i;
				*mj = j;
			}
		}
	}
	
	return min;
}

long unsigned initQ_MN(Matrix *D, Vector *sD, unsigned *N) {
	
	int i, j, mi, mj;
	unsigned N_i, *Ni, *Nj;
	long unsigned pos;
	double dist, max, Q, sD_i, *sDi, *sDj, *Dptr;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	mi = 0;
	mj = 0;
	max = -1 - sD->vec[1] - sD->vec[2];
	Dptr = *(D->mat) - 1;
	sDi = sD->vec;
	Ni = N;
	i = 0;
	while(++i < D->n) {
		N_i = *++Ni;
		Nj = N - 1;
		sD_i = *++sDi;
		sDj = sD->vec - 1;
		j = -1;
		while(++j < i) {
			if(0 <= (dist = *++Dptr)) {
				dist *= ((N_i + *++Nj - 4) >> 1);
				if(max < (Q = dist - sD_i - *++sDj)) {
					max = Q;
					mi = i;
					mj = j;
				}
			} /* else {
				// missing
				Q = 1;
			} */
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double initQ_MNchunk(Matrix *D, double *sD, unsigned *N, double max, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end;
	unsigned N_i;
	double sD_i, dist, Q, *Dptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	N_i = N[i];
	N += --j;
	sD_i = sD[i];
	sD += j;
	Dptr = D->mat[i] + j;
	
	/* get chunk */
	while(++j < end) {
		if(0 <= (dist = *++Dptr)) {
			dist *= ((N_i + *++N - 4) >> 1);
			if(max < (Q = dist - sD_i - *++sD)) {
				max = Q;
				*mi = i;
				*mj = j;
			}
		}
	}
	
	return max;
}

void * initQ_thread(void *arg) {
	
	const int chunk = 65536;
	static volatile int thread_begin = 1, thread_wait, lock[1] = {0};
	static int thread_num = 1, next_i, next_j, Mi, Mj;
	static double Min;
	NJthread *thread = arg;
	int i, j, mi, mj;
	unsigned *N;
	double min, *sDvec;
	Matrix *D;
	Vector *sD;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	if(!thread) {
		wait_atomic((thread_begin != thread_num));
		thread_num = 0;
		thread_begin = 0;
		return NULL;
	}
	D = thread->D;
	sD = thread->sD;
	N = thread->N;
	if(thread->min) {
		lock(lock);
		wait_atomic((thread_begin != thread_num));
		/* init start */
		Mi = 0;
		Mj = 0;
		Min = (initQchunkPtr == &initQchunk) ? 1 : -1 - sD->vec[1] - sD->vec[2];
		next_i = 1;
		next_j = 0;
		thread_num = thread->num;
		thread_wait = thread_num;
		thread_begin = 0;
		unlock(lock);
	}
	
	do {
		wait_atomic(thread_begin);
		
		/* terminate */
		if(!thread_num) {
			return NULL;
		}
		
		/* init */
		i = 0;
		mi = 0;
		mj = 0;
		min = (initQchunkPtr == &initQchunk) ? 1 : (-1 - sD->vec[1] - sD->vec[2]);
		sDvec = sD->vec;
		
		/* fill in chunks */
		while(i < D->n) {
			/* get next chunk */
			lockTime(lock, 2);
			i = next_i;
			j = next_j;
			if(i <= (next_j += chunk)) {
				++next_i;
				next_j = 0;
			}
			unlock(lock);
			
			/* init chunk */
			if(i < D->n) {
				min = initQchunkPtr(D, sDvec, N, min, &mi, &mj, i, j);
			}
		}
		
		/* check min */
		lockTime(lock, 2);
		if(initQchunkPtr == &initQ_MNchunk) {
			if(Min < min) {
				Min = min;
				Mi = mi;
				Mj = mj;
			} else if(Min == min && (mi < Mi || (mi == Mi && mj < Mj))) {
				/* ensure reproducibility */
				Mi = mi;
				Mj = mj;
			}
		} else if(min < Min) {
			Min = min;
			Mi = mi;
			Mj = mj;
		} else if(min == Min && (mi < Mi || (mi == Mi && mj < Mj))) {
			/* ensure reproducibility */
			Mi = mi;
			Mj = mj;
		}
		--thread_wait;
		unlock(lock);
		wait_atomic(thread_wait);
		__sync_add_and_fetch(&thread_begin, 1);
	} while(!thread->min);
	
	thread->min = Min;
	thread->mi = Mi;
	thread->mj = Mj;
	
	return NULL;
}

long unsigned minD(Matrix *D, Vector *sD, unsigned *N) {
	
	unsigned i, j, mi, mj;
	long unsigned pos;
	double min, *Dptr;
	
	mi = 0;
	mj = 0;
	Dptr = *(D->mat) - 1;
	min = (minDchunkPtr == &minDchunk) ? 2 * **(D->mat) : -1;
	for(i = 1; i < D->n; ++i) {
		for(j = 0; j < i; ++j) {
			if(0 <= *++Dptr && *Dptr < min) {
				min = *Dptr;
				mi = i;
				mj = j;
			}
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double minDchunk(Matrix *D, unsigned *N, double min, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end;
	double *Dptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	Dptr = D->mat[i] + --j;
	
	/* get chunk */
	while(++j < end) {
		if(0 <= *++Dptr && *Dptr < min) {
			min = *Dptr;
			*mi = i;
			*mj = j;
		}
	}
	
	return min;
}

long unsigned maxD(Matrix *D, Vector *sD, unsigned *N) {
	
	unsigned i, j, mi, mj;
	long unsigned pos;
	double max, *Dptr;
	
	mi = 0;
	mj = 0;
	max = -1;
	Dptr = *(D->mat) - 1;
	for(i = 1; i < D->n; ++i) {
		for(j = 0; j < i; ++j) {
			if(0 <= *++Dptr && max < *Dptr) {
				max = *Dptr;
				mi = i;
				mj = j;
			}
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double maxDchunk(Matrix *D, unsigned *N, double max, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end;
	double *Dptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	Dptr = D->mat[i] + --j;
	
	/* get chunk */
	while(++j < end) {
		if(0 <= *++Dptr && max < *Dptr) {
			max = *Dptr;
			*mi = i;
			*mj = j;
		}
	}
	
	return max;
}

void * minD_thread(void *arg) {
	
	const int chunk = 65536;
	static volatile int thread_begin = 1, thread_wait, lock[1] = {0};
	static int thread_num = 1, next_i, next_j, Mi, Mj;
	static double Min;
	NJthread *thread = arg;
	int i, j, mi, mj;
	unsigned *N;
	double min;
	Matrix *D;
	
	/* init */
	if(!thread) {
		wait_atomic((thread_begin != thread_num));
		thread_num = 0;
		thread_begin = 0;
		return NULL;
	}
	D = thread->D;
	N = thread->N;
	if(thread->min) {
		lock(lock);
		/* init start */
		Mi = 0;
		Mj = 0;
		Min = (minDchunkPtr == &minDchunk) ? 2 * **(D->mat) : -1;
		next_i = 1;
		next_j = 0;
		wait_atomic((thread_begin != thread_num));
		thread_num = thread->num;
		thread_wait = thread_num;
		thread_begin = 0;
		unlock(lock);
	}
	
	do {
		wait_atomic(thread_begin);
		/* terminate */
		if(!thread_num) {
			return NULL;
		}
		
		/* init */
		i = 0;
		mi = 0;
		mj = 0;
		min = (minDchunkPtr == &minDchunk) ? (2 * **(D->mat)) : -1;
		
		/* fill in chunks */
		while(i < D->n) {
			/* get next chunk */
			lockTime(lock, 2);
			i = next_i;
			j = next_j;
			if(i <= (next_j += chunk)) {
				++next_i;
				next_j = 0;
			}
			unlock(lock);
			
			/* init chunk */
			if(i < D->n) {
				min = minDchunkPtr(D, N, min, &mi, &mj, i, j);
			}
		}
		
		/* check min */
		lockTime(lock, 2);
		if(minDchunkPtr == &maxDchunk) {
			if(Min < min) {
				Min = min;
				Mi = mi;
				Mj = mj;
			} else if(Min == min && (mi < Mi || (mi == Mi && mj < Mj))) {
				/* ensure reproducibility */
				Mi = mi;
				Mj = mj;
			}
		} else if(min < Min) {
			Min = min;
			Mi = mi;
			Mj = mj;
		} else if(min == Min && (mi < Mi || (mi == Mi && mj < Mj))) {
			/* ensure reproducibility */
			Mi = mi;
			Mj = mj;
		}
		--thread_wait;
		unlock(lock);
		wait_atomic(thread_wait);
		__sync_add_and_fetch(&thread_begin, 1);
	} while(!thread->min);
	
	thread->min = Min;
	thread->mi = Mi;
	thread->mj = Mj;
	
	return NULL;
}

long unsigned minPair(Matrix *Q) {
	
	unsigned i, j, mi, mj;
	long unsigned pos;
	double min, *Qptr;
	
	mi = 0;
	mj = 0;
	min = 1;
	Qptr = *(Q->mat) - 1;
	for(i = 1; i < Q->n; ++i) {
		for(j = 0; j < i; j++) {
			if(*++Qptr < min) {
				min = *Qptr;
				mi = i;
				mj = j;
			}
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

void updateD(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj) {
	
	unsigned k, n, Dn;
	double dist, sd, D_ik, D_kj, D_ij, *D_i, *D_j, **Dmat, *sDvec;
	
	/*
	D(i,j) = D(j,i) <-> ltd -> j < i
	D(k,m) = D(i,k) + D(k,j) - D(i,j)
	! D(i,k) -> D(k,m) = D(k,j) - L(j)
	! D(k,j) -> D(k,m) = D(i,k) - L(i)
	*/
	
	/* init */
	Dmat = D->mat;
	D_i = Dmat[i];
	D_j = Dmat[j];
	D_ij = D_i[j];
	sDvec = sD->vec;
	
	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;
	
	/* update (first) row */
	for(k = 0; k < j; ++k) {
		D_ik = D_i[k];
		D_kj = D_j[k];
		if(0 <= D_ik && 0 <= D_kj) {
			dist = (D_ik + D_kj - D_ij) / 2;
			D_j[k] = dist;
			/* update N and sD */
			sDvec[k] -= (D_kj - dist);
			sd += dist;
			--N[k];
			++n;
		} else if(0 <= D_ik) {
			dist = D_ik - Li;
			D_j[k] = dist;
			/* update N and sD, sD(k) is new */
			sDvec[k] += dist;
			sd += dist;
			++n;
		} else if(0 <= D_kj) {
			dist = (D_j[k] -= Lj);
			/* update N and sD */
			sDvec[k] += (dist - D_kj);
			--N[k];
			sd += dist;
			++n;
		}
		sDvec[k] -= D_ik;
	}
	
	/* update jth column */
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			Dn = D->n;
		}
		while(++k < Dn) {
			D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
			D_kj = Dmat[k][j];
			if(0 <= D_ik && 0 <= D_kj) {
				dist = (D_kj + D_ik - D_ij) / 2;
				Dmat[k][j] = dist < 0 ? 0 : dist;
				
				/* update N and sD */
				sDvec[k] -= (D_kj - dist);
				sd += dist;
				--N[k];
				++n;
			} else if(0 <= D_ik) {
				dist = D_ik - Li;
				Dmat[k][j] = dist;
				/* update N and sD, sD(k) is new */
				sDvec[k] += dist;
				sd += dist;
				++n;
			} else if(0 <= D_kj) {
				dist = (Dmat[k][j] -= Lj);
				/* update N and sD */
				sDvec[k] += (dist - D_j[k]);
				--N[k];
				sd += dist;
				++n;
			}
			sDvec[k] -= D_ik;
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sDvec[j] = sd;
}

void updateD_UPGMA(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj) {
	
	unsigned k, n, Dn;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	
	/*
	UPGMA
	*/
	
	/* init */
	Dmat = D->mat;
	D_i = Dmat[i];
	D_j = Dmat[j];
	sDvec = sD->vec;
	
	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;
	
	/* update (first) row */
	for(k = 0; k < j; ++k) {
		D_ik = D_i[k];
		D_kj = D_j[k];
		if(0 <= D_ik && 0 <= D_kj) {
			dist = (D_ik + D_kj) / 2;
			D_j[k] = dist;
			/* update N and sD */
			sDvec[k] -= (D_kj - dist);
			sd += dist;
			--N[k];
			++n;
		} else if(0 <= D_ik) {
			D_j[k] = D_ik;
			/* update N and sD, sD(k) is new */
			sDvec[k] += D_ik;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			--N[k];
			sd += D_kj;
			++n;
		}
		sDvec[k] -= D_ik;
	}
	
	/* update jth column */
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			Dn = D->n;
		}
		while(++k < Dn) {
			D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
			D_kj = Dmat[k][j];
			if(0 <= D_ik && 0 <= D_kj) {
				dist = (D_kj + D_ik) / 2;
				Dmat[k][j] = dist;
				
				/* update N and sD */
				sDvec[k] -= (D_kj - dist);
				sd += dist;
				--N[k];
				++n;
			} else if(0 <= D_ik) {
				Dmat[k][j] = D_ik;
				/* update N and sD, sD(k) is new */
				sDvec[k] += D_ik;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				--N[k];
				sd += D_kj;
				++n;
			}
			sDvec[k] -= D_ik;
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sDvec[j] = sd;
}

void updateD_FF(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj) {
	
	unsigned k, n, Dn;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	
	/*
	Furthest First
	*/
	
	/* init */
	Dmat = D->mat;
	D_i = Dmat[i];
	D_j = Dmat[j];
	sDvec = sD->vec;
	
	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;
	
	/* update (first) row */
	for(k = 0; k < j; ++k) {
		D_ik = D_i[k];
		D_kj = D_j[k];
		if(0 <= D_ik && 0 <= D_kj) {
			dist = D_ik < D_kj ? D_kj : D_ik;
			D_j[k] = dist;
			/* update N and sD */
			sDvec[k] -= (D_kj - dist);
			sd += dist;
			--N[k];
			++n;
		} else if(0 <= D_ik) {
			D_j[k] = D_ik;
			/* update N and sD, sD(k) is new */
			sDvec[k] += D_ik;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			--N[k];
			sd += D_kj;
			++n;
		}
		sDvec[k] -= D_ik;
	}
	
	/* update jth column */
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			Dn = D->n;
		}
		while(++k < Dn) {
			D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
			D_kj = Dmat[k][j];
			if(0 <= D_ik && 0 <= D_kj) {
				dist = D_ik < D_kj ? D_kj : D_ik;
				Dmat[k][j] = dist;
				
				/* update N and sD */
				sDvec[k] -= (D_kj - dist);
				sd += dist;
				--N[k];
				++n;
			} else if(0 <= D_ik) {
				Dmat[k][j] = D_ik;
				/* update N and sD, sD(k) is new */
				sDvec[k] += D_ik;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				--N[k];
				sd += D_kj;
				++n;
			}
			sDvec[k] -= D_ik;
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sDvec[j] = sd;
}

void updateD_CF(Matrix *D, Vector *sD, unsigned *N, unsigned i, unsigned j, double Li, double Lj) {
	
	unsigned k, n, Dn;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	
	/*
	Closest First
	*/
	
	/* init */
	Dmat = D->mat;
	D_i = Dmat[i];
	D_j = Dmat[j];
	sDvec = sD->vec;
	
	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;
	
	/* update (first) row */
	for(k = 0; k < j; ++k) {
		D_ik = D_i[k];
		D_kj = D_j[k];
		if(0 <= D_ik && 0 <= D_kj) {
			dist = D_ik < D_kj ? D_ik : D_kj;
			D_j[k] = dist;
			/* update N and sD */
			sDvec[k] -= (D_kj - dist);
			sd += dist;
			--N[k];
			++n;
		} else if(0 <= D_ik) {
			D_j[k] = D_ik;
			/* update N and sD, sD(k) is new */
			sDvec[k] += D_ik;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			--N[k];
			sd += D_kj;
			++n;
		}
		sDvec[k] -= D_ik;
	}
	
	/* update jth column */
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			Dn = D->n;
		}
		while(++k < Dn) {
			D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
			D_kj = Dmat[k][j];
			if(0 <= D_ik && 0 <= D_kj) {
				dist = D_ik < D_kj ? D_ik : D_kj;
				Dmat[k][j] = dist < 0 ? 0 : dist;
				
				/* update N and sD */
				sDvec[k] -= (D_kj - dist);
				sd += dist;
				--N[k];
				++n;
			} else if(0 <= D_ik) {
				Dmat[k][j] = D_ik;
				/* update N and sD, sD(k) is new */
				sDvec[k] += D_ik;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				--N[k];
				sd += D_kj;
				++n;
			}
			sDvec[k] -= D_ik;
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sDvec[j] = sd;
}

unsigned * nj(Matrix *D, Vector *sD, unsigned *N, Qseqs **names) {
	
	unsigned i, j, mask, shift;
	long unsigned pair;
	double Li, Lj;
	Qseqs *tmp;
	
	/* init */
	N = initSummaD(sD, D, N);
	mask = UINT_MAX;
	shift = 8 * sizeof(unsigned);
	
	/* get pairs */
	while(D->n != 2 && (pair = minDist(D, sD, N))) {
		j = pair & mask;
		i = pair >> shift;
		
		/* get limbs */
		limbLengthPtr(&Li, &Lj, i, j, sD, N, D->mat[i][j]);
		
		/* form leaf */
		formNode(names[j], names[i], Lj, Li);
		
		/* update D and vectors */
		updateDptr(D, sD, N, i, j, Li, Lj);
		
		/* rearrange */
		ltdMatrix_popArrange(D, i);
		sD->vec[i] = sD->vec[--sD->n];
		N[i] = N[D->n];
		exchange(names[i], names[D->n], tmp);	
	}
	
	if(D->n == 2) {
		formLastNodePtr(*names, names[1], **(D->mat));
	} else {
		/* form remaining nodes with undefined distance */
		while(D->n != 1) {
			/* form leaf */
			formLastNodePtr(*names, names[--D->n], -1.0);
		}
	}
	
	/* verify valid newick format */
	if(*((*names)->seq) != '(') {
		byteshift((*names)->seq, (*names)->len, '(');
	}
	
	return N;
}

unsigned * nj_thread(Matrix *D, Vector *sD, unsigned *N, Qseqs **names, int thread_num) {
	
	unsigned i, j;
	double Li, Lj;
	NJthread *threads, *thread;
	Qseqs *tmp;
	
	/* here */
	if(0 && thread_num == 1) {
		return nj(D, sD, N, names);
	}
	
	/* init */
	N = initSummaD(sD, D, N);
	j = 0;
	
	/* start threads */
	threads = 0;
	thread = threads;
	i = thread_num;
	while(i) {
		thread = smalloc(sizeof(NJthread));
		thread->D = D;
		thread->sD = sD;
		thread->N = N;
		thread->min = 0;
		thread->mi = 0;
		thread->mj = 0;
		thread->num = thread_num;
		thread->next = threads;
		threads = thread;
		if(--i && (errno = pthread_create(&thread->id, NULL, minDist_thread, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", thread_num - i);
			i = 0;
		}
	}
	thread->id = 0;
	thread->min = 1;
	
	/* get pairs */
	while(D->n != 2 && (minDist_thread(thread) || ((i = thread->mi) | (j = thread->mj)))) {
		/* get limbs */
		limbLengthPtr(&Li, &Lj, i, j, sD, N, D->mat[i][j]);
		
		/* form leaf */
		formNode(names[j], names[i], Lj, Li);
		
		/* update D and vectors */
		updateDptr(D, sD, N, i, j, Li, Lj);
		
		/* rearrange */
		ltdMatrix_popArrange(D, i);
		sD->vec[i] = sD->vec[--sD->n];
		N[i] = N[D->n];
		exchange(names[i], names[D->n], tmp);
		
		threads->min = 1;
	}
	
	if(D->n == 2) {
		formLastNodePtr(*names, names[1], **(D->mat));
	} else {
		/* form remaining nodes with undefined distance */
		while(D->n != 1) {
			/* form leaf */
			formLastNodePtr(*names, names[--D->n], -1.0);
		}
	}
	
	/* verify valid newick format */
	if(*((*names)->seq) != '(') {
		byteshift((*names)->seq, (*names)->len, '(');
	}
	
	/* collect threads */
	minDist_thread(0);
	thread = threads->next;
	free(threads);
	while(thread) {
		threads = thread->next;
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
		free(thread);
		thread = threads;
	}
	
	return N;
}
