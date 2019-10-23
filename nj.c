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

#include <limits.h>
#include <stdlib.h>
#include "matrix.h"
#include "nj.h"
#include "nwck.h"
#include "qseqs.h"
#include "pherror.h"
#include "vector.h"

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
		/* j only shares similarity with i */
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

long unsigned initQ(Matrix *D, Matrix *Q, Vector *sD, unsigned *N) {
	
	int i, j, mi, mj;
	long unsigned pos;
	double dist, min, *Dptr, *Qptr, *sDvec;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* alloc */
	if(Q->size < (Q->n = D->n)) {
		ltdMatrix_realloc(Q, Q->n);
	}
	
	mi = 0;
	mj = 0;
	min = 1;
	Dptr = *(D->mat) - 1;
	Qptr = *(Q->mat) - 1;
	sDvec = sD->vec;
	for(i = 1; i < Q->n; ++i) {
		for(j = 0; j < i; j++) {
			if(0 <= (dist = *++Dptr)) {
				dist *= (N[i] + N[j] - 4) / 2;
				if((*++Qptr = dist - sDvec[i] - sDvec[j]) < min) {
					min = *Qptr;
					mi = i;
					mj = j;
				}
			} else {
				/* missing */
				*++Qptr = 1;
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
				Dmat[k][j] = dist;
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

unsigned * nj(Matrix *D, Matrix *Q, Vector *sD, unsigned *N, Qseqs **names) {
	
	unsigned i, j, mask, shift;
	long unsigned pair;
	double Li, Lj;
	Qseqs *tmp;
	
	/* init */
	N = initSummaD(sD, D, N);
	mask = UINT_MAX;
	shift = 8 * sizeof(unsigned);
	
	/* get pairs */
	while(D->n != 2 && (pair = initQ(D, Q, sD, N))) {
		j = pair & mask;
		i = pair >> shift;
		
		/* get limbs */
		limbLength(&Li, &Lj, i, j, sD, N, D->mat[i][j]);
		
		/* form leaf */
		formNode(names[j], names[i], Lj, Li);
		
		/* update D and vectors */
		updateD(D, sD, N, i, j, Li, Lj);
		
		/* rearrange */
		ltdMatrix_popArrange(D, i);
		sD->vec[i] = sD->vec[--sD->n];
		N[i] = N[D->n];
		tmp = names[i];
		names[i] = names[D->n];
		names[D->n] = tmp;
	}
	
	if(D->n == 2) {
		formLastNode(*names, names[1], **(D->mat));
	} else {
		/* form remaining nodes with undefined distance */
		while(D->n != 1) {
			/* form leaf */
			formLastNode(names[0], names[--D->n], -1.0);
		}
	}
	return N;
}
