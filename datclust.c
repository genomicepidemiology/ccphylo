/* Philip T.L.C. Clausen Jul 2021 plan@dtu.dk */

/*
 * Copyright (c) 2021, Philip Clausen, Technical University of Denmark
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

#include <float.h>
#include <limits.h>
#include "dat.h"
#include "datclust.h"
#include "distcmp.h"
#include "hclust.h"
#include "nwck.h"
#include "qseqs.h"
#include "vector.h"

void initQ_Dmat(Dat *Dmat, Vector *Q, int *P) {
	
	int m, n, i, j;
	double d, **Dptr, *Di, **Dj, *Qptr;
	float **Dfptr, *Dfi, **Dfj;
	short unsigned **Dsptr, *Dsi, **Dsj;
	unsigned char **Dbptr, *Dbi, **Dbj;
	
	/* init */
	m = Dmat->m;
	n = Dmat->n;
	Dptr = Dmat->mat;
	Dfptr = Dmat->fmat;
	Dsptr = Dmat->smat;
	Dbptr = Dmat->bmat;
	Q->n = m;
	Qptr = Q->vec;
	*Qptr = DBL_MAX;
	*P = -1;
	
	/* find min for row in D */
	i = 0;
	while(++i < m) {
		*++P = 0;
		j = 0;
		if(Dptr) {
			Di = Dptr[i];
			Dj = Dptr;
			*++Qptr = distcmp_d(Di, *Dj, n);
			while(++j < i) {
				if(0 <= (d = distcmp_d(Di, *++Dj, n)) && d <= *Qptr) {
					*Qptr = d;
					*P = j;
				}
			}
		} else if(Dfptr) {
			Dfi = Dfptr[i];
			Dfj = Dfptr;
			*++Qptr = distcmp_f(Dfi, *Dfj, n);
			while(++j < i) {
				if(0 <= (d = distcmp_f(Dfi, *++Dfj, n)) && d <= *Qptr) {
					*Qptr = d;
					*P = j;
				}
			}
		} else if(Dsptr) {
			Dsi = Dsptr[i];
			Dsj = Dsptr;
			*++Qptr = distcmp_s(Dsi, *Dsj, n);
			while(++j < i) {
				if(0 <= (d = distcmp_s(Dsi, *++Dsj, n)) && d <= *Qptr) {
					*Qptr = d;
					*P = j;
				}
			}
		} else {
			Dbi = Dbptr[i];
			Dbj = Dbptr;
			*++Qptr = distcmp_b(Dbi, *Dbj, n);
			while(++j < i) {
				if(0 <= (d = distcmp_b(Dbi, *++Dbj, n)) && d <= *Qptr) {
					*Qptr = d;
					*P = j;
				}
			}
		}
	}
}

void updateQP(Vector *Q, int *P, int m, int n) {
	
	int i;
	
	Q->vec[m] = DBL_MAX;
	P += m;
	*P = -1;
	i = Q->n;
	while(--i != m) {
		if(*++P == m) {
			*P = n;
		}
	}
}

long unsigned pairU(int *P, int n) {
	
	int i;
	long unsigned pos;
	
	pos = 0;
	i = 0;
	while(++i < n) {
		if(*++P != -1) {
			if(!pos) {
				pos = i;
			} else {
				pos <<= sizeof(unsigned) * 8;
				pos |= i;
				return pos;
			}
		}
	}
	
	return 0;
}

void tclust(Vector *Q, int *P, Qseqs **names) {
	
	int i, j, n, mask, shift;
	long unsigned pair;
	double limblength;
	Qseqs *tmp;
	
	/* init */
	j = 0;
	mask = UINT_MAX;
	shift = 8 * sizeof(unsigned);
	
	/* get pairs */
	n = Q->n;
	while(n != 1 && (pair = pairQ(Q, P))) {
		j = pair & mask;
		i = pair >> shift;
		
		/* get limblength */
		limblength = Q->vec[i] / 2;
		
		/* form leaf */
		formNode(names[j], names[i], limblength, limblength);
		
		/* update Q and P */
		updateQP(Q, P, i, j);
		--n;
	}
	
	if(n != 1) {
		/* form remaining nodes with undefined distance */
		while(n != 1 && (pair = pairU(P, Q->n))) {
			j = pair & mask;
			i = pair >> shift;
			
			/* form leaf */
			formLastNodePtr(names[j], names[i], -1.0);
			P[i] = -1;
			--n;
		}
	}
	exchange(*names, names[j], tmp);
	
}
