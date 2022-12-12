/* Philip T.L.C. Clausen Jul 2022 plan@dtu.dk */

/*
 * Copyright (c) 2022, Philip Clausen, Technical University of Denmark
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "filebuff.h"
#include "jobs.h"
#include "pherror.h"

void (*jobWeight)(Job *src, int n, double logbase) = &nullWeight;

Job * job_realloc(Job *src, int mv, int oldsize, int newsize) {
	
	int i;
	Job *ptr;
	
	/* realloc src */
	ptr = realloc(src, newsize * sizeof(Job));
	if(!ptr) {
		ERROR();
	}
	
	/* update links */
	if(ptr != src) {
		src = ptr;
		i = newsize < oldsize ? newsize : oldsize;
		++i;
		while(--i) {
			ptr->next = ptr + 1;
			++ptr;
		}
	} else {
		src = ptr;
	}
	
	/* init new part */
	if(oldsize < newsize) {
		ptr = src + oldsize;
		newsize -= --oldsize;
		while(--newsize) {
			ptr->num = ++oldsize;
			ptr->size = 0;
			ptr->weight = 0;
			if(mv) {
				ptr->Weights = calloc(mv, sizeof(double));
				if(!ptr->Weights) {
					ERROR();
				}
			} else {
				ptr->Weights = 0;
			}
			ptr->next = ptr + 1;
			++ptr;
		}
	}
	
	return src;
}

Job * jobWeights_realloc(Job *src, int mv, int mv_new, int n) {
	
	int i;
	double *Weights;
	Job *ptr;
	
	if(mv == mv_new) {
		return src;
	}
	
	/* init */
	--mv;
	ptr = src;
	++n;
	while(--n) {
		/* realloc weights */
		if(ptr->Weights) {
			ptr->Weights = realloc(ptr->Weights, mv_new * sizeof(double));
		} else {
			ptr->Weights = malloc(mv_new * sizeof(double));
		}
		if(!ptr->Weights) {
			ERROR();
		}
		
		if(mv < mv_new) {
			/* init new part */
			Weights = ptr->Weights + mv;
			i = mv_new - mv;
			while(--i) {
				*++Weights = 0;
			}
		}
	}
	
	return src;
}

Job * jobmerge(Job *L1, Job *L2) {
	
	Job *dest, *ptr;
	
	/* set dest to higher value */
	if(!L1) {
		return L2;
	} else if(!L2) {
		return L1;
	} else if(L1->weight < L2->weight) {
		dest = L2;
		L2 = L2->next;
	} else {
		dest = L1;
		L1 = L1->next;
	}
	ptr = dest;
	
	/* merge lists */
	while(L1 && L2) {
		if(L1->weight < L2->weight) {
			ptr->next = L2;
			L2 = L2->next;
		} else {
			ptr->next = L1;
			L1 = L1->next;
		}
		ptr = ptr->next;
	}
	
	/* link remaining */
	ptr->next = L1 ? L1 : L2;
	
	return dest;
}

Job * jobmerge_inc(Job *L1, Job *L2) {
	
	Job *dest, *ptr;
	
	/* set dest to higher value */
	if(!L1) {
		return L2;
	} else if(!L2) {
		return L1;
	} else if(L2->weight < L1->weight) {
		dest = L2;
		L2 = L2->next;
	} else {
		dest = L1;
		L1 = L1->next;
	}
	ptr = dest;
	
	/* merge lists */
	while(L1 && L2) {
		if(L2->weight < L1->weight) {
			ptr->next = L2;
			L2 = L2->next;
		} else {
			ptr->next = L1;
			L1 = L1->next;
		}
		ptr = ptr->next;
	}
	
	/* link remaining */
	ptr->next = L1 ? L1 : L2;
	
	return dest;
}

Job * jobsort(Job *src, int n) {
	
	int mid;
	Job *L1, *L2;
	
	if(n <= 1) {
		/* error if n is < 0 */
		if(n == 1) {
			src->next = 0;
		} else {
			src = 0;
		}
		return src;
	} else {
		/* split, sort and merge */
		mid = n >> 1;
		L1 = jobsort(src, mid);
		L2 = jobsort(src + mid, n - mid);
		return jobmerge(L1, L2);
	}
	
	return 0;
}

double totM(Job *J, int n) {
	
	double sum;
	
	sum = J->weight;
	while(--n) {
		sum += (++J)->weight;
	}
	
	return sum;
}

double * totMVM(Job *J, int n, int mv) {
	
	int i;
	double *targets, *targetsPtr, *targetsJ;
	
	if(!mv) {
		return 0;
	}
	
	targets = calloc(mv, sizeof(double));
	if(!targets) {
		ERROR();
	}
	++n;
	--J;
	while(--n) {
		targetsPtr = targets - 1;
		targetsJ = (++J)->Weights - 1;
		i = mv + 1;
		while(--i) {
			*++targetsPtr += *++targetsJ;
		}
	}
	
	return targets;
}

double optM(Job *J, int n, int m) {
	
	double sum, max;
	
	sum = J->weight;
	max = sum;
	while(--n) {
		sum += (++J)->weight;
		if(max < J->weight) {
			max = J->weight;
		}
	}
	
	sum /= m;
	
	return sum < max ? sum : max;
}

double targetM(Job *J, int n, int m) {
	
	double sum, max;
	
	sum = J->weight;
	max = sum;
	while(--n) {
		sum += (++J)->weight;
		if(max < J->weight) {
			max = J->weight;
		}
	}
	
	if(sum / m < max) {
		sum -= max;
		--m;
	}
	sum /= m;
	
	return sum;
}

void nullWeight(Job *src, int n, double logbase) {
	
	Job *dest;
	
	dest = src;
	++n;
	while(--n) {
		dest->weight = dest->size;
		++dest;
	}

}

void logWeight(Job *src, int n, double logbase) {
	
	Job *dest;
	
	if(!logbase || !(logbase = log(logbase))) {
		fprintf(stderr, "Invalid logbase\n");
		exit(1);
	}
	
	dest = src - 1;
	++n;
	while(--n) {
		if((++dest)->size) {
			dest->weight = 1 + log(dest->size) / logbase;
		} else {
			fprintf(stderr, "Invalid weight for log-transformation:\t%d\n", dest->size);
			exit(1);
		}
	}
}

void polWeight(Job *src, int n, double exponent) {
	
	Job *dest;
	
	dest = src;
	++n;
	while(--n) {
		dest->weight = pow(dest->size, exponent);
		++dest;
	}
}

void expWeight(Job *src, int n, double expobase) {
	
	Job *dest;
	
	dest = src;
	++n;
	while(--n) {
		dest->weight = pow(expobase, dest->size);
		++dest;
	}
}

int cleanJobs(Job *src, int n) {
	
	int i;
	Job *cur, *next;
	
	i = n + 1;
	cur = src;
	next = src;
	while(--i) {
		if(0 < cur->size) {
			*next = *cur;
			next->next = next + 1;
			++next;
		} else {
			free(cur->Weights);
			--n;
		}
		++cur;
	}
	(--next)->next = 0;
	
	return n;
}

int cmpJ(Job *Jm, Job *Jn, int m) {
	
	double *mWeigths, *nWeigths;
	
	if(Jm->weight != Jn->weight) {
		return Jm->weight < Jn->weight ? 1 : -1;
	} else if(m++) {
		mWeigths = Jm->Weights - 1;
		nWeigths = Jn->Weights - 1;
		while(--m) {
			if(*++mWeigths != *++nWeigths) {
				return *mWeigths < *nWeigths ? 1 : -1;
			}
		}
	}
	
	return 0;
}
