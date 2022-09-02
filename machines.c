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

#include "jobs.h"
#include "machines.h"
#include "pherror.h"

Machine * machinemerge(Machine *L1, Machine *L2) {
	
	Machine *dest, *ptr;
	
	/* set dest to higher value */
	if(!L1) {
		return L2;
	} else if(!L2) {
		return L1;
	} else if(L1->avail < L2->avail) {
		dest = L2;
		L2 = L2->next;
	} else {
		dest = L1;
		L1 = L1->next;
	}
	ptr = dest;
	
	/* merge lists */
	while(L1 && L2) {
		if(L2->avail < L1->avail) {
			ptr->next = L1;
			L1 = L1->next;
		} else {
			ptr->next = L2;
			L2 = L2->next;
		}
		ptr = ptr->next;
	}
	
	/* link remaining */
	ptr->next = L1 ? L1 : L2;
	
	return dest;
}

Machine * machinesort(Machine *src, int m) {
	
	int mid;
	Machine *L1, *L2;
	
	if(m <= 1) {
		/* error if n is < 0 */
		if(m == 1) {
			src->next = 0;
		} else {
			src = 0;
		}
		return src;
	} else {
		/* split, sort and merge */
		mid = m >> 1;
		L1 = machinesort(src, mid);
		L2 = machinesort(src + mid, m - mid);
		return machinemerge(L1, L2);
	}
	
	return 0;
}

Machine * initM(int m, int n, Job *J) {
	
	double m_target;
	Machine *dest, *ptr;
	
	dest = smalloc(m * sizeof(Machine));
	m_target = totM(J, n) / m;
	
	ptr = dest - 1;
	++m;
	while(--m) {
		(++ptr)->num = m;
		ptr->n = 0;
		ptr->avail = m_target;
		ptr->jobs = 0;
		ptr->next = ptr + 1;
	}
	ptr->next = 0;
	
	return dest;
}

Machine * initSkewM(int m, int n, Job *J, double *loads) {
	
	int i;
	double totL, m_target;
	Machine *dest, *ptr;
	
	i = m;
	totL = *loads;
	while(--i) {
		totL += *++loads;
	}
	loads -= m;
	
	dest = smalloc(m * sizeof(Machine));
	m_target = totM(J, n) / totL;
	
	ptr = dest - 1;
	i = m + 1;
	while(--i) {
		(++ptr)->num = i;
		ptr->n = 0;
		ptr->avail = m_target * *++loads;
		ptr->jobs = 0;
		ptr->next = ptr + 1;
	}
	ptr->next = 0;
	
	return dest;
}

double machineMSE(Machine *M) {
	
	int m;
	double mse;
	
	/* calculate mean square error of makespan */
	m = 1;
	mse = M->avail * M->avail;
	while((M = M->next)) {
		mse += M->avail * M->avail;
		++m;
	}
	
	return mse / m;
}

void print_stats(Machine *M) {
	
	int m;
	double mse, Cmax, Cmin, L1, OPT, Jmax;
	Job *J;
	
	m = 0;
	mse = 0;
	Cmax = M->avail;
	Cmin = M->avail;
	L1 = 0;
	Jmax = M->jobs->weight;
	OPT = 0;
	while(M) {
		/* Update Machine stats */
		if(Cmax < M->avail) {
			Cmax = M->avail;
		} else if(M->avail < Cmin) {
			Cmin = M->avail;
		}
		L1 += M->avail < 0 ? -M->avail : M->avail;
		mse += M->avail * M->avail;
		++m;
		
		/* Update OPT */
		J = M->jobs;
		while(J) {
			OPT += J->weight;
			if(Jmax < J->weight) {
				Jmax = J->weight;
			}
			J = J->next;
		}
		M = M->next;
	}
	mse /= m;
	OPT /= m;
	Cmax += OPT;
	Cmin += OPT;
	OPT = OPT < Jmax ? Jmax : OPT;
	
	fprintf(stderr, "## MSE:\t%f\n", mse);
	fprintf(stderr, "## L1:\t%f\n", L1);
	fprintf(stderr, "## OPT:\t%f\n", OPT);
	fprintf(stderr, "## Cmax:\t%f\n", Cmax);
	fprintf(stderr, "## Cmin:\t%f\n", Cmin);
}
