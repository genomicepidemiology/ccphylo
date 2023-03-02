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

Machine * initM(int m, int n, int mv, Job *J) {
	
	int i;
	double m_target, *m_targets, *Avails, *targetPtr;
	Machine *dest, *ptr;
	
	dest = smalloc(m * sizeof(Machine));
	m_target = totM(J, n) / m;
	m_targets = totMVM(J, n, mv);
	i = mv + 1;
	targetPtr = m_targets;
	while(--i) {
		*targetPtr = *targetPtr / m;
		++targetPtr;
	}
	
	ptr = dest - 1;
	++m;
	while(--m) {
		(++ptr)->num = m;
		ptr->n = 0;
		ptr->m = mv;
		ptr->buff = 0;
		ptr->avail = m_target;
		if(m_targets) {
			i = mv + 1;
			targetPtr = m_targets - 1;
			ptr->Avails = smalloc(mv * sizeof(double));
			Avails = ptr->Avails - 1;
			while(--i) {
				*++Avails = *++targetPtr;
			}
		} else {
			ptr->Avails = 0;
		}
		ptr->jobs = 0;
		ptr->next = ptr + 1;
	}
	ptr->next = 0;
	
	return dest;
}

Machine * initSkewM(int m, int n, int mv, Job *J, double *loads) {
	
	int i;
	double totL, m_target, *m_targets, *Avails, *targetPtr;
	Machine *dest, *ptr;
	
	i = m;
	totL = *loads;
	while(--i) {
		totL += *++loads;
	}
	loads -= m;
	
	dest = smalloc(m * sizeof(Machine));
	m_target = totM(J, n) / totL;
	m_targets = totMVM(J, n, mv);
	
	ptr = dest - 1;
	++m;
	while(--m) {
		(++ptr)->num = m;
		ptr->n = 0;
		ptr->m = mv;
		ptr->buff = 0;
		ptr->avail = m_target * *++loads;
		if(m_targets) {
			i = mv + 1;
			targetPtr = m_targets - 1;
			ptr->Avails = smalloc(mv * sizeof(double));
			Avails = ptr->Avails - 1;
			while(--i) {
				*++Avails = *++targetPtr * *loads / totL;
			}
		} else {
			ptr->Avails = 0;
		}
		ptr->jobs = 0;
		ptr->next = ptr + 1;
	}
	ptr->next = 0;
	free(m_targets);
	
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

double machineIMSE(Machine *M) {
	
	int i, m;
	double imse, *avails;
	
	/* calculate imbalance mean square error of makespan */
	m = 0;
	imse = 0;
	while(M) {
		i = M->m + 1;
		avails = M->Avails;
		while(--i) {
			imse += *avails * *avails;
			++avails;
		}
		++m;
		M = M->next;
	}
	
	return imse / m;
}

void print_stats(Machine *M) {
	
	int i, m;
	double mse, imse, Cmax, Cmin, L1, L1imse, OPT, Jmax, *weights;
	Job *J;
	
	m = 0;
	mse = 0;
	imse = 0;
	Cmax = M->avail;
	Cmin = M->avail;
	L1 = 0;
	L1imse = 0;
	Jmax = M->jobs ? M->jobs->weight : 0;
	weights = 0;
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
		
		/* imbalance error */
		i = M->m + 1;
		weights = M->Avails;
		while(--i) {
			imse += *weights * *weights;
			L1imse += *weights < 0 ? -*weights : *weights;
			++weights;
		}
		
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
	imse /= m;
	OPT /= m;
	Cmax += OPT;
	Cmin += OPT;
	OPT = OPT < Jmax ? Jmax : OPT;
	
	fprintf(stderr, "## MSE:\t%f\n", mse);
	if(weights) {
		fprintf(stderr, "## Imbalance MSE:\t%f\n", imse);
	}
	fprintf(stderr, "## L1:\t%f\n", L1);
	if(weights) {
		fprintf(stderr, "## Imbalance L1:\t%f\n", L1imse);
	}
	fprintf(stderr, "## OPT:\t%f\n", OPT);
	fprintf(stderr, "## Cmax:\t%f\n", Cmax);
	fprintf(stderr, "## Cmin:\t%f\n", Cmin);
}
