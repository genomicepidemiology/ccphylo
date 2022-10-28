/* Philip T.L.C. Clausen Oct 2022 plan@dtu.dk */

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
#include "mvmakespan.h"

void addMVDBE(Machine **Mdest, Machine **Edest, Job *J, int m, int n) {
	
	Machine *M, *E, *nextM;
	
	/* init */
	M = *Mdest;
	E = *Edest;
	
	/* put job on least loaded machine */
	M->n++;
	J->next = M->jobs;
	M->jobs = J;
	M->avail -= J->weight;
	
	/* move machine down the qeue */
	nextM = M->next;
	M->next = 0;
	if(M->n < n / m) {
		M = machinemerge(nextM, M);
	} else {
		/* remove machine from available list */
		E = machinemerge(E, M);
		M = nextM;
	}
	
	/* set pointers */
	*Mdest = M;
	*Edest = E;
}

Machine * addMVDBF(Machine *M, Job *J) {
	
	Machine *nextM;
	
	/* put job on least loaded machine */
	M->n++;
	J->next = M->jobs;
	M->jobs = J;
	M->avail -= J->weight;
	
	/* move machine down the qeue */
	nextM = M->next;
	M->next = 0;
	return machinemerge(nextM, M);
}

Machine * MVFirstFit(Machine *M, Job *J, int m) {
	
	double weight, best;
	Machine *F;
	
	weight = J->weight;
	best = M->avail;
	F = M;
	while(m) {
		/* here */
		/* test fit */
		if(weight <= M->avail) {
			/* job fits */
			M->n++;
			J->next = M->jobs;
			M->jobs = J;
			M->avail -= weight;
			return M;
		} else if(best < M->avail) {
			/* machine is less filled than the previous */
			best = M->avail;
			F = M;
		}
		
		/* go to next machine */
		M = M->next;
		--m;
	}
	
	/* job does not fit anywhere */
	F->n++;
	J->next = F->jobs;
	F->jobs = J;
	F->avail -= weight;
	
	return F;
}

Machine * MVFirstFet(Machine *M, Job *J) {
	
	double weight, best;
	Machine *F, *prev, *prevF;
	
	weight = J->weight;
	best = M->avail;
	F = M;
	prev = 0;
	prevF = 0;
	while(M) {
		/* here */
		/* test fit */
		if(weight <= M->avail) {
			/* job fits */
			M->n++;
			J->next = M->jobs;
			M->jobs = J;
			M->avail -= weight;
			return prev;
		} else if(best < M->avail) {
			/* machine is less filled than the previous */
			best = M->avail;
			prevF = prev;
			F = M;
		}
		
		/* go to next machine */
		prev = M;
		M = M->next;
	}
	
	/* job does not fit anywhere */
	F->n++;
	J->next = F->jobs;
	F->jobs = J;
	F->avail -= weight;
	
	return prevF;
}
