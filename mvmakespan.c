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
#include "mvjobs.h"
#include "mvmakespan.h"

void addMVDBE(Machine **Mdest, Machine **Edest, Job *J, int m, int n) {
	
	double max, test;
	Machine *M, *E, *B, *prev, *prevB;
	
	/* init */
	M = *Mdest;
	E = *Edest;
	B = M;
	prev = 0;
	prevB = 0;
	
	/* find best machine */
	max = M->avail < 0 ? M->avail - J->weight : -M->avail - J->weight;
	while(M) {
		test = addValue(M, J);
		if(max < test) {
			max = test;
			prevB = prev;
			B = M;
			if(max == J->weight) {
				break;
			}
		}
		
		prev = M;
		M = M->next;
	}
	
	/* put job on least loaded machine */
	addMVjobToMachine(B, J);
	
	/* isolate best machine */
	if(prevB) {
		M = prevB;
		M->next = B->next;
		M = *Mdest;
	} else {
		M = B->next;
	}
	B->next = 0;
	
	/* move machine down the qeue */
	if(B->n < n / m) {
		M = machinemerge(M, B);
	} else {
		/* remove machine from available list */
		E = machinemerge(E, B);
	}
	
	/* set pointers */
	*Mdest = M;
	*Edest = E;
}

Machine * addMVDBF(Machine *M, Job *J) {
	
	double max, test;
	Machine *B, *Mptr, *prev, *prevB;
	
	/* init */
	Mptr = M;
	B = M;
	prev = 0;
	prevB = 0;
	
	/* find best machine */
	max = M->avail < 0 ? M->avail - J->weight : -M->avail - J->weight;
	while(Mptr) {
		test = addValue(Mptr, J);
		if(max < test) {
			max = test;
			prevB = prev;
			B = Mptr;
			if(max == J->weight) {
				break;
			}
		}
		
		prev = Mptr;
		Mptr = Mptr->next;
	}
	
	/* put job on least loaded machine */
	addMVjobToMachine(B, J);
	
	/* isolate best machine */
	if(prevB) {
		Mptr = prevB;
		Mptr->next = B->next;
	} else {
		M = B->next;
	}
	B->next = 0;
	
	/* move machine down the qeue */
	return machinemerge(M, B);
}

Machine * MVFirstFit(Machine *M, Job *J, int m) {
	
	double weight, best, test;
	Machine *F;
	
	weight = J->weight;
	best = M->avail < 0 ? M->avail - weight : -M->avail - weight;
	F = M;
	while(m) {
		test = addValue(M, J);
		/* test fit */
		if(test == weight) {
			/* job fits */
			addMVjobToMachine(M, J);
			return M;
		} else if(best < test) {
			/* machine is less filled than the previous */
			best = test;
			F = M;
		}
		
		/* go to next machine */
		M = M->next;
		--m;
	}
	
	/* job does not fit anywhere */
	addMVjobToMachine(F, J);
	
	return F;
}

Machine * MVFirstFet(Machine *M, Job *J) {
	
	double weight, best, test;
	Machine *F, *prev, *prevF;
	
	weight = J->weight;
	best = M->avail < 0 ? M->avail - weight : -M->avail - weight;
	F = M;
	prev = 0;
	prevF = 0;
	while(M) {
		test = addValue(M, J);
		/* test fit */
		if(test == weight) {
			/* job fits */
			addMVjobToMachine(M, J);
			return prev;
		} else if(best < test) {
			/* machine is less filled than the previous */
			best = test;
			prevF = prev;
			F = M;
		}
		
		/* go to next machine */
		prev = M;
		M = M->next;
	}
	
	/* job does not fit anywhere */
	addMVjobToMachine(F, J);
	
	return prevF;
}
