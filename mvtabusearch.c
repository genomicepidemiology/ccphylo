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
#include "pherror.h"
#include "mvtabusearch.h"
#define ABS(num)((num) < 0 ? -(num) : (num))

/* here */
/* almost from scratch */
double negotiateMVM(Machine *Mm, Machine *Mn, Job **Jmbest, Job **Jnbest) {
	
	double Mmj, Mnj, w1, w2, min, test, best, base;
	Job *Jm, *Jn, *next, *Jmin, *JmPrev, *JnPrev;
	
	/* exchange example:
	M1: eeddd
	M2: ccbbaaa
	 to
	M1: aaaddd
	M2: ccbbee
	*/
	
	/* no chance of improvement */
	if(Mm->avail == Mn->avail || (Mm->n <= 1 && Mn->n <= 1)) {
		return 0;
	}
	
	/* init minimum trade value */
	w1 = Mm->avail < 0 ? -Mm->avail : Mm->avail;
	w2 = Mn->avail < 0 ? -Mn->avail : Mn->avail;
	base = w1 + w2;
	best = base;
	min = best;
	*Jmbest = 0;
	*Jnbest = 0;
	
	/* identify best trade option by minimizing:
	(Mm->avail + Mm->job_i - Mn->job_j)^2 + (Mn->avail - Mm->job_i + Mn->job_j)^2
	(Mmj - Mn->job_j)^2 + (Mnj + Mn->job_j)^2
	*/
	Jm = Mm->jobs;
	JmPrev = 0;
	Jn = Mn->jobs;
	JnPrev = 0;
	while(Jm) {
		/* set constants */
		Mmj = Mm->avail + Jm->weight;
		Mnj = Mn->avail - Jm->weight;
		
		/* set first option as best */
		w1 = (Mmj - Jn->weight);
		w1 = w1 < 0 ? -w1 : w1;
		w2 = (Mnj + Jn->weight);
		w2 = w2 < 0 ? -w2 : w2;
		min = Jm->weight != Jn->weight ? w1 + w2 : base; /* avoid double approximations */
		Jmin = JnPrev; /* set pointer to prev job to allow post-merging */
		next = Jn->next;
		
		/* since both list of jobs are sorted the best trade option can be 
		identified by a merge search in O(m + n / m) */
		while(next) {
			if(Jm->weight != next->weight) { /* avoid double approximations */
				/* test next trade option */
				w1 = (Mmj - next->weight);
				w1 = w1 < 0 ? -w1 : w1;
				w2 = (Mnj + next->weight);
				w2 = w2 < 0 ? -w2 : w2;
				test = w1 + w2;
				
				/* advance Jn if value-error is not increasing */
				if(test <= min) {
					min = test;
					Jmin = Jn; /* set pointer to prev job to allow post-merging */
					JnPrev = Jn;
					Jn = next;
					next = next->next;
				} else {
					/* continuing will decrease trade value */
					next = 0;
				}
			} else {
				JnPrev = Jn;
				Jn = next;
				next = next->next;
			}
		}
		
		/* test best */
		if(min < best) {
			best = min;
			/* set pointer to best exchange jobs */
			*Jmbest = JmPrev;
			*Jnbest = Jmin;
		}
		JmPrev = Jm;
		Jm = Jm->next;
	}
	
	/* avoid double approximations */
	if(*Jmbest) {
		Jm = (*Jmbest)->next;
	} else {
		Jm = Mm->jobs;
	}
	if(*Jnbest) {
		Jn = (*Jnbest)->next;
	} else {
		Jn = Mn->jobs;
	}
	
	/* avoid double approximations */
	Mmj = Mm->avail;
	Mmj += (Jm->weight - Jn->weight);
	Mmj = Mmj < 0 ? -Mmj : Mmj;
	Mnj = Mn->avail;
	Mnj += (Jn->weight - Jm->weight);
	Mnj = Mnj < 0 ? -Mnj : Mnj;
	if(base != Mmj + Mnj && Jm->weight != Jn->weight) {
		best -= base;
	} else {
		best = 0;
	}
	
	/* return absolute reduction in error */
	return best;
}

double testMVhandover(Machine *Mm, Machine *Mn, Job *J) {
	
	int m;
	double post, prev, *mAvails, *nAvails, *jWeights;
	
	/* init */
	m = Mm->m;
	mAvails = Mm->Avails;
	nAvails = Mn->Avails;
	jWeights = J->Weights;
	
	/* get error prior and post trade */
	prev = 0;
	post = 0;
	++m;
	while(--m) {
		prev = ABS(*mAvails) + ABS(*nAvails);
		post = ABS(*mAvails + *jWeights) + ABS(*nAvails - *jWeights);
		++mAvails;
		++nAvails;
		++jWeights;
	}
	
	return prev - post;
}

void rmMVjob(Machine *M, Job *J) {
	
	int m;
	double *Avails, *Weights;
	
	/* init */
	m = M->m + 1;
	Avails = M->Avails - 1;
	Weights = J->Weights - 1;
	
	/* update Avails */
	while(--m) {
		*++Avails += *++Weights;
	}
}

void addMVjob(Machine *M, Job *J) {
	
	int m;
	double *Avails, *Weights;
	
	/* init */
	m = M->m + 1;
	Avails = M->Avails - 1;
	Weights = J->Weights - 1;
	
	/* update Avails */
	while(--m) {
		*++Avails -= *++Weights;
	}
}

int mvhandover(Machine *Mm, Machine *Mn) {
	
	int handovers;
	Job *J;
	Machine *middleman;
	
	if(Mn->avail < Mm->avail) {
		middleman = Mn;
		Mn = Mm;
		Mm = middleman;
	}
	
	/* handover jobs as long as there is an improvement */
	/*
		(mma + j)**2 + (mna - j)**2 < mma**2 + mna**2
		mma**2 + j**2 + 2 mma j + mna**2 + j**2 - 2 mna j < mma**2 + mna**2
		2 j**2 + 2 mma j - 2 mna j < 0
		2 j (j + mma - mna) < 0 ->
		j + mma - mna < 0
		mma < mna - j
	*/
	handovers = 0;
	J = Mm->jobs;
	//while(testHandover(Mm, Mn, J)) {
	while(J && Mm->avail < Mn->avail - J->weight) {
		/* check mv improvement before handover */
		if(0 < testMVhandover(Mm, Mn, J)) {
			/* handover job */
			Mm->n--;
			Mn->n++;
			Mm->avail += J->weight;
			Mn->avail -= J->weight;
			rmMVjob(Mm, J);
			addMVjob(Mn, J);
			Mm->jobs = J->next;
			J->next = 0;
			Mn->jobs = jobmerge_inc(Mn->jobs, J);
			++handovers;
		}
		J = Mm->jobs;
	}
	
	return handovers;
}
