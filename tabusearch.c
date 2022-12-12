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
#include "mvjobs.h"
#include "mvtabusearch.h"
#include "pherror.h"
#include "tabusearch.h"
#define ABS(num)((num) < 0 ? -(num) : (num))

int (*tradeM)(Machine *) = &tradeBB;
double (*negotiatePtr)(Machine *, Machine *, Job **, Job **) = &negotiateM;
int (*handoverPtr)(Machine *, Machine *) = &handover;

Job ** sequenceJobs(Machine *M, Job *J, int m, int n) {
	
	Job **sJ, **sJptr, *Jptr;
	
	/* allocate sequential jobs */
	sJ = smalloc((m + n + 1) * sizeof(Job *));
	
	/* link jobs to sequential order */
	sJptr = sJ;
	*sJptr = 0;
	while(M) {
		Jptr = M->jobs;
		while(Jptr) {
			*++sJptr = Jptr;
			Jptr = Jptr->next;
		}
		*++sJptr = 0;
		M = M->next;
	}
	
	return sJ;
}

void moveUp(Job **sJ) {
	
	Job **next, *carry;
	
	next = sJ + 1;
	while(*next && (*sJ)->weight < (*next)->weight) {
		carry = *sJ;
		*sJ = *next;
		*++sJ = carry;
		++next;
	}
}

void moveDown(Job **sJ) {
	
	Job **next, *carry;
	
	next = --sJ;
	while(*sJ && (*sJ)->weight < (*next)->weight) {
		carry = *sJ;
		*sJ = *next;
		*--sJ = carry;
		--next;
	}
}

void jobexchange(Job **sJ, unsigned m, unsigned n) {
	
	int num, size;
	double weight, *Weights;
	Job *Jm, *Jn;
	
	Jm = sJ[m];
	Jn = sJ[n];
	
	/* exchange jobs */
	num = Jm->num;
	size = Jm->size;
	weight = Jm->weight;
	Weights = Jm->Weights;
	
	Jm->num = Jn->num;
	Jm->size = Jn->size;
	Jm->weight = Jn->weight;
	Jm->Weights = Jn->Weights;
	
	Jn->num = num;
	Jn->size = size;
	Jn->weight = weight;
	Jn->Weights = Weights;
	
	/* rearrange sequential jobs */
	if(Jm->weight < Jn->weight) {
		moveUp(sJ + m);
		moveDown(sJ + n);
	} else {
		moveUp(sJ + n);
		moveDown(sJ + m);
	}
}

void insertJob(Machine *M, Job *J) {
	
	Job *Ji;
	
	/* get job to insert */
	Ji = J->next;
	
	/* check if repositioning is necessary */
	if(J->weight != Ji->weight) {
		/* extract inserion job */
		J->next = Ji->next;
		Ji->next = 0;
		if(J->weight < Ji->weight) {
			/* insert from J */
			J->next = jobmerge(J->next, Ji);
		} else {
			/* insert from top */
			M->jobs = jobmerge(M->jobs, Ji);
		}
	}
}

int exchangeJobs(Machine *Mm, Machine *Mn, Job *Jm, Job *Jn) {
	
	Job *J;
	
	/* isolate jobs for exchange */
	if(Jm) {
		J = Jm->next;
		Jm->next = J->next;
	} else {
		J = Mm->jobs;
		Mm->jobs = J->next;
	}
	J->next = 0;
	Jm = J;
	if(Jn) {
		J = Jn->next;
		Jn->next = J->next;
	} else {
		J = Mn->jobs;
		Mn->jobs = J->next;
	}
	J->next = 0;
	Jn = J;
	
	/* insert jobs in new machines */
	Mn->jobs = jobmerge_inc(Mn->jobs, Jm);
	Mm->jobs = jobmerge_inc(Mm->jobs, Jn);
	
	/* adjust availability of machines */
	Mm->avail += Jm->weight - Jn->weight;
	Mn->avail += Jn->weight - Jm->weight;
	
	/* adjust class availability of machines */
	rmMVjob(Mm, Jm);
	addMVjob(Mm, Jn);
	rmMVjob(Mn, Jn);
	addMVjob(Mn, Jm);
	
	/* cmp jobs to validate exchange */
	return cmpJ(Jm, Jn, Mm->m);
}

double negotiateM(Machine *Mm, Machine *Mn, Job **Jmbest, Job **Jnbest) {
	
	int balance;
	double Mmj, Mnj, Jmw, w1, w2, min, test, best, base;
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
	
	/* set balance of machines */
	balance = (Mm->avail < 0 && 0 < Mn->avail) || (Mn->avail < 0 && 0 < Mm->avail);
	
	/* init minimum trade value */
	if(balance) {
		base = ABS(Mm->avail) + ABS(Mn->avail);
	} else {
		w1 = ABS(Mm->avail);
		w2 = ABS(Mn->avail);
		base = w1 < w2 ? w2 : w1;
	}
	best = base;
	min = best;
	*Jmbest = 0;
	*Jnbest = 0;
	
	/* identify best trade option by minimizing:
	(Mm->avail + Mm->job_i - Mn->job_j)^2 + (Mn->avail - Mm->job_i + Mn->job_j)^2
	(Mmj - Mn->job_j)^2 + (Mnj + Mn->job_j)^2 // Makes doulbe approximation
	*/
	Jm = Mm->jobs;
	JmPrev = 0;
	Jn = Mn->jobs;
	JnPrev = 0;
	while(Jm) {
		/* set constant */
		Jmw = Jm->weight;
		Mmj = Mm->avail + Jmw;
		Mnj = Mn->avail;
		
		/* set first option as best */
		w1 = Mmj - Jn->weight;
		w2 = Mnj + Jn->weight - Jmw;
		if(balance) {
			min = ABS(w1) + ABS(w2);
		} else {
			w1 = ABS(w1);
			w2 = ABS(w2);
			min = (w1 < w2 ? w2 : w1);
		}
		Jmin = JnPrev; /* set pointer to prev job to allow post-merging */
		next = Jn->next;
		
		/* since both list of jobs are sorted the best trade option can be 
		identified by a merge search in O(m + n / m) */
		while(next) {
			if(Jm->weight != next->weight) {
				/* test next trade option */
				w1 = Mmj - next->weight;
				w2 = Mnj + next->weight - Jmw;
				if(balance) {
					test = ABS(w1) + ABS(w2);
				} else {
					w1 = ABS(w1);
					w2 = ABS(w2);
					test = (w1 < w2 ? w2 : w1);
				}
				
				/* advance Jn if value-error is not increasing */
				if(test < min) {
					min = test;
					Jmin = Jn; /* set pointer to prev job to allow post-merging */
					JnPrev = Jn;
					Jn = next;
					next = next->next;
				} else if(test == min) {
					JnPrev = Jn;
					Jn = next;
					next = next->next;
				} else {
					/* continuing will decrease trade value */
					next = 0;
				}
				if(min == 0) {
					/* continuing will not increase trade value */
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
		Jm = best == 0 ? 0 : Jm->next;
	}
	
	/* set pointers to best jobs */
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
	
	/* test validity of trade */
	if(best != base && Jm->weight != Jn->weight) {
		best -= base;
	} else {
		best = 0;
	}
	
	/* return absolute reduction in error */
	return best;
}

int tradeDBEB(Machine *M) {
	
	/* Der Bor En Bager */
	int trades, nullTrades;
	double min, test;
	Job *Jn, *Jm, *JnBest, *JmBest;
	Machine *Mm, *Mn, *Mbest;
	
	/* check MSE */
	test = M->m ? machineIMSE(M) : machineMSE(M);
	fprintf(stderr, "## Pre-tabu MSE:\t%f\n", test);
	if(test == 0) {
		return 0;
	}
	
	/* trade as long as there is an improvement */
	trades = 0;
	do {
		Mbest = 0;
		nullTrades = trades;
		
		/* go through machines */
		Mm = M;
		while(Mm) {
			/* Negotiate trade options for Mm */
			min = 0;
			Jm = 0;
			Jn = 0;
			JmBest = 0;
			JnBest = 0;
			Mn = Mm->next;
			while(Mn) {
				/* Negotiate trade options between Mm and Mn */
				test = negotiatePtr(Mm, Mn, &Jm, &Jn);
				
				/* check negotiation between machines */
				if(test < min) {
					min = test;
					JmBest = Jm;
					JnBest = Jn;
					Mbest = Mn;
				}
				
				Mn = Mn->next;
			}
			
			/* exchange jobs */
			if(min < 0 && exchangeJobs(Mm, Mbest, JmBest, JnBest)) {
				++trades;
			} else {
				Mm = Mm->next;
			}
		}
	} while(nullTrades != trades);
	
	return trades;
}

int testHandover(Machine *Mm, Machine *Mn, Job *J) {
	
	double e, tmp;
	
	/* get error after trade */
	if(Mn->avail < Mm->avail) { /* Mm->avail < Mn->avail to make a valid trade */
		e = Mn->avail - Mm->avail;
	} else if(Mm->avail < 0 && 0 < Mn->avail) { /* machines are somewhat centered around zero */
		/* org L1 */
		/* calc: prev - post */
		e = ABS(Mm->avail) + ABS(Mn->avail);
		tmp = Mm->avail + J->weight;
		e -= ABS(tmp);
		tmp = Mn->avail - J->weight;
		e -= ABS(tmp);
	} else { /* both are over or under */
		e = Mn->avail - J->weight - Mm->avail;
	}
	
	return e;
}

int handover(Machine *Mm, Machine *Mn) {
	
	int handovers;
	Job *J;
	Machine *middleman;
	
	if(Mn->avail < Mm->avail) {
		middleman = Mn;
		Mn = Mm;
		Mm = middleman;
	} else if(Mm->avail == Mn->avail) {
		return 0;
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
	while(J && 0 < testHandover(Mm, Mn, J)) {
		/* handover job */
		Mm->n--;
		Mn->n++;
		Mm->avail += J->weight;
		Mn->avail -= J->weight;
		Mm->jobs = J->next;
		J->next = 0;
		Mn->jobs = jobmerge_inc(Mn->jobs, J);
		++handovers;
		J = Mm->jobs;
	}
	
	return handovers;
}

int tradeBB(Machine *M) {
	
	/* Babettes Buckets */
	int trades, nullTrades;
	double min, test;
	Job *Jn, *Jm, *JnBest, *JmBest;
	Machine *Mm, *Mn, *Mbest;
	
	/* check MSE */
	test = M->m ? machineIMSE(M) : machineMSE(M);
	fprintf(stderr, "## Pre-tabu MSE:\t%f\n", test);
	if(test == 0) {
		return 0;
	}
	
	/* trade as long as there is an improvement */
	trades = 0;
	do {
		Mbest = 0;
		nullTrades = trades;
		
		/* go through machines */
		Mm = M;
		while(Mm) {
			/* Negotiate trade options and handovers for Mm */
			min = 0;
			Jm = 0;
			Jn = 0;
			JmBest = 0;
			JnBest = 0;
			Mn = Mm->next;
			while(Mn) {
				/* check handover */
				trades += handoverPtr(Mm, Mn);
				
				/* Negotiate trade options between Mm and Mn */
				test = negotiatePtr(Mm, Mn, &Jm, &Jn);
				
				/* check negotiation between machines */
				if(test < min) {
					min = test;
					JmBest = Jm;
					JnBest = Jn;
					Mbest = Mn;
				}
				
				Mn = Mn->next;
			}
			
			/* exchange jobs */
			if(min < 0 && exchangeJobs(Mm, Mbest, JmBest, JnBest)) {
				++trades;
			} else {
				Mm = Mm->next;
			}
		}
	} while(nullTrades != trades);
	
	return trades;
}

Job * tradeMsequential(Machine *M, Job *J, int m, int n) {
	
	int trades;
	long unsigned pos;
	Job **sJ;
	Machine *Mm, *Mn;
	
	/* turn jobs into arrays to allow binary search */
	sJ = sequenceJobs(M, J, m, n);
	
	/* trade as long as there is an improvement */
	do {
		trades = 0;
		
		/* go through machines */
		Mm = M;
		while(Mm) {
			/* Negotiate trade options for Mm */
			pos = 0;
			Mn = Mm->next;
			while(Mn) {
				/* Negotiate trade options between Mm and Mn */
				//test = exchange(Mm, Mn);
				pos = 1;
				Mn = Mn->next;
			}
			
			/* exchange jobs */
			if(pos) {
				jobexchange(sJ, pos >> 32, pos & 4294967295);
				++trades;
			}
			if(!pos) {
				Mm = Mm->next;
			}
		}
	} while(trades);
	
	return J;
}
