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
#include "mvjobs.h"
#include "mvtabusearch.h"
#define ABS(num)((num) < 0 ? -(num) : (num))

double baseValue(Machine *Mm, Machine *Mn) {
	
	int m;
	double base, *Mma, *Mna;
	
	/* init */
	base = 0;
	Mma = Mm->Avails;
	Mna = Mn->Avails;
	
	/* get error by trade */
	m = Mm->m + 1;
	while(--m) {
		if((*Mma < 0 && 0 < *Mna) || (*Mna < 0 && 0 < *Mma)) { /* machines are somewhat centered around zero */
			base += ABS(*Mma) + ABS(*Mna);
		} else if(*Mma < 0) { /* both are over */
			base -= *Mma < *Mna ? *Mma : *Mna;
		} else { /* both are under */
			base += *Mma < *Mna ? *Mna : *Mma;
		}
		++Mma;
		++Mna;
	}
	
	return base;
}

double optValue(Machine *Mm, Machine *Mn) {
	
	int m;
	double opt, diff, *Mma, *Mna;
	
	/* init */
	opt = 0;
	Mma = Mm->Avails - 1;
	Mna = Mn->Avails - 1;
	
	/* get error by trade */
	m = Mm->m + 1;
	while(--m) {
		diff = *++Mma + *++Mna;
		if((*Mma < 0 && 0 < *Mna) || (*Mna < 0 && 0 < *Mma)) { /* machines are somewhat centered around zero */
			opt += ABS(diff);
		} else { /* both are over or under*/
			opt += 0.5 * ABS(diff);
		}
	}
	
	return opt;
}

double tradeValue(Machine *Mm, Machine *Mn, Job *Jm, Job *Jn) {
	
	int m;
	double post, tm, tn, *Mma, *Mna, *Jmw, *Jnw;
	
	/* init */
	post = 0;
	Mma = Mm->Avails - 1;
	Mna = Mn->Avails - 1;
	Jmw = Jm->Weights - 1;
	Jnw = Jn->Weights - 1;
	
	/* get error by trade */
	m = Mm->m + 1;
	while(--m) {
		/* new error */
		tm = *++Mma + *++Jmw - *++Jnw;
		tn = *++Mna + *Jnw - *Jmw;
		if((*Mma < 0 && 0 < *Mna) || (*Mna < 0 && 0 < *Mma)) { /* machines are somewhat centered around zero */
			post += ABS(tm) + ABS(tn);
		} else { /* both are over or under*/
			tm = ABS(tm);
			tn = ABS(tn);
			post += tm < tn ? tn : tm;
		}
	}
	
	return post;
}

double negotiateMVM(Machine *Mm, Machine *Mn, Job **Jmbest, Job **Jnbest) {
	
	double min, test, best, base, opt;
	Job *Jm, *Jn, *next, *Jmin, *JmPrev, *JnPrev;
	
	/* exchange example:
	M1: eeddd
	M2: ccbbaaa
	 to
	M1: aaaddd
	M2: ccbbee
	*/
	
	/* no chance of improvement */
	if(Mm->n <= 1 && Mn->n <= 1) {
		return 0;
	}
	
	/* init minimum trade value */
	base = baseValue(Mm, Mn);
	opt = optValue(Mm, Mn);
	best = base;
	*Jmbest = 0;
	*Jnbest = 0;
	
	/* identify best trade option by minimizing multivariate target */
	Jm = Mm->jobs;
	JmPrev = 0;
	while(Jm) {
		/* set first option as best */
		Jn = Mn->jobs;
		JnPrev = 0;
		min = tradeValue(Mm, Mn, Jm, Jn);
		Jmin = JnPrev;
		
		/* check the rest */
		JnPrev = Jn;
		next = Jn->next;
		while(next) {
			test = tradeValue(Mm, Mn, Jm, next);
			if(test < min) {
				min = test;
				Jmin = JnPrev;
			}
			
			JnPrev = next;
			next = min == opt ? 0 : next->next;
		}
		
		/* test best */
		if(min < best) {
			best = min;
			/* set pointer to best exchange jobs */
			*Jmbest = JmPrev;
			*Jnbest = Jmin;
		}
		JmPrev = Jm;
		Jm = best <= opt ? 0 : Jm->next; /* early stopping */
	}
	
	/* limit double approximations */
	if(best != base) {
		best -= base;
	} else {
		best = 0;
	}
	
	/* return absolute reduction in error */
	return best;
}

double testMVhandover(Machine *Mm, Machine *Mn, Job *J) {
	
	int m;
	double post, prev, t1, t2, *mAvails, *nAvails, *jWeights;
	
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
		if((*mAvails < 0 && 0 < *nAvails) || (*nAvails < 0 && 0 < *mAvails)) { /* machines are somewhat centered around zero */
			prev += ABS(*mAvails) + ABS(*nAvails);
			post += ABS((*mAvails + *jWeights)) + ABS((*nAvails - *jWeights));
		} else if(*mAvails < 0) { /* both are over */
			prev -= *mAvails < *nAvails ? *mAvails : *nAvails;
			t1 = *mAvails + *jWeights;
			t1 = t1 < 0 ? t1 : -t1;
			t2 = *nAvails - *jWeights;
			post -= t1 < t2 ? t1 : t2;
		} else { /* both are under */
			prev += *mAvails < *nAvails ? *nAvails : *mAvails;
			t1 = ABS(*nAvails - *jWeights);
			t2 = *mAvails + *jWeights;
			post += t1 < t2 ? t2 : t1;
		}
		++mAvails;
		++nAvails;
		++jWeights;
	}
	
	return prev - post;
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
	while(J && Mm->avail + J->weight < Mn->avail - J->weight) {
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
			J = Mm->jobs;
		} else {
			J = J->next;
		}
	}
	
	return handovers;
}
