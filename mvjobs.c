/* Philip T.L.C. Clausen Nov 2022 plan@dtu.dk */

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "jobs.h"
#include "machines.h"
#include "mvjobs.h"

void (*jobMVWeight)(Job *src, int m, int n, double logbase) = &nullMVWeight;

double addValue(Machine *M, Job *J) {
	
	int m;
	double e, *Ma, *Jw;
	
	/* init */
	e = 0;
	Ma = M->Avails;
	Jw = J->Weights;
	
	/* get absolute reduction in error */
	m = M->m + 1;
	while(--m) {
		if(*Jw <= *Ma) {
			e += *Jw;
		} else if(*Ma <= 0) {
			e -= *Jw;
		} else {
			e += *Ma + *Ma - *Jw; /* faster than 2 * *Ma - *jw */
		}
		++Ma;
		++Jw;
	}
	
	return e;
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

void addMVjobToMachine(Machine *M, Job *J) {
	
	M->n++;
	J->next = M->jobs;
	M->jobs = J;
	M->avail -= J->weight;
	addMVjob(M, J);
}

void nullMVWeight(Job *src, int m, int n, double logbase) {
	
	int i;
	double *weights;
	Job *dest;
	
	dest = src - 1;
	++n;
	while(--n) {
		(++dest)->weight = 0;
		
		i = m + 1;
		weights = dest->Weights - 1;
		while(--i) {
			dest->weight += *++weights;
		}
	}
}

void logMVWeight(Job *src, int m, int n, double logbase) {
	
	int i;
	double *weights;
	Job *dest;
	
	if(!logbase || !(logbase = log(logbase))) {
		fprintf(stderr, "Invalid logbase\n");
		exit(1);
	}
	
	dest = src - 1;
	++n;
	while(--n) {
		(++dest)->weight = 0;
		
		i = m + 1;
		weights = dest->Weights - 1;
		while(--i) {
			if(*++weights) {
				*weights = 1 + log(*weights) / logbase;
				dest->weight += *weights;
			}
		}
	}
}

void polMVWeight(Job *src, int m, int n, double exponent) {
	
	int i;
	double *weights;
	Job *dest;
	
	dest = src - 1;
	++n;
	while(--n) {
		(++dest)->weight = 0;
		
		i = m + 1;
		weights = dest->Weights - 1;
		while(--i) {
			if(*++weights) {
				*weights = pow(*weights, exponent);
				dest->weight += *weights;
			}
		}
	}
}

void expMVWeight(Job *src, int m, int n, double expobase) {
	
	int i;
	double *weights;
	Job *dest;
	
	dest = src - 1;
	++n;
	while(--n) {
		(++dest)->weight = 0;
		
		i = m + 1;
		weights = dest->Weights - 1;
		while(--i) {
			if(*++weights) {
				*weights = pow(expobase, *weights);
				dest->weight += *weights;
			}
		}
	}
}
