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

#ifndef JOB
typedef struct job Job;
struct job {
	int num;
	int size;
	double weight;
	double *Weights;
	struct job *next;
};
#define JOB 1
#endif

extern void (*jobWeight)(Job *src, int n, double logbase);

Job * job_realloc(Job *src, int mv, int oldsize, int newsize);
Job * jobWeights_realloc(Job *src, int mv, int mv_new, int n);
Job * jobmerge(Job *L1, Job *L2);
Job * jobmerge_inc(Job *L1, Job *L2);
Job * jobsort(Job *src, int n);
double totM(Job *J, int n);
double * totMVM(Job *J, int n, int mv);
double optM(Job *J, int n, int m);
double targetM(Job *J, int n, int m);
void nullWeight(Job *src, int n, double logbase);
void logWeight(Job *src, int len, double logbase);
void polWeight(Job *src, int len, double exponent);
void expWeight(Job *src, int len, double expobase);
int cleanJobs(Job *src, int n);
int cmpJ(Job *Jm, Job *Jn, int m);
