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

#ifndef MACHINES
typedef struct machine Machine;
struct machine {
	int num;
	int n;
	double avail;
	Job *jobs;
	struct machine *next;
};

#define MACHINES 1
#endif

Machine * machinemerge(Machine *L1, Machine *L2);
Machine * machinesort(Machine *src, int m);
Machine * initM(int m, int n, Job *J);
Machine * initSkewM(int m, int n, Job *J, double *loads);
double machineMSE(Machine *M);
void print_stats(Machine *M);