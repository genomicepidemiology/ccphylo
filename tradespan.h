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

Job ** sequenceJobs(Machine *M, Job *J, int m, int n);
void moveUp(Job **sJ);
void moveDown(Job **sJ);
void exchange(Job **sJ, unsigned m, unsigned n);
void insertJob(Machine *M, Job *J);
void exchangeJobs(Machine *Mm, Machine *Mn, Job *Jm, Job *Jn);
double negotiateM(Machine *Mm, Machine *Mn, Job **Jmbest, Job **Jnbest);
int tradeM(Machine *M);
Job * tradeMsequential(Machine *M, Job *J, int m, int n);
