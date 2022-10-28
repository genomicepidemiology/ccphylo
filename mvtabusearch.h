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

double baseValue(Machine *Mm, Machine *Mn);
double tradeValue(Machine *Mm, Machine *Mn, Job *Jm, Job *Jn);
double negotiateMVM(Machine *Mm, Machine *Mn, Job **Jmbest, Job **Jnbest);
double testMVhandover(Machine *Mm, Machine *Mn, Job *J);
void rmMVjob(Machine *M, Job *J);
void addMVjob(Machine *M, Job *J);
int mvhandover(Machine *Mm, Machine *Mn);
