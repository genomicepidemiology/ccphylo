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

#include "jobs.h"
#include "machines.h"

extern void (*jobMVWeight)(Job *src, int m, int n, double logbase);

double addValue(Machine *M, Job *J);
void rmMVjob(Machine *M, Job *J);
void addMVjob(Machine *M, Job *J);
void addMVjobToMachine(Machine *M, Job *J);
void nullMVWeight(Job *src, int m, int n, double logbase);
void logMVWeight(Job *src, int m, int n, double logbase);
void polMVWeight(Job *src, int m, int n, double exponent);
void expMVWeight(Job *src, int m, int n, double expobase);
