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

extern Machine * (*makespan_method)(Machine *, Job *, int, int);
extern void (*addDBEptr)(Machine **, Machine **, Job *, int, int);
extern Machine * (*addDBFptr)(Machine *, Job *);
extern Machine *(*FirstFitptr)(Machine *, Job *, int);
extern Machine *(*FirstFetptr)(Machine *, Job *);

void addDBE(Machine **Mdest, Machine **Edest, Job *J, int m, int n);
Machine * DBE(Machine *M, Job *J, int m, int n);
Machine * addDBF(Machine *M, Job *J);
Machine * DBF(Machine *M, Job *J, int m, int n);
Machine * FirstFit(Machine *M, Job *J, int m);
Machine * DFF(Machine *M, Job *J, int m, int n);
Machine * FirstFet(Machine *M, Job *J);
Machine * DFE(Machine *M, Job *J, int m, int n);
void print_makespan(Machine *M, FILE *out, FILE *mout);
void makespan(char *inputfilename, char *outputfilename, char *moutputfilename, int m, double *loads, int mv, int *MV, double base, unsigned char sep, int col);
int main_makespan(int argc, char **argv);
