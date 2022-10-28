/* Philip T.L.C. Clausen Apr 2021 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
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

#include <stdio.h>
#include "matrix.h"
#include "qseqs.h"

int dbscan(Matrix *D, int *N, int *C, double maxDist, int minN);
void print_dbscan(Qseqs **names, int *N, int *C, int Dn, int nClust, double maxDist, int minN, FILE *out);
void make_dbscan(char *inputfilename, char *outputfilename, double maxDist, int minN, char sep, char quotes);
int main_dbscan(int argc, char **argv);
