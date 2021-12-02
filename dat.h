/* Philip T.L.C. Clausen Jul 2021 plan@dtu.dk */

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

#ifndef DAT
typedef struct dat Dat;
struct dat {
	int m; /* # entries */
	int n; /* # used N */
	int M; /* # rows */
	int N; /* # columns */
	double **mat;
	float **fmat;
	short unsigned **smat;
	unsigned char **bmat;
	FILE *file;
};
#define DAT 1
#endif

extern Dat * (*Dat_init)(int, int);

Dat * DatInit(int M, int N);
Dat * DatMinit(int M, int N);
void Dat_mrealloc(Dat *src, int M);
void Dat_realloc(Dat *src, int M);
void Dat_mdestroy(Dat *src);
void Dat_destroy(Dat *src);
