/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

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

#include "qseqs.h"

#ifndef MATPARSE
typedef struct nucCount NucCount;
typedef struct matrixCounts MatrixCounts;
struct nucCount {
	unsigned char *name;
	unsigned size;
	unsigned total;
	unsigned char ref;
	short unsigned counts[6];
};
struct matrixCounts {
	Qseqs *name;
	unsigned len;
	unsigned size;
	unsigned nNucs;
	unsigned char *refs;
	short unsigned *counts;
};
#define MATPARSE 1
#endif

NucCount * initNucCount(unsigned size);
void destroyNucCount(NucCount *src);
int FileBuffGetRow(FileBuff *src, NucCount *dest);
int FileBuffSkipTemplate(FileBuff *src, NucCount *dest);
MatrixCounts * initMat(unsigned matSize, unsigned nameSize);
void destroyMat(MatrixCounts *src);
void setMatName(MatrixCounts *dest, NucCount *src);
int FileBuffLoadMat(MatrixCounts *dest, FileBuff *src, unsigned minDepth);
